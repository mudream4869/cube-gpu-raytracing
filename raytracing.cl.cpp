#if defined(cl_khr_fp64)
#  pragma OPENCL EXTENSION cl_khr_fp64: enable
#elif defined(cl_amd_fp64)
#  pragma OPENCL EXTENSION cl_amd_fp64: enable
#else
#  error double precision is not supported
#endif

__constant int MAX_X = 50, MAX_Y = 50, MAX_Z = 50;
__constant double cubeWidth = 50; // Cube size

__constant int CUBE_NULL  = 0;
__constant int CUBE_WATER = 1;
__constant int CUBE_GRASS = 2;
__constant int CUBE_MUD   = 3;
__constant int CUBE_LIGHT = 4;

struct _Ray{
    double3 o;
    double3 d;
    double t;
};

typedef struct _Ray Ray;

#define SAMPLE_COUNT 4

int getRand(int* prng){
    (*prng) = (*prng)*8763;
    (*prng) %= 1000000007;
    (*prng) ^= 0xdeadbeef;
    if((*prng) < 0) (*prng) = -(*prng);
    return *prng;
}

// Return Random in [0, 1)
double sample(int* prng){
    int a = getRand(prng)%12345678;
    return a/12345678.;
}

// |Sample| = n*n
void samplesSquare(int* prng, local double* x, local double* y, int n){
    double delta = 1./n;
    for(int lx = 0;lx < n;lx++){
        for(int ly = 0;ly < n;ly++){
            x[lx*n + ly] = delta*(lx + sample(prng));
            y[lx*n + ly] = delta*(ly + sample(prng));
        }
    }
    return;
}

int isZero(double x){
    return fabs(x) < 0.000001;
}

double3 callRay(Ray r, double tt){
    return r.o + r.d*tt;
}

int getCube(const global int* Cubes, uint3 i){
    return Cubes[i.x*MAX_Y*MAX_Z + i.y*MAX_Z + i.z];
}

// Return True if pos in [0, MAX_X)*[0, MAX_Y)*[0, MAX_Z)
int isOutOfWorld(double3 pos){
    if(pos.x < 0 || pos.y < 0 || pos.z < 0) return 1;
    if(pos.x >= MAX_X*cubeWidth ||
       pos.y >= MAX_Y*cubeWidth ||
       pos.z >= MAX_Z*cubeWidth){
        return 1;
    }
    return 0;
}

// Return True if Intersect a Cube
int intersectCubes(const global int* Cubes, Ray ray, uint3* index, double* now_t){
    double width = cubeWidth;
    (*now_t) = ray.t;
    for(;;){

        double3 curr = callRay(ray, (*now_t));

        if(isOutOfWorld(curr)){
            return 0;
        }

        int d_x = floor(curr.x/width),
            d_y = floor(curr.y/width),
            d_z = floor(curr.z/width);
        
        if(getCube(Cubes, (uint3)(d_x, d_y, d_z)) != CUBE_NULL){
            index->x = d_x;
            index->y = d_y;
            index->z = d_z;
            return 1;
        }
        
        double n_x = d_x*width, n_y = d_y*width, n_z = d_z*width;

        if(ray.d.x > 0) n_x += width;
        if(ray.d.y > 0) n_y += width;
        if(ray.d.z > 0) n_z += width;

        double min_t = 100000000;

        if(!isZero(ray.d.x)) min_t = min((n_x - ray.o.x)/ray.d.x, min_t);
        if(!isZero(ray.d.y)) min_t = min((n_y - ray.o.y)/ray.d.y, min_t);
        if(!isZero(ray.d.z)) min_t = min((n_z - ray.o.z)/ray.d.z, min_t);

        (*now_t) = min_t + 0.0001;
    }
    // Always return in loops
}

double3 getNormal(Ray ray, uint3 index, double t){
    double nx = (index.x + .5)*cubeWidth,
           ny = (index.y + .5)*cubeWidth,
           nz = (index.z + .5)*cubeWidth;

    double3 r = callRay(ray, t) - (double3)(nx, ny, nz);

    double ax = fabs(r.x),
           ay = fabs(r.y),
           az = fabs(r.z);

    if(ax >= max(ay, az)) return (double3)(r.x > 0 ? 1 : -1, 0, 0);
    if(ay >= max(ax, az)) return (double3)(0, r.y > 0 ? 1 : -1, 0);
    return (double3)(0, 0, r.z > 0 ? 1 : -1);
}

Ray reflect(int* prng, Ray ray, uint3 cubeIndex, double t){
    double3 I = ray.d;
    double3 N = getNormal(ray, cubeIndex, t);

    double sx = sample(prng), sy = sample(prng);
    sx -= 0.5, sy -= 0.5;
    sx *= 0.3, sy *= 0.3;
    if(isZero(fabs(N.x) - 1)) N.y = sx, N.z = sy;
    else if(isZero(fabs(N.y) - 1)) N.x = sx, N.z = sy;
    else if(isZero(fabs(N.z) - 1)) N.y = sx, N.x = sy;

    N = normalize(N);

    double dt = -dot(I, N);

    Ray ret = {callRay(ray, t), normalize(N*dt*2 + I), 0.01};
    return ret; 
}

double3 rayTracing(const global int* Cubes, int* prng, Ray ray, int level){
    const double3 white = (double3)(1, 1, 1);
    const double3 black = (double3)(0, 0, 0);
    const double3 red = (double3)(1, 0, 0);
    const double3 green = (double3)(0, 1, 0);
    const double3 brown = (double3)(139/255.,69/255.,19/255.);

    if(level >= 4){
        return black;
    }

    uint3 cubeIndex;
    double tt;

    if(intersectCubes(Cubes, ray, &cubeIndex, &tt)){
        int cube_id = getCube(Cubes, cubeIndex);
        
        if(cube_id == CUBE_LIGHT){
            return white;
        }

        double3 cube_color;
        if(cube_id == CUBE_MUD)   cube_color = brown;
        if(cube_id == CUBE_GRASS) cube_color = green;

        const int RAY_COUNT = 4;

        double3 ret = (double3)(0, 0, 0);
        for(int lx = 0;lx < RAY_COUNT;lx++){
            ret = ret + rayTracing(Cubes, prng, reflect(prng, ray, cubeIndex, tt), level+1);
        }

        ret = ret*cube_color;
        ret = ret/RAY_COUNT;
        return ret;
    }

    return black;
}


kernel void draw(
    uint H, uint W,
    global unsigned char *pic_r,
    global unsigned char *pic_g,
    global unsigned char *pic_b,
    const global int* Cubes
){
    size_t i = get_global_id(0);
    int prng = i;
    
    if(i >= H*W){
        return;
    }

    local double sx[SAMPLE_COUNT*SAMPLE_COUNT], sy[SAMPLE_COUNT*SAMPLE_COUNT];

    int x = i%W, y = i/W;

    double3 col = (double3)(0, 0, 0);

    samplesSquare(&prng, sx, sy, SAMPLE_COUNT);

    for(int lx = 0;lx < SAMPLE_COUNT*SAMPLE_COUNT;lx++){
        Ray ray = {(double3)(x+sx[lx] + 1300,y+sy[lx],0),
                   normalize((double3)(x+sx[lx] - W/2, y+sy[lx] - H/2, 512)),
                   0};
        col = col + rayTracing(Cubes, &prng, ray, 0);
    }

    col /= SAMPLE_COUNT*SAMPLE_COUNT;

    pic_r[y*W + x] = (int)(col.x*255);
    pic_g[y*W + x] = (int)(col.y*255);
    pic_b[y*W + x] = (int)(col.z*255);
    return;
}
