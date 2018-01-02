#include <cmath>
#include <cassert>

#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <utility>
#include <vector>

double clamp_(double x, double l, double u){
    return (x < l) ? l : ((x > u) ? u : x);
}

template<typename T>
struct Vec3 {
    T x,y,z;
    Vec3<T>(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}
    Vec3<T> operator + (const Vec3<T>& v) const { return Vec3<T>(x+v.x, y+v.y, z+v.z); }
    Vec3<T> operator - (const Vec3<T>& v) const { return Vec3<T>(x-v.x, y-v.y, z-v.z); }
    Vec3<T> operator * (T d) const { return Vec3<T>(x*d, y*d, z*d); }
    Vec3<T> operator / (T d) const { return Vec3<T>(x/d, y/d, z/d); }

    bool operator < (const Vec3<T>& v) const {
        return std::tie(x, y, z) < std::tie(v.x, v.y, v.z);
    }

    double len() const {
        return sqrt(x*x + y*y + z*z);
    }

    Vec3<T> normalize() const {
        double mg = sqrt(x*x + y*y + z*z);
        return Vec3<T>(x/mg,y/mg,z/mg);
    }

    void clamp(double lower_bound, double upper_bound) {
        x = clamp_(x, lower_bound, upper_bound); 
        y = clamp_(y, lower_bound, upper_bound); 
        z = clamp_(z, lower_bound, upper_bound); 
        return;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const Vec3<T>& v){
    o << "(" << v.x << "," << v.y << "," << v.z << ")";
    return o;
}

template<typename T>
T dot(const Vec3<T>& a, const Vec3<T>& b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;

struct Ray {
    Vec3d o,d;
    double t;
    Ray(const Vec3d& o, const Vec3d& d, double t = 0) : o(o), d(d.normalize()), t(t) {}

    Vec3d operator()(double _t) const{
        return o + d*_t;
    }
};


const int MAX_X = 50, MAX_Y = 50, MAX_Z = 50;
const double cubeWidth = 50; // Cube size
bool Cubes[MAX_X][MAX_Y][MAX_Z] = {0};

// Return True if pos in [0, MAX_X)*[0, MAX_Y)*[0, MAX_Z)
bool isOutOfWorld(Vec3d pos){
    if(pos.x < 0 or pos.y < 0 or pos.z < 0) return true;
    if(pos.x >= MAX_X*cubeWidth or
       pos.y >= MAX_Y*cubeWidth or
       pos.z >= MAX_Z*cubeWidth){
        return true;
    }
    return false;
}

// Return True if |x| < eps
bool isZero(double x){
    return std::abs(x) < std::numeric_limits<double>::min();
}

// Return True if Intersect a Cube
bool intersectCubes(Ray& ray, Vec3i& index, double& now_t){
    double width = cubeWidth;
    now_t = ray.t;
    for(;;){

        Vec3d curr = ray(now_t);

        //std::cout << curr << std::endl;
        //std::cout << ray.d << std::endl;

        if(isOutOfWorld(curr)){
            return false;
        }

        int d_x = floor(curr.x/width),
            d_y = floor(curr.y/width),
            d_z = floor(curr.z/width);
        
        if(Cubes[d_x][d_y][d_z]){
            index = Vec3i(d_x, d_y, d_z);
            return true;
        }
        
        double n_x = d_x*width, n_y = d_y*width, n_z = d_z*width;

        if(ray.d.x > 0) n_x += width;
        if(ray.d.y > 0) n_y += width;
        if(ray.d.z > 0) n_z += width;

        double min_t = std::numeric_limits<double>::infinity();

        if(not isZero(ray.d.x)) min_t = std::min((n_x - ray.o.x)/ray.d.x, min_t);
        if(not isZero(ray.d.y)) min_t = std::min((n_y - ray.o.y)/ray.d.y, min_t);
        if(not isZero(ray.d.z)) min_t = std::min((n_z - ray.o.z)/ray.d.z, min_t);

        now_t = min_t + 0.0001; 
    }
    // Always return in loops
}

Vec3d getNormal(Ray& ray, Vec3i index, double t){
    std::vector< std::pair<double, Vec3d> > normal_list;
    
    double nx = index.x*cubeWidth,
           ny = index.y*cubeWidth,
           nz = index.z*cubeWidth;
           
    if(not isZero(ray.d.x)){
        normal_list.push_back({(nx - ray.o.x)/ray.d.x, Vec3d(-1, 0, 0)});
        normal_list.push_back({(nx+cubeWidth - ray.o.x)/ray.d.x, Vec3d(1, 0, 0)});
    }

    if(not isZero(ray.d.y)){
        normal_list.push_back({(ny - ray.o.y)/ray.d.y, Vec3d(0, -1, 0)});
        normal_list.push_back({(ny+cubeWidth - ray.o.y)/ray.d.y, Vec3d(0, 1, 0)});
    }

    if(not isZero(ray.d.z)){
        normal_list.push_back({(nz - ray.o.z)/ray.d.z, Vec3d(0, 0, -1)});
        normal_list.push_back({(nz+cubeWidth - ray.o.z)/ray.d.z, Vec3d(0, 0, 1)});
    }


    for(auto& p : normal_list){
        p.first = std::abs(p.first - t);
    }

    std::sort(normal_list.begin(), normal_list.end());

    assert(normal_list.size());
    return normal_list[0].second;
}

int main(){

    const int H = 500;
    const int W = 500;

    const Vec3d white(255, 255, 255);
    const Vec3d black(0, 0, 0);
    const Vec3d red(255, 0, 0);

    Cubes[8][8][1] = 1;

    const Vec3d light(500, 0, 0);

    std::ofstream out("out.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";

    double t;
    
    Vec3d eye(W/2, H/2, -500);

    for (int y = 0; y < H; ++y){
        for (int x = 0; x < W; ++x){
            Vec3d pix_col(black);
            Ray ray(Vec3d(x,y,0), Vec3d(x, y, 0) - eye);

            Vec3i cubeIndex;
            double tt;
            if(intersectCubes(ray, cubeIndex, tt)){
                Vec3d L = light - ray(tt);
                Vec3d N = getNormal(ray, cubeIndex, tt); 
                double dt = dot(L.normalize(), N.normalize());
                double dis = L.len();

                pix_col = (red + white*dt);
                pix_col.clamp(0, 255);
            }
            out << (int)pix_col.x << ' '
                << (int)pix_col.y << ' '
                << (int)pix_col.z << '\n';
        }
    }
}
