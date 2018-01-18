#include <cmath>
#include <cassert>

#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <utility>
#include <vector>
#include <thread>
#include <functional>

#include "basic.h"
#include "sampler.h"

#include "bitmap_image.hpp"

const int MAX_X = 50, MAX_Y = 50, MAX_Z = 50;
const double cubeWidth = 50; // Cube size
int Cubes[MAX_X][MAX_Y][MAX_Z] = {0};

const int CUBE_NULL  = 0;
const int CUBE_WATER = 1;
const int CUBE_GRASS = 2;
const int CUBE_MUD   = 3;
const int CUBE_LIGHT = 4;

int& getCube(const Vec3i& i){
    return Cubes[i.x][i.y][i.z];
}

// Return True if pos in [0, MAX_X)*[0, MAX_Y)*[0, MAX_Z)
bool isOutOfWorld(const Vec3d& pos){
    if(pos.x < 0 or pos.y < 0 or pos.z < 0) return true;
    if(pos.x >= MAX_X*cubeWidth or
       pos.y >= MAX_Y*cubeWidth or
       pos.z >= MAX_Z*cubeWidth){
        return true;
    }
    return false;
}

// Return True if Intersect a Cube
bool intersectCubes(const Ray& ray, Vec3i& index, double& now_t){
    double width = cubeWidth;
    now_t = ray.t;
    for(;;){

        Vec3d curr = ray(now_t);

        if(isOutOfWorld(curr)){
            return false;
        }

        int d_x = floor(curr.x/width),
            d_y = floor(curr.y/width),
            d_z = floor(curr.z/width);
        
        if(Cubes[d_x][d_y][d_z] != CUBE_NULL){
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

Vec3d getNormal(const Ray& ray, const Vec3i& index, double t){
    double nx = (index.x + .5)*cubeWidth,
           ny = (index.y + .5)*cubeWidth,
           nz = (index.z + .5)*cubeWidth;

    Vec3d r = ray(t) - Vec3d(nx, ny, nz);

    double ax = std::abs(r.x),
           ay = std::abs(r.y),
           az = std::abs(r.z);

    if(ax >= std::max(ay, az)) return Vec3d(r.x > 0 ? 1 : -1, 0, 0);
    if(ay >= std::max(ax, az)) return Vec3d(0, r.y > 0 ? 1 : -1, 0);
    if(az >= std::max(ax, ay)) return Vec3d(0, 0, r.z > 0 ? 1 : -1);
    
    assert(0);
}

Ray reflect(int* prng, const Ray& ray, const Vec3i& cubeIndex, double t){
    Vec3d I = ray.d;
    Vec3d N = getNormal(ray, cubeIndex, t);
    
    // Scattering
    double sx = sample(prng), sy = sample(prng);
    sx -= 0.5, sy -= 0.5;
    sx *= 0.3, sy *= 0.3;

    if(isZero(std::abs(N.x) - 1)) N.y = sx, N.z = sy;
    else if(isZero(std::abs(N.y) - 1)) N.x = sx, N.z = sy;
    else if(isZero(std::abs(N.z) - 1)) N.y = sx, N.x = sy;

    N = N.normalize();

    double dt = -dot(I, N);

    return Ray(ray(t), N*dt*2 + I, 0.01);
}

Vec3d rayTracing(int* prng, const Ray& ray, int level = 0){
    const Vec3d white(1, 1, 1);
    const Vec3d black(0, 0, 0);
    const Vec3d red(1, 0, 0);
    const Vec3d green(0, 1, 0);
    const Vec3d brown(139/255.,69/255.,19/255.);

    if(level >= 4){
        return black;
    }

    Vec3i cubeIndex;
    double tt;

    // Only Intersect with Cubes
    if(intersectCubes(ray, cubeIndex, tt)){
        int cube_id = getCube(cubeIndex);
        
        if(cube_id == CUBE_LIGHT){
            return white;
        }

        Vec3d cube_color;
        if(cube_id == CUBE_MUD)   cube_color = brown;
        if(cube_id == CUBE_GRASS) cube_color = green;

        const int RAY_COUNT = 3;

        Vec3d ret;
        for(int lx = 0;lx < RAY_COUNT;lx++){
            ret = ret + rayTracing(prng, reflect(prng, ray, cubeIndex, tt), level+1);
        }

        ret = ret*cube_color;
        ret = ret/RAY_COUNT;
        return ret;
    }

    return black;
}

int main(int argc, char** argv){

    const int H = 512;
    const int W = 512;

    for(int lx = 0;lx < MAX_X;lx++)
        for(int lz = 0;lz < MAX_Z;lz++){
            Cubes[lx][10][lz] = CUBE_MUD;
            if(lx%2 == 1 and lz%2 == 1)
                Cubes[lx][9][lz] = CUBE_GRASS,
                Cubes[lx][0][lz] = CUBE_LIGHT;
        }
    
    double t;
    
    Vec3d eye(W/2, H/2, -512);
    Vec3d pic[H][W];

    int SAMPLE_COUNT = 16;
    int THREAD_COUNT = 1;

    assert(W%THREAD_COUNT == 0);

    auto sample_function = [&pic, eye, SAMPLE_COUNT, THREAD_COUNT](int id){
        int prng = id;
        for(int y = 0; y < H; ++y){
            std::cout << y*100/H << "%\r";
            std::cout << std::flush;

            for(int _x = 0; _x < W/THREAD_COUNT; ++_x){

                int x = _x*THREAD_COUNT + id;

                double sx[SAMPLE_COUNT*SAMPLE_COUNT], sy[SAMPLE_COUNT*SAMPLE_COUNT];
                samplesSquare(&prng, sx, sy, SAMPLE_COUNT);

                Vec3d col;

                for(int lx = 0;lx < SAMPLE_COUNT*SAMPLE_COUNT;lx++){
                    int cid = lx;
                    Ray ray(Vec3d(x+sx[cid] + 1300,y+sy[cid],0), Vec3d(x+sx[cid], y+sy[cid], 0) - eye);
                    col = col + rayTracing(&prng, ray);
                }

                pic[y][x] = col/SAMPLE_COUNT/SAMPLE_COUNT;
            }
        }
    };
    
    std::vector<std::thread> ths(THREAD_COUNT);

    for(int lx = 0;lx < THREAD_COUNT;lx++){
        ths[lx] = std::thread(sample_function, lx); 
    }

    for(int lx = 0;lx < THREAD_COUNT;lx++){
        ths[lx].join();
    }

    bitmap_image image(W, H);

    for(int ly = 0;ly < H;ly++){
        for(int lx = 0;lx < W;lx++){
            auto& col = pic[ly][lx];
            col = col*255;
            rgb_t c;
            c.red = col.x; 
            c.green = col.y; 
            c.blue = col.z;
 
            image.set_pixel(lx, ly, c);
        }
    }

    image.save_image("out.bmp");

    std::cout << std::endl;
    return 0;
}
