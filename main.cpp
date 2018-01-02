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

        //std::cout << curr << std::endl;
        //std::cout << ray.d << std::endl;

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

Ray reflect(const Ray& ray, const Vec3i& cubeIndex, double t){
    Vec3d I = ray.d;
    Vec3d N = getNormal(ray, cubeIndex, t);

    double sx, sy;
    sampleSquare(sx, sy);

    if(isZero(std::abs(N.x) - 1)) N.y = sx, N.z = sy;
    else if(isZero(std::abs(N.y) - 1)) N.x = sx, N.z = sy;
    else if(isZero(std::abs(N.z) - 1)) N.y = sx, N.x = sy;

    double dt = dot(I.normalize()*(-1), N.normalize());

    return Ray(ray(t), N*dt*2 + ray.d, 0.01);
}

Ray patrenate(const Ray& ray, const Vec3i& cubeIndex, double t){
    Vec3d curr = ray(t);

    int d_x = floor(curr.x/cubeWidth),
        d_y = floor(curr.y/cubeWidth),
        d_z = floor(curr.z/cubeWidth);

    double n_x = d_x*cubeWidth, n_y = d_y*cubeWidth, n_z = d_z*cubeWidth;

    if(ray.d.x > 0) n_x += cubeWidth;
    if(ray.d.y > 0) n_y += cubeWidth;
    if(ray.d.z > 0) n_z += cubeWidth;

    double min_t = std::numeric_limits<double>::infinity();

    if(not isZero(ray.d.x)) min_t = std::min((n_x - ray.o.x)/ray.d.x, min_t);
    if(not isZero(ray.d.y)) min_t = std::min((n_y - ray.o.y)/ray.d.y, min_t);
    if(not isZero(ray.d.z)) min_t = std::min((n_z - ray.o.z)/ray.d.z, min_t);

    return Ray(ray.o, ray.d, min_t + 0.01); 
}

Vec3d rayTracing(const Ray& ray, int level = 0){
    const Vec3d white(1, 1, 1);
    const Vec3d black(0, 0, 0);
    const Vec3d red(1, 0, 0);
    const Vec3d green(0, 1, 0);
    const Vec3d brown(139/255.,69/255.,19/255.);

    if(level >= 10){
        return black;
    }

    Vec3i cubeIndex;
    double tt;

    if(intersectCubes(ray, cubeIndex, tt)){
        int cube_id = getCube(cubeIndex);
        
        if(cube_id == CUBE_LIGHT){
            if(level)
                return white;
            else
                return rayTracing(patrenate(ray, cubeIndex, tt), level);
        }

        Vec3d cube_color;
        if(cube_id == CUBE_MUD)   cube_color = brown;
        if(cube_id == CUBE_GRASS) cube_color = green;
        return rayTracing(reflect(ray, cubeIndex, tt), level+1)*cube_color;
    }

    return black;
}

int main(){
    const int H = 512;
    const int W = 512;

    for(int lx = 0;lx < MAX_X;lx++)
        for(int lz = 0;lz < MAX_Z;lz++){
            Cubes[lx][9][lz] = CUBE_MUD;
            if(lx%2 == 1 and lz%2 == 1)
                Cubes[lx][8][lz] = CUBE_GRASS;
            if(lx%4 == 0 and lz%4 == 0)
                Cubes[lx][8][lz] = CUBE_LIGHT;
        }

    double t;
    
    Vec3d eye(W/2, H/2, -512);
    Vec3d pic[H][W];

    auto sample_function = [&pic, eye](int id){
        for(int y = 0; y < H; ++y){
            std::cout << y*100/H << "%\r";
            std::cout << std::flush;

            for(int _x = 0; _x < W/4; ++_x){

                int x = _x*4 + id;

                const int SAMPLE_COUNT = 128;

                double sx[SAMPLE_COUNT], sy[SAMPLE_COUNT];
                samplesSquare(sx, sy, SAMPLE_COUNT);

                Vec3d col;

                for(int lx = 0;lx < SAMPLE_COUNT;lx++){
                    int cid = lx;
                    Ray ray(Vec3d(x+sx[cid] + 1300,y+sy[cid],0), Vec3d(x+sx[cid], y+sy[cid], 0) - eye);
                    col = col + rayTracing(ray);
                }

                pic[y][x] = col/SAMPLE_COUNT;
            }
        }
    };
    
    std::thread ths[4];

    for(int lx = 0;lx < 4;lx++){
        ths[lx] = std::thread(sample_function, lx); 
    }

    for(int lx = 0;lx < 4;lx++){
        ths[lx].join();
    }

    std::ofstream out("out.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";

    for(int ly = 0;ly < H;ly++){
        for(int lx = 0;lx < W;lx++){
            auto& col = pic[ly][lx];
            col = col*255;
            out << (int)col.x << ' '
                << (int)col.y << ' '
                << (int)col.z << '\n';
        }
    }

    out.close();

    std::cout << std::endl;
    return 0;
}
