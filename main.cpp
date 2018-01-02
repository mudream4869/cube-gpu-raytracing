#include <cmath>
#include <cassert>

#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <utility>
#include <vector>

#include "basic.h"
#include "sampler.h"

struct Sphere {
    Vec3d c;
    double r;
    Sphere(const Vec3d& c, double r) : c(c), r(r) {}
    Vec3d getNormal(const Vec3d& pi) const { return (pi - c) / r; }
    bool intersect(const Ray& ray, double &t) const {
        const Vec3d o = ray.o;
        const Vec3d d = ray.d;
        const Vec3d oc = o - c;
        const double b = 2 * dot(oc, d);
        const double c = dot(oc, oc) - r*r;
        double disc = b*b - 4 * c;
        if (disc < 1e-4) return false;
        disc = sqrt(disc);
        const double t0 = -b - disc;
        const double t1 = -b + disc;
        t = (t0 < t1) ? t0 : t1;
        return true;
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

Vec3d getNormal(const Ray& ray, Vec3i index, double t){
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

Ray reflect(const Ray& ray, Vec3i cubeIndex, double t){
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

Vec3d rayTracing(const Ray& ray, int level = 0){
    const Vec3d white(1, 1, 1);
    const Vec3d black(0, 0, 0);
    const Vec3d red(1, 0, 0);

    if(level >= 10){
        return black;
    }

    const Sphere light(Vec3d(500, 0, 500), 100);

    Vec3i cubeIndex;
    double tt;

    if(light.intersect(ray, tt)){
        return white;
    }

    if(intersectCubes(ray, cubeIndex, tt)){
        return rayTracing(reflect(ray, cubeIndex, tt), level+1)*red*0.95;
    }

    return black;
}

int main(){
    const int H = 500;
    const int W = 500;

    for(int lx = 0;lx < MAX_X;lx++)
        for(int lz = 0;lz < MAX_Z;lz++){
            Cubes[lx][7][lz] = 1;
            if(lx&1 and lz&1)
                Cubes[lx][6][lz] = 1;
        }

    std::ofstream out("out.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";

    double t;
    
    Vec3d eye(W/2, H/2, -500);

    for (int y = 0; y < H; ++y){
        for (int x = 0; x < W; ++x){
            std::cout << "\r" << y*W + x << "/" << W*H;

            Vec3d col;

            for(int s = 0;s < 64;s++){ 
                double sx, sy;
                sampleSquare(sx, sy);
                Ray ray(Vec3d(x+sx,y+sy,0), Vec3d(x+sx, y+sy, 0) - eye);
                col = col + rayTracing(ray);
            }

            col = col/32.;
            col = col*255;

            out << (int)col.x << ' '
                << (int)col.y << ' '
                << (int)col.z << '\n';
        }
    }
    std::cout << std::endl;
    return 0;
}
