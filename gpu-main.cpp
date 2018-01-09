#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include "bitmap_image.hpp"
#include "sampler.h"

typedef unsigned char uchar;

const int CUBE_NULL  = 0;
const int CUBE_WATER = 1;
const int CUBE_GRASS = 2;
const int CUBE_MUD   = 3;
const int CUBE_LIGHT = 4;

int main() {
    std::string source;

    // Read Threshold.cl code
    {
        std::ifstream ifs ("raytracing.cl.cpp");
        getline (ifs, source, (char) ifs.eof());
    }

    try {
	// Get list of OpenCL platforms.
	std::vector<cl::Platform> platform;
	cl::Platform::get(&platform);

	if (platform.empty()) {
	    std::cerr << "OpenCL platforms not found." << std::endl;
	    return 1;
	}

	// Get first available GPU device which supports double precision.
	cl::Context context;
	std::vector<cl::Device> device;
	for(auto p = platform.begin(); device.empty() && p != platform.end(); p++) {
	    std::vector<cl::Device> pldev;

	    try {
		p->getDevices(CL_DEVICE_TYPE_GPU, &pldev);

		for(auto d = pldev.begin(); device.empty() && d != pldev.end(); d++) {
		    if (!d->getInfo<CL_DEVICE_AVAILABLE>()) continue;

		    std::string ext = d->getInfo<CL_DEVICE_EXTENSIONS>();

		    if (
			    ext.find("cl_khr_fp64") == std::string::npos &&
			    ext.find("cl_amd_fp64") == std::string::npos
		       ) continue;

		    device.push_back(*d);
		    context = cl::Context(device);
		}
	    } catch(...) {
		device.clear();
	    }
	}

	if (device.empty()) {
	    std::cerr << "GPUs with double precision not found." << std::endl;
	    return 1;
	}

	std::cout << device[0].getInfo<CL_DEVICE_NAME>() << std::endl;

	// Create command queue.
	cl::CommandQueue queue(context, device[0]);

	// Compile OpenCL program for found device.
	cl::Program program(context, cl::Program::Sources(
        1, std::make_pair(source.c_str(), source.size()))
    );

	try {
	    program.build(device);
	} catch (const cl::Error&) {
	    std::cerr
		<< "OpenCL compilation error" << std::endl
		<< program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device[0])
		<< std::endl;
	    return 1;
	}

    const int H = 512;
    const int W = 512;

    double t;
    
    char* pic[H][W];

    const int SAMPLE_COUNT = 2;
    const int RS = SAMPLE_COUNT*SAMPLE_COUNT;
    double sx[SAMPLE_COUNT*SAMPLE_COUNT], sy[SAMPLE_COUNT*SAMPLE_COUNT];
    samplesSquare(sx, sy, SAMPLE_COUNT);

	std::vector<uchar> pic_r(H*W), pic_g(H*W), pic_b(H*W);

    int MAX_X = 50, MAX_Y = 50, MAX_Z = 50;
    std::vector<int> Cubes(MAX_X*MAX_Y*MAX_Z);

#define CINDEX(x, y, z) ((x)*MAX_Y*MAX_Z + (y)*MAX_Z + (z))
    for(int lx = 0;lx < MAX_X;lx++)
        for(int lz = 0;lz < MAX_Z;lz++){
            Cubes[CINDEX(lx, 10, lz)] = CUBE_MUD;
            if(lx%2 == 1 and lz%2 == 1)
                Cubes[CINDEX(lx, 9, lz)] = CUBE_GRASS,
                Cubes[CINDEX(lx, 0, lz)] = CUBE_LIGHT;
        }
#undef CINDEX

	cl::Buffer R(context, CL_MEM_READ_WRITE, pic_r.size() * sizeof(uchar), pic_r.data());
	cl::Buffer G(context, CL_MEM_READ_WRITE, pic_g.size() * sizeof(uchar), pic_g.data());
	cl::Buffer B(context, CL_MEM_READ_WRITE, pic_b.size() * sizeof(uchar), pic_b.data());

	cl::Buffer SX(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, RS * sizeof(double), sx);
	cl::Buffer SY(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, RS * sizeof(double), sy);

	cl::Buffer CUBES(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Cubes.size() * sizeof(int), Cubes.data());

	cl::Kernel draw(program, "draw");

	draw.setArg(0, static_cast<cl_uint>(H));
	draw.setArg(1, static_cast<cl_uint>(W));
	draw.setArg(2, R);
	draw.setArg(3, G);
	draw.setArg(4, B);
	draw.setArg(5, CUBES);
	draw.setArg(6, SX);
	draw.setArg(7, SY);

	// Launch kernel on the compute device.
	queue.enqueueNDRangeKernel(draw, cl::NullRange, H*W, cl::NullRange);

    std::cout << "Run Ok" << std::endl;

	// Get result back to host.
	queue.enqueueReadBuffer(R, CL_TRUE, 0, pic_r.size() * sizeof(uchar), pic_r.data());
	queue.enqueueReadBuffer(G, CL_TRUE, 0, pic_g.size() * sizeof(uchar), pic_g.data());
	queue.enqueueReadBuffer(B, CL_TRUE, 0, pic_b.size() * sizeof(uchar), pic_b.data());

    bitmap_image image(W, H);

    for(int ly = 0;ly < H;ly++){
        for(int lx = 0;lx < W;lx++){
            rgb_t c;
            c.red = pic_r[ly*W + lx];
            c.green = pic_g[ly*W + lx];
            c.blue = pic_b[ly*W + lx];
            image.set_pixel(lx, ly, c);
        }
    }

    image.save_image("out.bmp");

    } catch (const cl::Error &err) {
	std::cerr
	    << "OpenCL error: "
	    << err.what() << "(" << err.err() << ")"
	    << std::endl;
	return 1;
    }
}
