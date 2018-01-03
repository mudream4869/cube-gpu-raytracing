render.out: main.cpp
	g++ main.cpp -std=c++11 -O2 -o render.out

gpu-render.out: gpu-main.cpp
	g++ gpu-main.cpp -I/usr/local/cuda/include/ -lOpenCL -std=c++11 -o gpu-render.out
