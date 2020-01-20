/*
 ============================================================================
 Name        : CUDANewton.cu
 Author      : Jakob Vokac
 Version     :
 Copyright   : Your copyright notice
 Description : Bezier Surface Distance algorithm implemented in CUDA

 This is the CUDA translation of the Bezier Surface Distance algorithm on:
 https://github.com/JakobVokac/BezierSurfaceDistance

 For more details on the algorithm, vision the repo above. 

 The main differences between this implementation and the C++ version are:

	- this implementation only includes the Newton-Raphson method
	  with bisection and quadratic interpolation for preprocessing,
	  which has been tested to perform the most reliably and the fastest on a CPU,

 	- interfaces have been taken out due to limited CUDA support for C++ and
	  the implementations and optimizer classes have been altered to fit the
	  new implementation accordingly,

	- all classes can now function as host or device code, so the optimizer
	  can be tested on both the CPU and GPU from this code,

	- the main function invokes a kernel function for the GPU to compute the input;
	  currently all surfaces are passed to each block and placed into shared memory,
	  which has proven to be suboptimal, this should be fixed in future work

	- the measurement and plotting classes have been left out as they only serve to
	  test the reliabiliy of the algorithm, which should perform identically in this
	  version and the original C++ version

 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include "model.h"
#include "geometry/vector/vec3d.h"
#include "geometry/curve/cubiccrv.h"
#include "geometry/surface/TopParametric.h"
#include "optimizer/optimizer.h"
static void CheckCudaErrorAux (const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

using namespace std;

//This function runs in parallel for each thread. The function takes in the split collections of points and
//gives them to the allocated blocks. Each block takes care of one collection and outputs the distances to global memory.
__global__ void kernel(
		vec3d *points1,
		vec3d *points2,
		vec3d *points3,
		vec3d *points4,
		vec3d *points5,
		vec3d *points6,
		vec3d *points7,
		vec3d *points8,
		vec3d *points9,
		vec3d *points10,
		vec3d *points11,
		vec3d *points12,
		int idx1,
		int idx2,
		int idx3,
		int idx4,
		int idx5,
		int idx6,
		int idx7,
		int idx8,
		int idx9,
		int idx10,
		int idx11,
		int idx12,
		double *distances, int size,
		TopParametric *sur1,
		TopParametric *sur2,
		TopParametric *sur3,
		TopParametric *sur4,
		TopParametric *sur5,
		TopParametric *sur6,
		BottomParametric *sur7,
		BottomParametric *sur8,
		BottomParametric *sur9,
		BottomParametric *sur10,
		BottomParametric *sur11,
		BottomParametric *sur12
	){

	optimizer op = optimizer(
				bisection(8),
				quadraticInterpolation(8),
				Newton(1.0),
				sur1,
				sur7,
				vec3d(0,0,0),
				0.00000001,
				20,
				2
	);

	if(blockIdx.x < 6){
		TopParametric top;
		vec3d *points;
		int size;
		if(blockIdx.x == 0){
			top = *sur1;
			points = points1;
			size = idx1;
		}
		if(blockIdx.x == 1){
			top = *sur2;
			points = points2;
			size = idx2;
		}
		if(blockIdx.x == 2){
			top = *sur3;
			points = points3;
			size = idx3;
		}
		if(blockIdx.x == 3){
			top = *sur4;
			points = points4;
			size = idx4;
		}
		if(blockIdx.x == 4){
			top = *sur5;
			points = points5;
			size = idx5;
		}
		if(blockIdx.x == 5){
			top = *sur6;
			points = points6;
			size = idx6;
		}

		op.setTopOrBot(1);
		op.setTop(&top);

		for(int i = 0; i < size; i+= blockDim.x){
			int arrId = threadIdx.x + i;

			if(arrId >= size)
				break;
			vec3d P = points[arrId];

			OptState2D loc = op.optimizeForPoint(P);

			distances[arrId] = loc.dist;
		}

	}else{
		BottomParametric bot;
		vec3d *points;
		int size;
		if(blockIdx.x == 6){
			bot = *sur7;
			points = points7;
			size = idx7;
		}
		if(blockIdx.x == 7){
			bot = *sur8;
			points = points8;
			size = idx8;
		}
		if(blockIdx.x == 8){
			bot = *sur9;
			points = points9;
			size = idx9;
		}
		if(blockIdx.x == 9){
			bot = *sur10;
			points = points10;
			size = idx10;
		}
		if(blockIdx.x == 10){
			bot = *sur11;
			points = points11;
			size = idx11;
		}
		if(blockIdx.x == 11){
			bot = *sur12;
			points = points12;
			size = idx12;
		}

		op.setTopOrBot(0);
		op.setBot(&bot);

		for(int i = 0; i < size; i+= blockDim.x){
			int arrId = threadIdx.x + i;

			if(arrId >= size)
				break;
			vec3d P = points[arrId];

			OptState2D loc = op.optimizeForPoint(P);

			distances[arrId] = loc.dist;
		}

	}

}


//The main function takes in the input from the input.txt file and creates the model.
//The specifications of the current model are passed as parameters to the model class.
//It then gets all of the surfaces and builds the optimizer, then passes both to the GPU.
//It also splits the collection of points into 12 arrays, one for each block.
//It then calls the kernel function.
int main(void)
{
	std::vector<double> inputPoints;
	ifstream inputFile("input.txt");        // Input file stream object

	// Check if exists and then open the file.
	if (inputFile.good()) {
		// Push items into a vector
		double current_number = 0;
		while (inputFile >> current_number){
			inputPoints.push_back(current_number);
		}

		// Close the file.
		inputFile.close();

		cout << endl;
	}else {
		cout << "Error reading input file!" << endl;
		cout << "Input file must have name \"input.txt\" and must consist of only numbers, no text!" << endl;
		cout << "The program will read the numbers in triples and each 3 numbers will be interpreted as one point." << endl;

		exit(0);
	}

	if(inputPoints.size() % 3 != 0){
		cout << "Number of input numbers should be divisible by 3 (3 dimensions per point)!" << endl;
		cout << "Ignoring last few numbers." << endl;

		int size = inputPoints.size();
		int truncate = size % 3;
		for (int i = 0; i < truncate; ++i) {
			inputPoints.pop_back();
		}
	}


	int N = inputPoints.size()/3;
	vec3d *points1, *points2, *points3, *points4, *points5, *points6, *points7, *points8, *points9, *points10, *points11, *points12;
	int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9, idx10, idx11, idx12;
	idx1 = idx2 = idx3 = idx4 = idx5 = idx6 = idx7 = idx8 = idx9 = idx10 = idx11 = idx12 = 0;
	double *distances;
	size_t freeMem, totalMem;
	TopParametric *sur1, *sur2, *sur3, *sur4, *sur5, *sur6;
	BottomParametric *sur7, *sur8, *sur9, *sur10, *sur11, *sur12;

	Model model = Model(
			12.0,
            59.5/180 * M_PI,
            11.4,
            14.4,
            10,
            1.2,
            50.0/180*M_PI,
            7.2,
            16.8,
            3.5,
            1.35,
            -0.2,
            -0.2,
            0.01,
            1.0,
            6.5);

    Model p0 = Model::getPart(model,0),
          p1 = Model::getPart(model,1),
          p2 = Model::getPart(model,2),
          p3 = Model::getPart(model,3),
          p4 = Model::getPart(model,4),
          p5 = Model::getPart(model,5);

	cudaMemGetInfo(&freeMem,&totalMem);
	printf("free memory before data load: %d, total memory: %d\n",freeMem,totalMem);

	// Allocate Unified Memory – accessible from CPU or GPU
	CUDA_CHECK_RETURN(cudaMallocManaged(&points1, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points2, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points3, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points4, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points5, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points6, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points7, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points8, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points9, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points10, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points11, N*sizeof(vec3d)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&points12, N*sizeof(vec3d)));


	CUDA_CHECK_RETURN(cudaMallocManaged(&distances, N*sizeof(double)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur1, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur2, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur3, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur4, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur5, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur6, sizeof(TopParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur7,  sizeof(BottomParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur8,  sizeof(BottomParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur9,  sizeof(BottomParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur10, sizeof(BottomParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur11, sizeof(BottomParametric)));
	CUDA_CHECK_RETURN(cudaMallocManaged(&sur12, sizeof(BottomParametric)));

	cudaMemGetInfo(&freeMem,&totalMem);
	printf("free memory after data load: %d, total memory: %d\n",freeMem,totalMem);

	*sur1 = p0.getTopParametric();
	*sur2 = p3.getTopParametric();
	*sur3 = p2.getTopParametric();
	*sur4 = p5.getTopParametric();
	*sur5 = p4.getTopParametric();
	*sur6 = p1.getTopParametric();

	*sur7 = p0.getBottomParametric();
	*sur8 = p3.getBottomParametric();
	*sur9 = p2.getBottomParametric();
	*sur10 = p5.getBottomParametric();
	*sur11 = p4.getBottomParametric();
	*sur12 = p1.getBottomParametric();

	srand(1);
	vec3d h = sur1->at(0,0);
	double divHeight = h.z;
	// initialize x and y arrays on the host
	for (int i = 0; i < N; i++) {

        vec3d P = {double(rand())/RAND_MAX*20 + 5,double(rand())/RAND_MAX*20 + 5,double(rand())/RAND_MAX*20 + 5};

		if(P.x == 0){
			if(P.y == 0){
				if(P.z > divHeight){
					points1[idx1] = P;
					idx1++;
				}else{
					points7[idx7] = P;
					idx7++;
				}
			}
			else if(P.y > 0){
				if(P.z > divHeight){
					points2[idx2] = P;
					idx2++;
				}else{
					points8[idx8] = P;
					idx8++;
				}
			}
			else if(P.y < 0){
				if(P.z > divHeight){
					points5[idx5] = P;
					idx5++;
				}else{
					points11[idx11] = P;
					idx11++;
				}
			}
		}

		if(P.x > 0 && P.y > 0){
			if(P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
				if(P.z > divHeight){
					points1[idx1] = P;
					idx1++;
				}else{
					points7[idx7] = P;
					idx7++;
				}
			}else if(P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				if(P.z > divHeight){
					points1[idx1] = P;
					idx1++;
				}else{
					points7[idx7] = P;
					idx7++;
				}
			}else{
				if(P.z > divHeight){
					points1[idx2] = P;
					idx2++;
				}else{
					points8[idx8] = P;
					idx8++;
				}
			}
		}else if(P.x < 0 && P.y < 0){
			if(P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
				if(P.z > divHeight){
					points3[idx3] = P;
					idx3++;
				}else{
					points9[idx9] = P;
					idx9++;
				}
			}else if(P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				if(P.z > divHeight){
					points3[idx3] = P;
					idx3++;
				}else{
					points9[idx9] = P;
					idx9++;
				}
			}else{
				if(P.z > divHeight){
					points1[idx2] = P;
					idx2++;
				}else{
					points8[idx8] = P;
					idx8++;
				}
			}
		}else if(P.x < 0 && P.y < 0){
			if(-P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
				if(P.z > divHeight){
					points4[idx4] = P;
					idx4++;
				}else{
					points10[idx10] = P;
					idx10++;
				}
			}else if(-P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				if(P.z > divHeight){
					points4[idx4] = P;
					idx4++;
				}else{
					points10[idx10] = P;
					idx10++;
				}
			}else{
				if(P.z > divHeight){
					points5[idx5] = P;
					idx5++;
				}else{
					points11[idx11] = P;
					idx11++;
				}
			}
		}else if(P.x > 0 && P.y < 0){
			if(-P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
				if(P.z > divHeight){
					points6[idx6] = P;
					idx6++;
				}else{
					points12[idx12] = P;
					idx12++;
				}
			}else if(-P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				if(P.z > divHeight){
					points6[idx6] = P;
					idx6++;
				}else{
					points12[idx12] = P;
					idx12++;
				}
			}else{
				if(P.z > divHeight){
					points5[idx5] = P;
					idx5++;
				}else{
					points11[idx11] = P;
					idx11++;
				}
			}
		}
	}

//	*sur = TopParametric(*c,*c,*c,*c,x[0],x[1],x[2],x[3]);
//	printf("dividing height: %f\n",h.z);

	printf("currently working with %d blocks with %d threads each\n",12,4);
    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;
    t_start = chrono::high_resolution_clock::now();

	kernel<<<12,4>>>(points1, points2, points3, points4, points5, points6, points7, points8, points9, points10, points11, points12,
					idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9, idx10, idx11, idx12,
					distances,N,sur1,sur2,sur3,sur4,sur5,sur6,sur7,sur8,sur9,sur10,sur11,sur12);
	cudaDeviceSynchronize();
	t_stop = chrono::high_resolution_clock::now();
	t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);

	cout << "Computation time: " << t_duration.count() << " microseconds" << endl;
	double sum = 0;
	for(int i = 0; i < N; i++){
		sum += distances[i];
	}
	printf("total distance: %lf\n", sum);
	printf("average distance: %lf\n", sum/N);
	printf("number of points: %d\n", N);

	CUDA_CHECK_RETURN(cudaFree(points1));
	CUDA_CHECK_RETURN(cudaFree(points2));
	CUDA_CHECK_RETURN(cudaFree(points3));
	CUDA_CHECK_RETURN(cudaFree(points4));
	CUDA_CHECK_RETURN(cudaFree(points5));
	CUDA_CHECK_RETURN(cudaFree(points6));
	CUDA_CHECK_RETURN(cudaFree(points7));
	CUDA_CHECK_RETURN(cudaFree(points8));
	CUDA_CHECK_RETURN(cudaFree(points9));
	CUDA_CHECK_RETURN(cudaFree(points10));
	CUDA_CHECK_RETURN(cudaFree(points11));
	CUDA_CHECK_RETURN(cudaFree(points12));

	CUDA_CHECK_RETURN(cudaFree(distances));
	CUDA_CHECK_RETURN(cudaFree(sur1));
	CUDA_CHECK_RETURN(cudaFree(sur2));
	CUDA_CHECK_RETURN(cudaFree(sur3));
	CUDA_CHECK_RETURN(cudaFree(sur4));
	CUDA_CHECK_RETURN(cudaFree(sur5));
	CUDA_CHECK_RETURN(cudaFree(sur6));
	CUDA_CHECK_RETURN(cudaFree(sur7));
	CUDA_CHECK_RETURN(cudaFree(sur8));
	CUDA_CHECK_RETURN(cudaFree(sur9));
	CUDA_CHECK_RETURN(cudaFree(sur10));
	CUDA_CHECK_RETURN(cudaFree(sur11));
	CUDA_CHECK_RETURN(cudaFree(sur12));
}

/**
 * Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 */
static void CheckCudaErrorAux (const char *file, unsigned line, const char *statement, cudaError_t err)
{
	if (err == cudaSuccess)
		return;
	std::cerr << statement<<" returned " << cudaGetErrorString(err) << "("<<err<< ") at "<<file<<":"<<line << std::endl;
	exit (1);
}

