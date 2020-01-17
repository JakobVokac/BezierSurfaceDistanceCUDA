/*
 ============================================================================
 Name        : CUDANewton.cu
 Author      : Jakob Vokac
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
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
/**
 * CUDA kernel that computes reciprocal values for a given vector
 */


__global__ void kernel(vec3d *points, double *distances, int size,
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
		BottomParametric *sur12,
		double divHeight
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
	__shared__ TopParametric top1, top2, top3, top4, top5, top6;
	__shared__ BottomParametric bot1, bot2, bot3, bot4, bot5, bot6;

	if(threadIdx.x == 0)
		top1 = *sur1;
	if(threadIdx.x == 1)
		top2 = *sur2;
	if(threadIdx.x == 2)
		top3 = *sur3;
	if(threadIdx.x == 3)
		top4 = *sur4;
	if(threadIdx.x == 4)
		top5 = *sur5;
	if(threadIdx.x == 5)
		top6 = *sur6;
	if(threadIdx.x == 6)
		bot1 = *sur7;
	if(threadIdx.x == 7)
		bot2 = *sur8;
	if(threadIdx.x == 8)
		bot3 = *sur9;
	if(threadIdx.x == 9)
		bot4 = *sur10;
	if(threadIdx.x == 10)
		bot5 = *sur11;
	if(threadIdx.x == 11)
		bot6 = *sur12;

	__syncthreads();
	for(int i = 0; i < size; i+= 1024){
		int arrId = threadIdx.x + blockDim.x * blockIdx.x + i;

		if(arrId >= size)
			break;
		vec3d P = points[arrId];

		op.setTopOrBot(P.z > divHeight);

		if(P.x == 0){
			if(P.y == 0){
				op.setTop(&top1);
				op.setBot(&bot1);
			}
			else if(P.y > 0){
				op.setTop(&top2);
				op.setBot(&bot2);
			}
			else if(P.y < 0){
				op.setTop(&top5);
				op.setBot(&bot5);
			}
		}

		if(P.x > 0 && P.y > 0){
			if(P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
				op.setTop(&top1);
				op.setBot(&bot1);
			}else if(P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				op.setTop(&top1);
				op.setBot(&bot1);
			}else{
				op.setTop(&top2);
				op.setBot(&bot2);
			}
		}else if(P.x < 0 && P.y < 0){
			if(P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
				op.setTop(&top3);
				op.setBot(&bot3);
			}else if(P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				op.setTop(&top3);
				op.setBot(&bot3);
			}else{
				op.setTop(&top2);
				op.setBot(&bot2);
			}
		}else if(P.x < 0 && P.y < 0){
			if(-P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
				op.setTop(&top4);
				op.setBot(&bot4);
			}else if(-P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				op.setTop(&top4);
				op.setBot(&bot4);
			}else{
				op.setTop(&top5);
				op.setBot(&bot5);
			}
		}else if(P.x > 0 && P.y < 0){
			if(-P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
				op.setTop(&top6);
				op.setBot(&bot6);
			}else if(-P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
				op.setTop(&top6);
				op.setBot(&bot6);
			}else{
				op.setTop(&top5);
				op.setBot(&bot5);
			}
		}
		OptState2D loc = op.optimizeForPoint(P);
//		printf("threadidx.x: %d, blockDim.x: %d, blockIdx.x: %d, i: %d, arrId: %d, dist: %lf\n",threadIdx.x, blockDim.x, blockIdx.x, i, arrId, dist);

		distances[arrId] = loc.dist;
	}
}



//
///**
// * Host function that copies the data and launches the work on GPU
// */
//float *gpuReciprocal(float *data, unsigned size)
//{
//	float *rc = new float[size];
//	float *gpuData;
//
//	CUDA_CHECK_RETURN(cudaMalloc((void **)&gpuData, sizeof(float)*size));
//	CUDA_CHECK_RETURN(cudaMemcpy(gpuData, data, sizeof(float)*size, cudaMemcpyHostToDevice));
//
//	static const int BLOCK_SIZE = 256;
//	const int blockCount = (size+BLOCK_SIZE-1)/BLOCK_SIZE;
//	reciprocalKernel<<<blockCount, BLOCK_SIZE>>> (gpuData, size);
//
//	CUDA_CHECK_RETURN(cudaMemcpy(rc, gpuData, sizeof(float)*size, cudaMemcpyDeviceToHost));
//	CUDA_CHECK_RETURN(cudaFree(gpuData));
//	return rc;
//}
//
//float *cpuReciprocal(float *data, unsigned size)
//{
//	float *rc = new float[size];
//	for (unsigned cnt = 0; cnt < size; ++cnt) rc[cnt] = 1.0/data[cnt];
//	return rc;
//}
//
//
//void initialize(float *data, unsigned size)
//{
//	for (unsigned i = 0; i < size; ++i)
//		data[i] = .5*(i+1);
//}


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

	int nBlocks = 4;
	int nThreads = 32;
	int N = inputPoints.size()/3;
	double dist = 0.5;
	vec3d *points;
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
	CUDA_CHECK_RETURN(cudaMallocManaged(&points, N*sizeof(vec3d)));
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

	// initialize x and y arrays on the host
	for (int i = 0; i < N; i++) {

        vec3d a = {inputPoints[i*3 + 0],inputPoints[i*3 + 1],inputPoints[i*3 + 2]};

		points[i] = a;
	}

//	*sur = TopParametric(*c,*c,*c,*c,x[0],x[1],x[2],x[3]);
	vec3d h = sur1->at(0,0);
//	printf("dividing height: %f\n",h.z);

	printf("currently working with %d blocks with %d threads each\n",nBlocks,nThreads);

	kernel<<<nBlocks,nThreads>>>(points,distances,N,sur1,sur2,sur3,sur4,sur5,sur6,sur7,sur8,sur9,sur10,sur11,sur12,h.z);
	cudaDeviceSynchronize();

	double sum = 0;
	for(int i = 0; i < N; i++){
		sum += distances[i];
	}
	printf("total distance: %lf\n", sum);
	printf("average distance: %lf\n", sum/N);
	printf("number of points: %d\n", N);

	CUDA_CHECK_RETURN(cudaFree(points));
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

