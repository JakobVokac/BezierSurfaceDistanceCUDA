//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the binary-search/bisection implementation for preprocessing.
 *
 * The 1D version uses basic bisection.
 *
 * The 2D version uses a binary search, which works better for surfaces. It splits the surface into two parts per parameter,
 * takes the median of each side and continues with the closer side to point P. It alternates between u and v.
 *
 * Along further work, the binary search should be cut out of this class and placed in its own class and a proper 2D
 * bisection method should be implemented here.
 */
#ifndef HEARTVALVEMODEL_BISECTION_H
#define HEARTVALVEMODEL_BISECTION_H

#include "../states.h"
#include "../../geometry/surface/TopParametric.h"
#include "../../geometry/surface/BottomParametric.h"

class bisection{
private:
    int iterations;
public:
    __host__ __device__ bisection(int iter){
        this->iterations = iter;
    };
    ~bisection() = default;
    __host__ __device__ OptState2D preprocess(TopParametric *sur, const vec3d &P);
    __host__ __device__ OptState2D preprocess(BottomParametric *sur, const vec3d &P);
    __host__ __device__ OptState1D preprocess(cubiccrv &crv, const vec3d &P);
};


#endif //HEARTVALVEMODEL_BISECTION_H
