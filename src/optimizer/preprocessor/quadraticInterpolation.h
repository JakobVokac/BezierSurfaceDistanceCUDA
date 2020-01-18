//
// Created by s3179222 on 12/10/19.
//
/*
 * This is the implementation of the quadratic interpolation method.
 *
 * The 1D version simply interpolates from three starting locations (t = 0,0.5,1) on the edge and replaces the worst
 * point by the minimum of the quadratic approximation of the curve. It continues this for a specified number of iterations.
 *
 * The 2D method alternates between parameters u and v and makes an interpolation along one, while using the current best
 * solution for the other parameter as a constant. It continues this for the specified number of iterations for each parameter.
 */
#ifndef HEARTVALVEMODEL_QUADRATICINTERPOLATION_H
#define HEARTVALVEMODEL_QUADRATICINTERPOLATION_H


#include "../states.h"
#include "../../geometry/surface/TopParametric.h"
#include "../../geometry/surface/BottomParametric.h"

class quadraticInterpolation{
    int iterations;
public:
    __host__ __device__ quadraticInterpolation(int iter){
        this->iterations = iter;
    };
    ~quadraticInterpolation() = default;
    __host__ __device__ OptState2D preprocess(TopParametric &sur, const vec3d &P);
    __host__ __device__ OptState2D preprocess(BottomParametric &sur, const vec3d &P);

    __host__ __device__ OptState1D preprocess(cubiccrv &crv, const vec3d &P);
};


#endif //HEARTVALVEMODEL_QUADRATICINTERPOLATION_H
