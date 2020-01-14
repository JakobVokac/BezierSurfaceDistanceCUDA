//
// Created by s3179222 on 12/9/19.
//

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
