//
// Created by s3179222 on 12/10/19.
//

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
