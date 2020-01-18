//
// Created by s3179222 on 12/9/19.
//
/*
 * This class is the implementation of the Newton-Raphson method for the numerical optimizer.
 */
#ifndef HEARTVALVEMODEL_NEWTON_H
#define HEARTVALVEMODEL_NEWTON_H

#include "../states.h"
#include "../../geometry/surface/TopParametric.h"
#include "../../geometry/surface/BottomParametric.h"

class Newton{
private:
    double sigma;
    double Gu, Gv;
    double Huu, Huv, Hvv;
public:
    __host__ __device__ Newton(double sigma){
        this->sigma = sigma;
        Gu = 0;
        Gv = 0;
        Huu = 0;
        Huv = 0;
        Hvv = 0;
    }
    ~Newton() = default;
    __host__ __device__ OptState2D doStep(TopParametric *sur, vec3d &P, OptState2D &loc);
    __host__ __device__ OptState2D doStep(BottomParametric *sur, vec3d &P, OptState2D &loc);
    __host__ __device__ OptState1D doStep(cubiccrv &crv, vec3d &P, OptState1D &loc);
};


#endif //HEARTVALVEMODEL_NEWTON_H
