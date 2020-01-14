//
// Created by s3179222 on 12/4/19.
//

#ifndef HEARTVALVEMODEL_CUBICCRV_H
#define HEARTVALVEMODEL_CUBICCRV_H

#include <iostream>
#include "../vector/vec3d.h"
class cubiccrv{
private:
    vec3d c0,c1,c2,c3;

public:
    cubiccrv() = default;

    __host__ __device__ cubiccrv(vec3d c0, vec3d c1, vec3d c2, vec3d c3){
        this->c0 = c0;
        this->c1 = c1;
        this->c2 = c2;
        this->c3 = c3;
    }

    __host__ __device__ cubiccrv(vec3d c0, vec3d c1, vec3d c2){
        this->c0 = c0;
        this->c1 = c0 * 1.0/3.0 + c1 * 2.0/3.0;
        this->c2 = c1 * 2.0/3.0 + c2 * 1.0/3.0;
        this->c3 = c2;
    }
    ~cubiccrv() = default;

    __host__ __device__ vec3d f(double t);
    __host__ __device__ vec3d df(double t);
    __host__ __device__ vec3d ddf(double t);

    __host__ __device__ double sqDistTo(double t, vec3d A);
    __host__ __device__ double sqDistToDer1(double t, vec3d A);
    __host__ __device__ double sqDistToDer2(double t, vec3d A);

    __host__ __device__ double distTo(double t, vec3d A);
    __host__ __device__ double distToDer1(double t, vec3d A);
    __host__ __device__ double distToDer2(double t, vec3d A);

    __host__ __device__ vec3d curvePlaneNormal();
    __host__ __device__ vec3d getCtrlP(int i);

    __host__ __device__ cubiccrv RotAxisZ(double rad);
    __host__ __device__ cubiccrv MirrorY();

    std::string asString();

    __host__ __device__ bool hasConvexPolygon();
    __host__ __device__ bool closestPointInCurve(vec3d P);
    __host__ __device__ void subdivide(double t, cubiccrv &crv1, cubiccrv &crv2);
};



#endif //HEARTVALVEMODEL_CUBICCRV_H
