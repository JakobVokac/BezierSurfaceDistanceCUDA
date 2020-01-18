//
// Created by s3179222 on 12/4/19.
//
/*
 * This is an implementation of a 3rd degree Bezier surface. Which contains all of the 16 control points, the actual functions,
 * and some additional function for detection whether its control net is valid as well as a function for checking
 * whether the closest point on the surface to a point P is inside the surface or on an edge.
 *
 * The net and edge detection algorithms were taken and implemented from:
 *      Yingliang Ma & W. Terry Hewitt,
 *      "Point inversion and projection for NURBS curve: Control polygon approach",
 *      https://www.researchgate.net/publication/232638188_Point_inversion_and_projection_for_NURBS_curve_Control_polygon_approach
 */
#ifndef HEARTVALVEMODEL_BICUBICSRF_H
#define HEARTVALVEMODEL_BICUBICSRF_H

#include "mySurface.h"

class bicubicsrf : public mySurface {
private:
    vec3d ctrl[16]{{0,0,0}};
    cubiccrv U0, V0, U1, V1;
public:
    bicubicsrf() = default;
    __host__ __device__ bicubicsrf(vec3d ctrl[16]){
        for (int i = 0; i < 16; ++i) {
            this->ctrl[i] = ctrl[i];
        }
        this->U0 = cubiccrv(ctrl[0],ctrl[1],ctrl[2],ctrl[3]);
        this->V0 = cubiccrv(ctrl[0],ctrl[4],ctrl[8],ctrl[12]);
        this->V1 = cubiccrv(ctrl[3],ctrl[7],ctrl[11],ctrl[15]);
        this->U1 = cubiccrv(ctrl[12],ctrl[13],ctrl[14],ctrl[15]);
    }
    ~bicubicsrf() = default;

    __host__ __device__ vec3d at(double u, double v);
    __host__ __device__ vec3d atDerU(double u, double v);
    __host__ __device__ vec3d atDerV(double u, double v);
    __host__ __device__ vec3d atDerUU(double u, double v);
    __host__ __device__ vec3d atDerVV(double u, double v);
    __host__ __device__ vec3d atDerUV(double u, double v);

    __host__ __device__ cubiccrv & edgeU0();
    __host__ __device__ cubiccrv & edgeU1();
    __host__ __device__ cubiccrv & edgeV0();
    __host__ __device__ cubiccrv & edgeV1();

    vec3d ctrlP(int i);
    bool hasValidControlNet();
    bool closestPointInPatch(vec3d P);
    void subdivideInDir(bool dir, double t, bicubicsrf &srf1, bicubicsrf &srf2);
    void subdivide(bicubicsrf &tl, bicubicsrf &tr, bicubicsrf &bl, bicubicsrf &br);

};


#endif //HEARTVALVEMODEL_BICUBICSRF_H
