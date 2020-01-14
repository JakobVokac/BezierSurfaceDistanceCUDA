//
// Created by s3179222 on 12/8/19.
//

#ifndef HEARTVALVEMODEL_BOTTOMPARAMETRIC_H
#define HEARTVALVEMODEL_BOTTOMPARAMETRIC_H

#include "../curve/cubiccrv.h"

class BottomParametric {
public:
    BottomParametric() = default;
    BottomParametric(
            cubiccrv sinCurve,
            cubiccrv symCurve,
            cubiccrv bendCurve,
            vec3d Q_r,
            vec3d Q_b_sin,
            vec3d Q_b_sym
            );
    ~BottomParametric() = default;

    __host__ __device__ vec3d at(double u, double v);
    __host__ __device__ vec3d atDerU(double u, double v);
    __host__ __device__ vec3d atDerV(double u, double v);
    __host__ __device__ vec3d atDerUU(double u, double v);
    __host__ __device__ vec3d atDerVV(double u, double v);
    __host__ __device__ vec3d atDerUV(double u, double v);

    __host__ __device__ double sqDistTo(double u, double v, vec3d A);
    __host__ __device__ double sqDistToDerU(double u, double v, vec3d A);
    __host__ __device__ double sqDistToDerV(double u, double v, vec3d A);
    __host__ __device__ double sqDistToDerUU(double u, double v, vec3d A);
    __host__ __device__ double sqDistToDerVV(double u, double v, vec3d A);
    __host__ __device__ double sqDistToDerUV(double u, double v, vec3d A);

    __host__ __device__ double distTo(double u, double v, vec3d A);
    __host__ __device__ double distToDerU(double u, double v, vec3d A);
    __host__ __device__ double distToDerV(double u, double v, vec3d A);
    __host__ __device__ double distToDerUU(double u, double v, vec3d A);
    __host__ __device__ double distToDerVV(double u, double v, vec3d A);
    __host__ __device__ double distToDerUV(double u, double v, vec3d A);

    __host__ __device__ cubiccrv & edgeU0();
    __host__ __device__ cubiccrv & edgeU1();
    __host__ __device__ cubiccrv & edgeV0();
    __host__ __device__ cubiccrv & edgeV1();

    __host__ __device__ vec3d getBottomCorner(){
    	return Q_r;
    }
private:
    cubiccrv sinCurve, symCurve, bendCurve, bottomCorner;
    vec3d Q_r, Q_b_sin, Q_b_sym;

};


#endif //HEARTVALVEMODEL_BOTTOMPARAMETRIC_H
