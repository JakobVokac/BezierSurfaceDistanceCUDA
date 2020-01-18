//
// Created by s3179222 on 12/8/19.
//
/*
 * This is the implementation of the top part of each leaf in the model. The class contains all of the corners and edges
 * (see model.h), which include the bending curve, the leaf curve, the top half of the symmetry curve and the top half
 * of the sinusoidal curve.
 *
 * The surface itself is a bilinear interpolation of the 4 edges.
 */

#ifndef HEARTVALVEMODEL_TOPPARAMETRIC_H
#define HEARTVALVEMODEL_TOPPARAMETRIC_H


#include "../curve/cubiccrv.h"

class TopParametric {
public:
    TopParametric() = default;
    __host__ __device__ TopParametric(
            cubiccrv leafCurve,
            cubiccrv bendCurve,
            cubiccrv symCurve,
            cubiccrv sinCurve,
            vec3d Q_l_sin,
            vec3d Q_l_sym,
            vec3d Q_b_sin,
            vec3d Q_b_sym
    );

    ~TopParametric() = default;

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

private:
    cubiccrv leafCurve, bendCurve, symCurve, sinCurve;

    vec3d Q_l_sin, Q_l_sym, Q_b_sin, Q_b_sym;
};


#endif //HEARTVALVEMODEL_TOPPARAMETRIC_H
