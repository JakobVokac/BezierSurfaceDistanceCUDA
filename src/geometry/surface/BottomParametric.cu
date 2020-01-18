//
// Created by s3179222 on 12/8/19.
//
/*
 * This is the implementation of the top part of each leaf in the model. The class contains all of the corners and edges
 * (see model.h), which include the bending curve, the bottom half of the symmetry curve and the bottom half
 * of the sinusoidal curve.
 *
 * The surface is an interpolation between the 3 edges. In the case of the bottom edge, the class returns an edge,
 * where all control points are Q_r.
 */
#include "BottomParametric.h"
__host__ __device__
BottomParametric::BottomParametric(cubiccrv sinCurve, cubiccrv symCurve, cubiccrv bendCurve, vec3d Q_r, vec3d Q_b_sin,
                                   vec3d Q_b_sym) {
    this->sinCurve = sinCurve;
    this->symCurve = symCurve;
    this->bendCurve = bendCurve;
    this->Q_r = Q_r;
    this->Q_b_sym = Q_b_sym;
    this->Q_b_sin = Q_b_sin;
    this->bottomCorner = cubiccrv(Q_r,Q_r,Q_r,Q_r);
}
__host__ __device__
vec3d BottomParametric::at(double u, double v) {
    vec3d res = v * Q_r + (1 - v) * bendCurve.f(u) + u * sinCurve.f(1 - v) + (1 - u) * symCurve.f(1 - v) -
                (u * v * Q_r + (1 - u) * (1 - v) * Q_b_sym + (1 - u) * v * Q_r + u * (1 - v) * Q_b_sin);
    return res;
}
__host__ __device__
vec3d BottomParametric::atDerU(double u, double v) {
    vec3d res = (1 - v) * bendCurve.df(u) + sinCurve.f(1 - v) - symCurve.f(1 - v) -
                (-(1 - v) * Q_b_sym + (1 - v) * Q_b_sin);
    return res;
}
__host__ __device__
vec3d BottomParametric::atDerV(double u, double v) {
    vec3d res = Q_r - bendCurve.f(u) - u * sinCurve.df(1 - v) - (1 - u) * symCurve.df(1 - v) -
                (u * Q_r - (1 - u) * Q_b_sym + (1 - u) * Q_r - u * Q_b_sin);
    return res;
}
__host__ __device__
vec3d BottomParametric::atDerUU(double u, double v) {
    vec3d res = (1 - v) * bendCurve.ddf(u);
    return res;
}
__host__ __device__
vec3d BottomParametric::atDerVV(double u, double v) {
    vec3d res = u * sinCurve.ddf(1 - v) + (1 - u) * symCurve.ddf(1 - v);
    return res;
}
__host__ __device__
vec3d BottomParametric::atDerUV(double u, double v) {
    vec3d res = - bendCurve.df(u) - sinCurve.f(1 - v) + symCurve.f(1 - v) -
                (Q_b_sym - Q_b_sin);
    return res;
}
__host__ __device__
cubiccrv & BottomParametric::edgeU0() {
    return symCurve;
}
__host__ __device__
cubiccrv & BottomParametric::edgeU1() {
    return sinCurve;
}
__host__ __device__
cubiccrv & BottomParametric::edgeV0() {
    return bendCurve;
}
__host__ __device__
cubiccrv & BottomParametric::edgeV1() {
    return bottomCorner;
}

__host__ __device__
double BottomParametric::sqDistTo(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    double dist = top.sqdist(A);
    return dist;
}
__host__ __device__
double BottomParametric::sqDistToDerU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerU(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double BottomParametric::sqDistToDerV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerV(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double BottomParametric::sqDistToDerUU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerU(u,v);
    vec3d topDer2 = atDerUU(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double BottomParametric::sqDistToDerVV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerV(u,v);
    vec3d topDer2 = atDerVV(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double BottomParametric::sqDistToDerUV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDerU = atDerU(u,v);
    vec3d topDerV = atDerV(u,v);
    vec3d topDerUV = atDerUV(u,v);

    double dist = 2*(((top - A)*topDerUV).sum() + (topDerU*topDerV).sum());

    return dist;
}
__host__ __device__
double BottomParametric::distTo(double u, double v, class vec3d A) {
    return sqrt(sqDistTo(u,v,A));
}
__host__ __device__
double BottomParametric::distToDerU(double u, double v, vec3d A){

    return sqDistToDerU(u,v,A) / distTo(u,v,A);
}
__host__ __device__
double BottomParametric::distToDerV(double u, double v, vec3d A) {

    return sqDistToDerV(u, v, A) / distTo(u, v, A);
}
__host__ __device__
double BottomParametric::distToDerUU(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerU(u,v);
    vec3d f2 = atDerUU(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double BottomParametric::distToDerVV(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerV(u,v);
    vec3d f2 = atDerVV(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double BottomParametric::distToDerUV(double u, double v, vec3d A) {
    vec3d f = at(u, v);
    vec3d f1u = atDerU(u, v);
    vec3d f1v = atDerV(u, v);
    vec3d f2 = atDerUV(u, v);

    double dist = ((((f - A)*(f - A)).sum())*((f2*(f-A)).sum() + (f1u*f1v).sum()) - (f1u*(f-A)).sum()*(f1v*(f-A)).sum())/(pow(((f-A)*(f-A)).sum(),1.5));

    return dist;
}
