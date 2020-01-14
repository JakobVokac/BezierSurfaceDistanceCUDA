//
// Created by s3179222 on 12/8/19.
//

#include "TopParametric.h"
__host__ __device__
TopParametric::TopParametric(cubiccrv leafCurve, cubiccrv bendCurve, cubiccrv symCurve, cubiccrv sinCurve,
                             vec3d Q_l_sin, vec3d Q_l_sym, vec3d Q_b_sin, vec3d Q_b_sym) {
    this->leafCurve = leafCurve;
    this->bendCurve = bendCurve;
    this->sinCurve = sinCurve;
    this->symCurve = symCurve;

    this->Q_l_sin = Q_l_sin;
    this->Q_l_sym = Q_l_sym;
    this->Q_b_sin = Q_b_sin;
    this->Q_b_sym = Q_b_sym;
}
__host__ __device__
vec3d TopParametric::at(double u, double v) {
    vec3d point = (1 - v) * bendCurve.f(u) + v * leafCurve.f(u) +
                  (1 - u) * symCurve.f(1 - v) + u * sinCurve.f(1 - v) -
                  ((1 - u) * (1 - v) * Q_b_sym + u * v * Q_l_sin +
                   u * (1 - v) * Q_b_sin + (1 - u) * v * Q_l_sym);

    return point;
}
__host__ __device__
vec3d TopParametric::atDerU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.df(u) + v * leafCurve.df(u)
                  - symCurve.f(1 - v) + sinCurve.f(1 - v) -
                  (- (1 - v) * Q_b_sym + v * Q_l_sin +
                   (1 - v) * Q_b_sin - v * Q_l_sym);

    return point;
}
__host__ __device__
vec3d TopParametric::atDerV(double u, double v) {
    vec3d point = - bendCurve.f(u) + leafCurve.f(u) -
                  (1 - u) * symCurve.df(1 - v) - u * sinCurve.df(1 - v) -
                  (- (1 - u) * Q_b_sym + u * Q_l_sin
                   - u * Q_b_sin + (1 - u) * Q_l_sym);


    return point;
}
__host__ __device__
vec3d TopParametric::atDerUU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.ddf(u) + v * leafCurve.ddf(u);

    return point;
}
__host__ __device__
vec3d TopParametric::atDerVV(double u, double v) {
    vec3d point = (1 - u) * symCurve.ddf(1 - v) + u * sinCurve.ddf(1 - v);

    return point;
}
__host__ __device__
vec3d TopParametric::atDerUV(double u, double v) {
    vec3d point = - bendCurve.df(u) + leafCurve.df(u)
                  + symCurve.df(1 - v) - sinCurve.df(1 - v) -
                  (Q_b_sym + Q_l_sin - Q_b_sin - Q_l_sym);

    return point;
}
__host__ __device__
cubiccrv & TopParametric::edgeU0() {
    return symCurve;
}
__host__ __device__
cubiccrv & TopParametric::edgeU1() {
    return sinCurve;
}
__host__ __device__
cubiccrv & TopParametric::edgeV0() {
    return bendCurve;
}
__host__ __device__
cubiccrv & TopParametric::edgeV1() {
    return leafCurve;
}

__host__ __device__
double TopParametric::sqDistTo(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    double dist = top.sqdist(A);
    return dist;
}
__host__ __device__
double TopParametric::sqDistToDerU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerU(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double TopParametric::sqDistToDerV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerV(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double TopParametric::sqDistToDerUU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerU(u,v);
    vec3d topDer2 = atDerUU(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double TopParametric::sqDistToDerVV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerV(u,v);
    vec3d topDer2 = atDerVV(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double TopParametric::sqDistToDerUV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDerU = atDerU(u,v);
    vec3d topDerV = atDerV(u,v);
    vec3d topDerUV = atDerUV(u,v);

    double dist = 2*(((top - A)*topDerUV).sum() + (topDerU*topDerV).sum());

    return dist;
}
__host__ __device__
double TopParametric::distTo(double u, double v, class vec3d A) {
    return sqrt(sqDistTo(u,v,A));
}
__host__ __device__
double TopParametric::distToDerU(double u, double v, vec3d A){

    return sqDistToDerU(u,v,A) / distTo(u,v,A);
}
__host__ __device__
double TopParametric::distToDerV(double u, double v, vec3d A) {

    return sqDistToDerV(u, v, A) / distTo(u, v, A);
}
__host__ __device__
double TopParametric::distToDerUU(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerU(u,v);
    vec3d f2 = atDerUU(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double TopParametric::distToDerVV(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerV(u,v);
    vec3d f2 = atDerVV(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double TopParametric::distToDerUV(double u, double v, vec3d A) {
    vec3d f = at(u, v);
    vec3d f1u = atDerU(u, v);
    vec3d f1v = atDerV(u, v);
    vec3d f2 = atDerUV(u, v);

    double dist = ((((f - A)*(f - A)).sum())*((f2*(f-A)).sum() + (f1u*f1v).sum()) - (f1u*(f-A)).sum()*(f1v*(f-A)).sum())/(pow(((f-A)*(f-A)).sum(),1.5));

    return dist;
}


