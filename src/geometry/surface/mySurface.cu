//
// Created by s3179222 on 12/8/19.
//
/*
 * This class is a translation of the surface interface. Currently it has no use.
 */
#include "mySurface.h"

__host__ __device__
double mySurface::sqDistTo(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    double dist = top.sqdist(A);
    return dist;
}
__host__ __device__
double mySurface::sqDistToDerU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerU(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double mySurface::sqDistToDerV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerV(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
__host__ __device__
double mySurface::sqDistToDerUU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerU(u,v);
    vec3d topDer2 = atDerUU(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double mySurface::sqDistToDerVV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerV(u,v);
    vec3d topDer2 = atDerVV(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
__host__ __device__
double mySurface::sqDistToDerUV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDerU = atDerU(u,v);
    vec3d topDerV = atDerV(u,v);
    vec3d topDerUV = atDerUV(u,v);

    double dist = 2*(((top - A)*topDerUV).sum() + (topDerU*topDerV).sum());

    return dist;
}
__host__ __device__
double mySurface::distTo(double u, double v, class vec3d A) {
    return sqrt(sqDistTo(u,v,A));
}
__host__ __device__
double mySurface::distToDerU(double u, double v, vec3d A){

    return sqDistToDerU(u,v,A) / distTo(u,v,A);
}
__host__ __device__
double mySurface::distToDerV(double u, double v, vec3d A) {

    return sqDistToDerV(u, v, A) / distTo(u, v, A);
}
__host__ __device__
double mySurface::distToDerUU(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerU(u,v);
    vec3d f2 = atDerUU(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double mySurface::distToDerVV(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerV(u,v);
    vec3d f2 = atDerVV(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}
__host__ __device__
double mySurface::distToDerUV(double u, double v, vec3d A) {
    vec3d f = at(u, v);
    vec3d f1u = atDerU(u, v);
    vec3d f1v = atDerV(u, v);
    vec3d f2 = atDerUV(u, v);

    double dist = ((((f - A)*(f - A)).sum())*((f2*(f-A)).sum() + (f1u*f1v).sum()) - (f1u*(f-A)).sum()*(f1v*(f-A)).sum())/(pow(((f-A)*(f-A)).sum(),1.5));

    return dist;
}
