//
// Created by s3179222 on 12/16/19.
//

#include "compositeBicubicsrf.h"

__host__ __device__
bool compositeBicubicsrf::closestPointInPatch(vec3d P) {
    if(hasValidControlNet()) {

        if(topLeft == nullptr)
            return srf.closestPointInPatch(P);

        return topLeft->closestPointInPatch(P) ||
               topRight->closestPointInPatch(P) ||
               bottomLeft->closestPointInPatch(P) ||
               bottomRight->closestPointInPatch(P);
    }
    return false;
}
__host__ __device__
bool compositeBicubicsrf::hasValidControlNet() {
    if(topLeft == nullptr)
        return srf.hasValidControlNet();

    return topLeft->hasValidControlNet() ||
           topRight->hasValidControlNet() ||
           bottomLeft->hasValidControlNet() ||
           bottomRight->hasValidControlNet();
}
__host__ __device__
vec3d compositeBicubicsrf::at(double u, double v) {
    if(topLeft == nullptr)
        return srf.at(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->at(u*2,v*2);
        }else{
            return topLeft->at(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->at((u-0.5)*2,v*2);
        }else{
            return topRight->at((u-0.5)*2,(v-0.5)*2);
        }
    }
}
__host__ __device__
vec3d compositeBicubicsrf::atDerU(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerU(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerU(u*2,v*2);
        }else{
            return topLeft->atDerU(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerU((u-0.5)*2,v*2);
        }else{
            return topRight->atDerU((u-0.5)*2,(v-0.5)*2);
        }
    }
}
__host__ __device__
vec3d compositeBicubicsrf::atDerV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerV(u*2,v*2);
        }else{
            return topLeft->atDerV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerV((u-0.5)*2,(v-0.5)*2);
        }
    }
}
__host__ __device__
vec3d compositeBicubicsrf::atDerUU(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerUU(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerUU(u*2,v*2);
        }else{
            return topLeft->atDerUU(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerUU((u-0.5)*2,v*2);
        }else{
            return topRight->atDerUU((u-0.5)*2,(v-0.5)*2);
        }
    }
}
__host__ __device__
vec3d compositeBicubicsrf::atDerVV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerVV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerVV(u*2,v*2);
        }else{
            return topLeft->atDerVV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerVV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerVV((u-0.5)*2,(v-0.5)*2);
        }
    }
}
__host__ __device__
vec3d compositeBicubicsrf::atDerUV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerUV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerUV(u*2,v*2);
        }else{
            return topLeft->atDerUV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerUV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerUV((u-0.5)*2,(v-0.5)*2);
        }
    }
}


__host__ __device__
cubiccrv &compositeBicubicsrf::edgeU0() {
    return srf.edgeU0();
}
__host__ __device__
cubiccrv &compositeBicubicsrf::edgeU1() {
    return srf.edgeU1();
}
__host__ __device__
cubiccrv &compositeBicubicsrf::edgeV0() {
    return srf.edgeV0();
}
__host__ __device__
cubiccrv &compositeBicubicsrf::edgeV1() {
    return srf.edgeV1();
}

