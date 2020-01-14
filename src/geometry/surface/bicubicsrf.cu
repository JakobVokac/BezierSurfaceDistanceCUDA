//
// Created by s3179222 on 12/4/19.
//

#include "bicubicsrf.h"

int factorial(int n){
    int f = 1;
    for (int i = 1; i <= n; ++i) {
        f *= i;
    }
    return f;
}
int binom(int n, int i){
    if(n < 2){
        return 1;
    }
    return factorial(n)/(factorial(i)*factorial(n-i));
}

vec3d bicubicsrf::at(double u, double v){
    vec3d res = {0,0,0};
    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            double param = binom(3,i)*pow(u,i)*pow(1-u,3-i)*binom(3,j)*pow(v,j)*pow(1-v,3-j);
            res += ctrl[i*4 + j]*param;
        }
    }
    return res;
};
vec3d bicubicsrf::atDerU(double u, double v){
    vec3d res = {0,0,0};
    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j <= 3; ++j) {
            double param = binom(3, j) * pow(v, j) * pow(1 - v, 3 - j) *
                           binom(2, i) * pow(u, i) * pow(1 - u, 2 - i);
            res += (ctrl[(i + 1) * 4 + j] - ctrl[i * 4 + j]) * 3 * param;
        }
    }
    return res;
};
vec3d bicubicsrf::atDerV(double u, double v){
    vec3d res = {0,0,0};
    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j <= 3; ++j) {
            double param = binom(2,j)*pow(v,j)*pow(1-v,2-j)*
                           binom(3,i)*pow(u,i)*pow(1-u,3-i);
            res += (ctrl[i*4 + j+1] - ctrl[i*4 + j]) * 3 * param;
        }
    }
    return res;
};
vec3d bicubicsrf::atDerUU(double u, double v){
    vec3d res = {0,0,0};
    for (int i = 0; i <= 1; ++i) {
        for (int j = 0; j <= 3; ++j) {
            double param = binom(3, j) * pow(v, j) * pow(1 - v, 3 - j) *
                           binom(1, i) * pow(u, i) * pow(1 - u, 1 - i);
            res += (ctrl[(i + 2) * 4 + j] - ctrl[(i + 1) * 4 + j]*2 + ctrl[i * 4 + j]) * 6 * param;
        }
    }
    return res;
};
vec3d bicubicsrf::atDerVV(double u, double v){
    vec3d res = {0,0,0};
    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j <= 1; ++j) {
            double param = binom(1, j) * pow(v, j) * pow(1 - v, 1 - j) *
                           binom(3, i) * pow(u, i) * pow(1 - u, 3 - i);
            res += (ctrl[i * 4 + j + 2] - ctrl[i * 4 + j + 1]*2 + ctrl[i * 4 + j]) * 6 * param;
        }
    }
    return res;
};
vec3d bicubicsrf::atDerUV(double u, double v){
    vec3d res = {0, 0, 0};
    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j <= 2; ++j) {
            double param = binom(2, j) * pow(v, j) * pow(1 - v, 2 - j) *
                           binom(2, i) * pow(u, i) * pow(1 - u, 2 - i);
            res += (ctrl[(i + 1) * 4 + j + 1] - ctrl[(i + 1) * 4 + j] - ctrl[i * 4 + j + 1] + ctrl[i * 4 + j]) * 9 * param;
        }
    }
    return res;
}
__host__ __device__
cubiccrv &bicubicsrf::edgeU0() {
    return U0;
}
__host__ __device__
cubiccrv &bicubicsrf::edgeU1() {
    return U1;
}
__host__ __device__
cubiccrv &bicubicsrf::edgeV0() {
    // TODO: check orientation because of parametric
    return V0;
}
__host__ __device__
cubiccrv &bicubicsrf::edgeV1() {
    // TODO: check orientation because of parametric
    return V1;
};

bool bicubicsrf::hasValidControlNet(){
    for (int i = 0; i < 3; ++i) {
        cubiccrv c = {ctrl[i*4 + 0],ctrl[i*4 + 1],ctrl[i*4 + 2],ctrl[i*4 + 3]};
        if(!c.hasConvexPolygon())
            return false;
    }
    for (int i = 0; i < 3; ++i) {
        cubiccrv c = {ctrl[i + 0],ctrl[i + 4],ctrl[i + 8],ctrl[i + 12]};
        if(!c.hasConvexPolygon())
            return false;
    }
    return true;
}

bool bicubicsrf::closestPointInPatch(vec3d P){
    bool flag = false;
    for (int i = 0; i < 3; ++i) {
        cubiccrv c = {ctrl[i*4 + 0],ctrl[i*4 + 1],ctrl[i*4 + 2],ctrl[i*4 + 3]};
        if(c.closestPointInCurve(P)) {
            flag = true;
            break;
        }
    }
    if(!flag)
        return false;

    flag = false;

    for (int i = 0; i < 3; ++i) {
        cubiccrv c = {ctrl[i + 0],ctrl[i + 4],ctrl[i + 8],ctrl[i + 12]};
        if(c.closestPointInCurve(P)) {
            flag = true;
            break;
        }
    }

    if(!flag)
        return false;

    return true;
}

void bicubicsrf::subdivideInDir(bool dir, double t, bicubicsrf &srf1, bicubicsrf &srf2) {
    if(dir){
        // direction u
        for (int i = 0; i < 4; ++i) {
            cubiccrv c = {ctrl[i + 0], ctrl[i + 4], ctrl[i + 8], ctrl[i + 12]};
            cubiccrv c1, c2;

            c.subdivide(t, c1, c2);

            for (int j = 0; j < 4; ++j) {
                srf1.ctrl[i + 4 * j] = c1.getCtrlP(j);
                srf2.ctrl[i + 4 * j] = c2.getCtrlP(j);
            }
        }
    }else{
        // direction v
        for (int i = 0; i < 4; ++i) {
            cubiccrv c = {ctrl[i*4 + 0], ctrl[i*4 + 1], ctrl[i*4 + 2], ctrl[i*4 + 3]};
            cubiccrv c1, c2;
            c.subdivide(t, c1, c2);
            for (int j = 0; j < 4; ++j) {
                srf1.ctrl[i + 4 * j] = c1.getCtrlP(j);
                srf2.ctrl[i + 4 * j] = c2.getCtrlP(j);
            }
        }
    }

}

vec3d bicubicsrf::ctrlP(int i) {
    if(i < 0 || i > 15){
        std::cout << "bicubic mySurface control point index out of bounds" << std::endl;
        return {NAN,NAN,NAN};
    }
    return ctrl[i];
}

void
bicubicsrf::subdivide(bicubicsrf &tl, bicubicsrf &tr, bicubicsrf &bl, bicubicsrf &br) {
    bicubicsrf s1{},s2{};

    subdivideInDir(true, 0.5, s1, s2);
    s1.subdivideInDir(false, 0.5, tl, bl);

    s2.subdivideInDir(false, 0.5, tr, br);

}


