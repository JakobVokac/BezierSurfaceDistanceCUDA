//
// Created by s3179222 on 12/9/19.
//

#include "Newton.h"

__host__ __device__
OptState2D Newton::doStep(TopParametric *sur, vec3d &P, OptState2D &loc) {
    double u = loc.u;
    double v = loc.v;
    Gu = sur->sqDistToDerU(u,v,P);
    Gv = sur->sqDistToDerV(u,v,P);

    Huu = sur->sqDistToDerUU(u, v, P);
    Huv = sur->sqDistToDerUV(u, v, P);
    Hvv = sur->sqDistToDerVV(u, v, P);

    double invConst = 1 / (Huu * Hvv - Huv * Huv);

    double Htuu, Htuv, Htvv;
    Htuu = Huu;
    Htuv = Huv;
    Htvv = Hvv;

    Huu = invConst * Htvv;
    Huv = invConst * (-Htuv);
    Hvv = invConst * Htuu;

    u -= sigma * (Huu * Gu + Huv * Gv);
    v -= sigma * (Huv * Gu + Hvv * Gv);

    return {u,v,sur->distTo(u,v,P)};
}

__host__ __device__
OptState2D Newton::doStep(BottomParametric *sur, vec3d &P, OptState2D &loc) {
    double u = loc.u;
    double v = loc.v;
    Gu = sur->sqDistToDerU(u,v,P);
    Gv = sur->sqDistToDerV(u,v,P);

    Huu = sur->sqDistToDerUU(u, v, P);
    Huv = sur->sqDistToDerUV(u, v, P);
    Hvv = sur->sqDistToDerVV(u, v, P);

    double invConst = 1 / (Huu * Hvv - Huv * Huv);

    double Htuu, Htuv, Htvv;
    Htuu = Huu;
    Htuv = Huv;
    Htvv = Hvv;

    Huu = invConst * Htvv;
    Huv = invConst * (-Htuv);
    Hvv = invConst * Htuu;

    u -= sigma * (Huu * Gu + Huv * Gv);
    v -= sigma * (Huv * Gu + Hvv * Gv);

    return {u,v,sur->distTo(u,v,P)};
}

__host__ __device__
OptState1D Newton::doStep(cubiccrv &crv, vec3d &P, OptState1D &loc) {

    double dt = -crv.sqDistToDer1(loc.t,P)/crv.sqDistToDer2(loc.t,P);

    loc.t += dt;

    loc.dist = crv.distTo(loc.t,P);
    return loc;
}
