//
// Created by s3179222 on 12/8/19.
//

#ifndef HEARTVALVEMODEL_MYSURFACE_H
#define HEARTVALVEMODEL_MYSURFACE_H

#include "../curve/cubiccrv.h"

class mySurface {
private:
	cubiccrv emptycrv = cubiccrv();
public:
    __host__ __device__ vec3d at(double u, double v){return vec3d(0,0,0);};
    __host__ __device__ vec3d atDerU(double u, double v){return vec3d(0,0,0);};
    __host__ __device__ vec3d atDerV(double u, double v){return vec3d(0,0,0);};
    __host__ __device__ vec3d atDerUU(double u, double v){return vec3d(0,0,0);};
    __host__ __device__ vec3d atDerVV(double u, double v){return vec3d(0,0,0);};
    __host__ __device__ vec3d atDerUV(double u, double v){return vec3d(0,0,0);};

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

    __host__ __device__ cubiccrv & edgeU0(){
    	return emptycrv;
    };
    __host__ __device__ cubiccrv & edgeU1(){
    	return emptycrv;
    };
    __host__ __device__ cubiccrv & edgeV0(){
    	return emptycrv;
    };
    __host__ __device__ cubiccrv & edgeV1(){
    	return emptycrv;
    };
protected:
private:
};


#endif //HEARTVALVEMODEL_MYSURFACE_H
