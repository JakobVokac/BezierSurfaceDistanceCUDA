//
// Created by s3179222 on 12/4/19.
//
/*
 * This is a simple implementation of a vector with doubles for higher precision.
 * It overrides all of the basic operators and adds the dot product, cross product, an output function for cout.
 *
 * There is also an implementation for solving a system of 3 linear equations, which is useful for finding the
 * parameters in a linear combination of vectors for another vector.
 */

#ifndef HEARTVALVEMODEL_VEC3D_CUH
#define HEARTVALVEMODEL_VEC3D_CUH


class vec3d {
private:
public:
    double x,y,z;

    vec3d() = default;
    __host__ __device__ vec3d(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    __host__ __device__ vec3d& operator += (const vec3d& b);

    __host__ __device__ vec3d& operator -= (const vec3d& b);

    __host__ __device__ vec3d& operator *= (const vec3d& b);

    __host__ __device__ vec3d& operator /= (vec3d b);

    __host__ __device__ vec3d& operator += (const double& b);

    __host__ __device__ vec3d& operator -= (const double& b);

    __host__ __device__ vec3d& operator *= (const double& b);

    __host__ __device__ vec3d& operator /= (const double& b);
    __host__ __device__ vec3d operator + ();
    __host__ __device__ vec3d operator - ();
    __host__ __device__ double dot (vec3d b);
    __host__ __device__ vec3d cross (vec3d b);
    __host__ __device__ double mag ();
    __host__ __device__ double sum ();
    __host__ __device__ double sqdist (vec3d b);
    __host__ __device__ double dist (vec3d b);

    __host__ __device__ vec3d RotAxisZ(double rad);
    __host__ __device__ vec3d MirrorY();

    __host__ __device__ double getx();
    __host__ __device__ double gety();
    __host__ __device__ double getz();

    __host__ __device__ bool operator == (const vec3d& b) const;
    __host__ __device__ bool operator != (const vec3d& b) const;
};

__host__ __device__ vec3d operator + (vec3d a, const vec3d& b);
__host__ __device__ vec3d operator - (vec3d a, const vec3d& b);
__host__ __device__ vec3d operator * (vec3d a, const vec3d& b);
__host__ __device__ vec3d operator / (vec3d a, const vec3d& b);
__host__ __device__ vec3d operator + (vec3d a, const double& b);
__host__ __device__ vec3d operator - (vec3d a, const double& b);
__host__ __device__ vec3d operator * (vec3d a, const double& b);
__host__ __device__ vec3d operator / (vec3d a, const double& b);
__host__ __device__ vec3d operator + (const double& b, vec3d a);
__host__ __device__ vec3d operator - (const double& b, vec3d a);
__host__ __device__ vec3d operator * (const double& b, vec3d a);
__host__ __device__ vec3d operator / (const double& b, vec3d a);

__host__ __device__ double curvature(vec3d c1, vec3d c2);
__host__ __device__ vec3d circleCenterDir(vec3d c1,vec3d c2);
__host__ __device__ double sign(vec3d P, vec3d Q);
__host__ __device__ vec3d solve3Dlinear(vec3d v1, vec3d v2, vec3d v3, vec3d v4);


#endif //HEARTVALVEMODEL_VEC3D_CUH
