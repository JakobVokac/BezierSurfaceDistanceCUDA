//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the main numerical optimizer class. It implements the skeleton of the numerical optimization, which includes:
 *  - calling the preprocessing method for the surface search (prep2D)
 *  - iterating through the surface search with a specified step (step2D),
 *  - checking if the optimization has been outsize of the domain ((u,v) e [0,1]^2) for 2 or more iterations (rangecount)
 *    in which case an edge or corner search is called,
 *  - calling the preprocessing method for the edge search (prep1D)
 *  - iterating through the edge with a specified step (step1D)
 *  - choosing the most appropirate solution between the surface and edge searches,
 *  - keeping track of the number of iterations
 */
#ifndef HEARTVALVEMODEL_OPTIMIZER_H
#define HEARTVALVEMODEL_OPTIMIZER_H


#include "../geometry/surface/TopParametric.h"
#include "../geometry/surface/BottomParametric.h"
#include "states.h"
#include "preprocessor/bisection.h"
#include "preprocessor/quadraticInterpolation.h"
#include "step/Newton.h"
#include <cstdlib>

class optimizer {
protected:
    bisection prep2D;
    quadraticInterpolation prep1D;
    Newton step;
    TopParametric *TopSur;
    BottomParametric *BottomSur;
    vec3d P;
    double eps;
    double distOld = 0;
    int iter = 0;
    int iterMax;
    int rangecount = 0;
    int rcMax;
    OptState2D loc{};
    OptState1D eLoc{};
    OptState2D tLoc{};
    bool edgeLoc = false;
    bool cornerLoc = false;
    bool cornerSearched = false;
    bool edgesSearched = false;
    bool useTopOrBot = true;
public:
    __host__ __device__ optimizer(
    		bisection prep2D,
			quadraticInterpolation prep1D,
			Newton step,
			TopParametric *TopSur,
			BottomParametric *BottomSur,
            vec3d P,
            double eps,
            int iterMax,
            int outOfRangeMax
    )
            : prep2D(prep2D), prep1D(prep1D), step(step){
    	this->TopSur = TopSur;
    	this->BottomSur = BottomSur;
        this->eps = eps;
        this->iterMax = iterMax;
        this->rcMax = outOfRangeMax;
        this->P = P;
        this->loc = {0.5, 0.5, distOld};
    };

    __host__ __device__ OptState2D optimize();

    __host__ __device__ OptState2D optimizeForPoint(vec3d P);

    __host__ __device__ OptState1D optimizeEdge(cubiccrv &crv);

    ~optimizer() = default;

    __host__ __device__ void checkEdges();

    __host__ __device__ void cornerSearch(double u, double v);

    __host__ __device__ void checkCorner();

    __host__ __device__ void initialize();

    __host__ __device__ TopParametric &getSurface();

    __host__ __device__ int getIterations();

    __host__ __device__ bool edgeSolution();

    __host__ __device__ bool cornerSolution();

    __host__ __device__ bool cornerSearchSolution();

    __host__ __device__ void setTopOrBot(bool b){
    	this->useTopOrBot = b;
    }
    __host__ __device__ void setTop(TopParametric *s){
    	this->TopSur = s;
    }
    __host__ __device__ void setBot(BottomParametric *s){
    	this->BottomSur = s;
    }
};

#endif //HEARTVALVEMODEL_OPTIMIZER_H
