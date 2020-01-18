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
#include "optimizer.h"

__host__ __device__ double myabs(double a){
	if(a < 0)
		return -a;
	return a;
}

__host__ __device__
OptState2D optimizer::optimize() {
    initialize();
    if(useTopOrBot)
    	tLoc = loc = prep2D.preprocess(TopSur, P);
    else
    	tLoc = loc = prep2D.preprocess(BottomSur, P);
    if(isnan(loc.u) || isnan(loc.v))
    	printf("Preprocessed u: %lf, v %lf\n",loc.u, loc.v);
    do{
        distOld = loc.dist;
        if(useTopOrBot)
        	loc = step.doStep(TopSur,P,loc);
        else
        	loc = step.doStep(BottomSur,P,loc);

        if(!edgesSearched) {
            checkEdges();
        }
        iter++;
        if(isnan(loc.u) || isnan(loc.v))
        	printf("iter u: %lf, v %lf\n",loc.u, loc.v);

    }while(myabs(distOld-loc.dist) > eps && iter < iterMax);
    if(isnan(tLoc.u) || isnan(tLoc.v))
		printf("2ter u: %lf, v %lf, %d, %d\n",tLoc.u, tLoc.v, cornerSearched, edgesSearched);
    if(loc.u < 0 || loc.v < 0 || loc.u > 1 || loc.v > 1)
        return tLoc;
    if(isnan(tLoc.u) || isnan(tLoc.v))
		printf("3ter u: %lf, v %lf\n",tLoc.u, tLoc.v);
    if(tLoc.dist < loc.dist)
        return tLoc;
    if(isnan(loc.u) || isnan(loc.v))
		printf("4ter u: %lf, v %lf\n",loc.u, loc.v);

    return loc;
}

__host__ __device__
OptState1D optimizer::optimizeEdge(cubiccrv &crv) {
    eLoc = prep1D.preprocess(crv, P);
    eLoc.dist = crv.distTo(eLoc.t,P);
    do{
        distOld = eLoc.dist;
        eLoc = step.doStep(crv, P, eLoc);
        checkCorner();
        iter++;
    }while(myabs(distOld - eLoc.dist) > eps && iter < iterMax);

    OptState1D lc, rc;
    lc = {0,crv.distTo(0,P)};
    rc = {1,crv.distTo(1,P)};
    if(eLoc.t < 0 || eLoc.t > 1)
        return (lc.dist < rc.dist ? lc : rc);

    if(lc.dist < eLoc.dist)
        eLoc = lc;
    if(rc.dist < eLoc.dist)
        eLoc = rc;

    return eLoc;
}

__host__ __device__
void optimizer::checkEdges() {

    if (loc.u < 0 && rangecount > rcMax) {

        if (loc.v < 0) {
            cornerSearch(0, 0);
        }else if (loc.v > 1) {
            cornerSearch(0, 1);
        }else {
        	OptState1D res;
        	if(useTopOrBot)
        		res = optimizeEdge(TopSur->edgeU0());
        	else
                res = optimizeEdge(BottomSur->edgeU0());

            tLoc = {0, res.t, res.dist};

//            plotEdge(TopSur.edgeU0());
        }

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.u > 1 && rangecount > rcMax) {

        if(loc.v < 0) {
            cornerSearch(1, 0);
        }else if (loc.v > 1) {
            cornerSearch(1, 1);
        }else {
        	OptState1D res;
        	if(useTopOrBot)
        		res = optimizeEdge(TopSur->edgeU1());
        	else
                res = optimizeEdge(BottomSur->edgeU1());

            tLoc = {1, res.t, res.dist};

//            plotEdge(TopSur.edgeU1());
        }
        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.v < 0 && rangecount > rcMax) {

    	OptState1D res;
    	if(useTopOrBot)
    		res = optimizeEdge(TopSur->edgeV0());
    	else
            res = optimizeEdge(BottomSur->edgeV0());

        tLoc = {res.t, 0, res.dist};

//        plotEdge(TopSur.edgeV0());

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.v > 1 && rangecount > rcMax) {

    	OptState1D res;
    	if(useTopOrBot)
    		res = optimizeEdge(TopSur->edgeV1());
    	else{
			vec3d Q_r = BottomSur->getBottomCorner();
			res = {0.5,Q_r.dist(P)};
		}
    	tLoc = {res.t, 1, res.dist};

//        plotEdge(TopSur.edgeV1());

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;
}

__host__ __device__
void optimizer::cornerSearch(double u, double v) {
    OptState1D res1{}, res2{};

    if(useTopOrBot){
		if(u == 0)
			res1 = optimizeEdge(TopSur->edgeU0());
		else
			res1 = optimizeEdge(TopSur->edgeU1());
		if(v == 0)
			res2 = optimizeEdge(TopSur->edgeV0());
		else
			res2 = optimizeEdge(TopSur->edgeV1());
    }else{
    	if(u == 0)
			res1 = optimizeEdge(BottomSur->edgeU0());
		else
			res1 = optimizeEdge(BottomSur->edgeU1());
		if(v == 0)
			res2 = optimizeEdge(BottomSur->edgeV0());
		else{
			vec3d Q_r = BottomSur->getBottomCorner();
			res2 = {0.5,Q_r.dist(P)};
		}
    }

//    plotEdge(TopSur.edgeU0());
//    plotEdge(TopSur.edgeV0());

    if (res1.dist < res2.dist)
        tLoc = {res1.t, v, res1.dist};
    else
        tLoc = {u, res2.t, res2.dist};

    edgeLoc = true;
    cornerSearched = true;
}

__host__ __device__
void optimizer::checkCorner() {
    if(eLoc.t < 0){
        eLoc.t = 0;
        cornerLoc = true;
    }else if(eLoc.t > 1){
        eLoc.t = 1;
        cornerLoc = true;
    }
}

__host__ __device__
void optimizer::initialize() {
    edgeLoc = false;
    cornerLoc = false;
    cornerSearched = false;
    edgesSearched = false;
    iter = 0;
    rangecount = 0;
}

__host__ __device__
OptState2D optimizer::optimizeForPoint(vec3d P) {
    this->P = P;

    return optimize();
}

__host__ __device__
TopParametric &optimizer::getSurface() {
    return *TopSur;
}

__host__ __device__
bool optimizer::edgeSolution() {
    return edgeLoc;
}

__host__ __device__
bool optimizer::cornerSearchSolution() {
    return cornerSearched;
}

__host__ __device__
bool optimizer::cornerSolution() {
    return cornerLoc;
}

__host__ __device__
int optimizer::getIterations() {
    return iter;
}
