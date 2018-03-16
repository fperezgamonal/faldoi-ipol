#ifndef NLTVL1_MODEL
#define NLTVL1_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"

//OPTICAL FLOW PARAMETERS
#define NLTV_LAMBDA  40//40
#define NLTV_THETA   0.3
#define NLTV_TAU     0.125 //0.25
#define NLTV_NWARPS  1  //5
#define NLTV_TOL_D   0.01
#define NLTV_VERBOSE 0  //0

void  intialize_stuff_nltvl1(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore);


void  free_stuff_nltvl1(SpecificOFStuff *ofStuff);

void eval_nltvl1(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        NonLocalTVL1Stuff *nltvl1,
        float *ener_N,
        const PatchIndexes index,
        const float lambda,  // weight of the data term
        const float theta
        );



/*
 * - Name: getP

 *
*/

void nltvl1_getP(
        float *v1,
        float *v2,
        float *div_p1,
        float *div_p2,
        int *mask,
        float theta,
        float tau,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int w,
        float *u1,
        float *u2,
        float *err
        );

void nltvl1_getD(
        float *u1,
        float *u2,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int w,
        int n_d,
        float tau,
        DualVariables *p1,
        DualVariables *p2
        );

void guided_nltvl1(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        NonLocalTVL1Stuff *nltvl1,
        float *ener_N,
        const PatchIndexes index, // end row
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   warps,   // number of warpings per scale
        const bool  verbose  // enable/disable the verbose mode
        );

#endif //TVL2-L1 functional
