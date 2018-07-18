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

void  initialize_stuff_nltvl1(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore, int w, int h);


void  free_stuff_nltvl1(SpecificOFStuff *ofStuff);

void eval_nltvl1(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        NonLocalTVL1Stuff *nltvl1,
        float *ener_N,
        PatchIndexes index,
        float lambda,  // weight of the data term
        float theta
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
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej, // end row
        int w,
        float *u1,
        float *u2,
        float *err
        );

void nltvl1_getD(
        float *u1,
        float *u2,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej, // end row
        int w,
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
        PatchIndexes index, // end row
        float lambda,  // weight of the data term
        float theta,   // weight of the data term
        float tau,     // time step
        float tol_OF,  // tol max allowed
        int   warps,   // number of warpings per scale
        bool  verbose, // enable/disable the verbose mode
        int w,         // width of I0 (and I1)
        int h          // height of I0 (and I1)
        );

#endif //TVL2-L1 functional
