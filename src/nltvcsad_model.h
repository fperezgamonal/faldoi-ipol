#ifndef NLTVCSAD_MODEL_H
#define NLTVCSAD_MODEL_H

#include "energy_structures.h"
#include "energy_model.h"

//OPTICAL FLOW PARAMETERS
#define NLTVCSAD_LAMBDA  40//40
#define NLTVCSAD_THETA   0.3
#define NLTVCSAD_TAU     0.125 //0.25
#define NLTVCSAD_NWARPS  1  //5
#define NLTVCSAD_TOL_D   0.01
#define NLTVCSAD_VERBOSE 0  //0

void  initialize_stuff_nltvcsad(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore, int w, int h);

void  free_stuff_nltvcsad(SpecificOFStuff *ofStuff);




void eval_nltvcsad(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff *nltvcsad,
    float *ener_N,
    int ii, // initial column
    int ij, // initial row
    int ei, // end column
    int ej, // end row
    float lambda,  // weight of the data term
    float theta
    );



void guided_nltvcsad(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff *nltvcsad,
    float *ener_N,
    int ii, // initial column
    int ij, // initial row
    int ei, // end column
    int ej, // end row
    float lambda,  // weight of the data term
    float theta,   // weight of the data term
    float tau,     // time step
    float tol_OF,  // tol max allowed
    int   warps,   // number of warpings per scale
    bool  verbose, // enable/disable the verbose mode
    int w,         // width of I0 (and I1)
    int h          // height of I0 (and I1)
    );
#endif
