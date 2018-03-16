#ifndef NLTVCSAD_MODEL_H
#define NLTVCSAD_MODEL_H

//OPTICAL FLOW PARAMETERS
#define NLTVCSAD_LAMBDA  40//40
#define NLTVCSAD_THETA   0.3
#define NLTVCSAD_TAU     0.125 //0.25
#define NLTVCSAD_NWARPS  1  //5
#define NLTVCSAD_TOL_D   0.01
#define NLTVCSAD_VERBOSE 0  //0

void  intialize_stuff_nltvcsad(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore);

void  free_stuff_nltvcsad(SpecificOFStuff *ofStuff);




void eval_nltvcsad(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff *nltvcsad,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta
    );



void guided_nltvcsad(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff *nltvcsad,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta,   // weight of the data term
    const float tau,     // time step
    const float tol_OF,  // tol max allowed
    const int   warps,   // number of warpings per scale
    const bool  verbose  // enable/disable the verbose mode
    );
#endif
