#ifndef NLTVCSAD_MODEL_W_H
#define NLTVCSAD_MODEL_W_H

void  intialize_stuff_nltvcsad_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore, int w, int h);

void  free_stuff_nltvcsad_w(SpecificOFStuff *ofStuff);



void eval_nltvcsad_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff_W *nltvcsadw,
    float *ener_N,
    int ii, // initial column
    int ij, // initial row
    int ei, // end column
    int ej, // end row
    float lambda,  // weight of the data term
    float theta
    );

void guided_nltvcsad_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTvCsadStuff_W *nltvcsadw,
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
