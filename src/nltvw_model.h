#ifndef NLTVL1_MODEL_W_H
#define NLTVL1_MODEL_W_H



void  initialize_stuff_nltvl1_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore, int w, int h);



void  free_stuff_nltvl1_w(SpecificOFStuff *ofStuff);


void eval_nltvl1_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
    float *ener_N,
    int ii, // initial column
    int ij, // initial row
    int ei, // end column
    int ej, // end row
    float lambda,  // weight of the data term
    float theta
    );


void guided_nltvl1_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
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

#endif //NLTVL1_MODEL_W
