#ifndef TVL2W_MODEL_H
#define TVL2W_MODEL_H


////INITIALIZATION OF EACH METHOD
void  initialize_stuff_tvl2coupled_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore, int w, int h);

void  free_stuff_tvl2coupled_w(SpecificOFStuff *ofStuff);


void eval_tvl2coupled_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    Tvl2CoupledOFStuff_W *tvl2w,
    float *ener_N,
    int ii, // initial column
    int ij, // initial row
    int ei, // end column
    int ej, // end row
    float lambda,  // weight of the data term
    float theta
    );

void guided_tvl2coupled_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    Tvl2CoupledOFStuff_W *tvl2w,
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
    int nx,               // width of I0 (and I1)
    int ny                // height of I0 (and I1)
    );
#endif //TVL2-L1 functional
