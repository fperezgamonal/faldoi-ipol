#ifndef TVCSAD_MODEL_W_H
#define TVCSAD_MODEL_W_H


void  initialize_stuff_tvcsad_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore, int w, int h);


void  free_stuff_tvcsad_w(SpecificOFStuff *ofStuff);

void eval_tvcsad_w(
    const float *I0,
    const float *I1,
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
    float *ener_N,
    int ii,
    int ij,
    int ei,
    int ej,
    float lambda,
    float theta
    );


void guided_tvcsad_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
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
    int nx,        // width of I0 (and I1)
    int ny         // height of I0 (and I1)
  );

#endif //TVCSAD 
