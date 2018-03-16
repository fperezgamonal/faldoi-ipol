#ifndef TVCSAD_MODEL_W_H
#define TVCSAD_MODEL_W_H


void  intialize_stuff_tvcsad_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore);


void  free_stuff_tvcsad_w(SpecificOFStuff *ofStuff);

void eval_tvcsad_w(
    const float *I0,
    const float *I1,
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
    float *ener_N,
    const int ii,
    const int ij,
    const int ei,
    const int ej,
    const float lambda,
    const float theta
    );


void guided_tvcsad_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
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

#endif //TVCSAD 
