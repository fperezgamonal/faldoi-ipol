#ifndef TVL2_MODEL_OCC_H
#define TVL2_MODEL_OCC_H

#include "energy_structures.h"
#include "aux_energy_model.h"


//OPTICAL FLOW PARAMETERS
#define TVL2_LAMBDA  40//40
#define TVL2_THETA   0.3
#define TVL2_TAU     0.125 //0.25
#define TVL2_NWARPS  1  //5
#define TVL2_TOL_D   0.01
#define TVL2_VERBOSE 0  //0

////INITIALIZATION OF EACH METHOD
void intialize_stuff_tvl2coupled_occ(
        SpecificOFStuff& ofStuff,
        const OpticalFlowData& ofCore);


void free_stuff_tvl2coupled_occ(SpecificOFStuff *ofStuff);

void eval_tvl2coupled_occ(
        const float *I0,           // source image
        const float *I1,           // forward image
        const float *I_1,           // backward image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
        float *ener_N,
        const PatchIndexes index, // end row
        const float lambda,  // weight of the data term
        const float theta,
        const float alpha,
        const float beta
        );

// Variational Optical flow method based on initial fixed values
// It minimizes the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)|
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)


void guided_tvl2coupled_occ(
        const float *I1,           // source image
        const float *I2,           // forward image
        const float *I0,           // backward image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
        float *ener_N,
        const PatchIndexes index
        /*,
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau_u,     // time step for u
        const float tau_eta,     // time step for eta
        const float tau_chi,     // time step for chi
        const float beta,
        const float alpha, //weight of the norm term
        const float tol_OF,  // tol max allowed
        const int   warps,   // number of warpings
        const bool  verbose  // enable/disable the verbose mode*/
        );

void tvl2OF_occ(
        const float *I0,           // source image
        const float *I1,           // target image
        const float *I_1,
        float *u1,           // x component of the optical flow
        float *u2,           // y component of the optical flow
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        float *chi,
        const Parameters params
        );
#endif //TVL2-L1 functional with occlusions
