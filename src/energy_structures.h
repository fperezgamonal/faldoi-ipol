#ifndef MATCH_VARIATIONAL_CORE_DATA_H
#define MATCH_VARIATIONAL_CORE_DATA_H

#include <vector>
#include "parameters.h"
#include <iostream>
#include <queue>

#define MAX(x,y) ((x)>(y)?(x):(y))


//Struct
struct SparseOF {
    int i; // column
    int j; // row
    float u; //x- optical flow component
    float v; //y- optical flow component
    float sim_node; //similarity measure for the actual pixel
    float occluded; //similarity measure for the accumulated path.
};

class CompareSparseOF {
public:
    bool operator()(SparseOF& e1, SparseOF& e2){
        return e1.sim_node > e2.sim_node;
    }
};

//Priority queue
typedef std::priority_queue<SparseOF, std::vector<SparseOF>, CompareSparseOF> pq_cand;

//Empty priority queue
template< typename T >
void makeEmpty( std::priority_queue<T>& q )
{
    std::priority_queue<T>empty;
    std::swap( q, empty );
}

struct BilateralWeight{
    float wp[NL_DUAL_VAR]; // weight of non local
    int   apj[NL_DUAL_VAR]; //absolute positon of p(y,x) (row)
    int   api[NL_DUAL_VAR]; //absolute position of p(y,x) (colum)
    float wt = 0.0;
};

struct PatchIndexes{
    int i; //center column
    int j; //center row
    int ii; // initial column
    int ij; // initial row
    int ei; // end column
    int ej; // end row
};

struct Weights_Bilateral{
    float* weight;
};

struct Parameters{
    float lambda;
    float theta;
    float tau;
    float beta;
    float alpha;
    float tau_u;
    float tau_eta;
    float tau_chi;
    float tol_OF;
    float mu;
    bool verbose;
    int warps;
    int w;
    int h;
    int pd;
    int w_radio;
    int val_method;
    int step_algorithm;
    int iterations_of;
    int max_iter_patch;
    int split_img;
    int h_parts;
    int v_parts;
	float epsilon;
	int part_res;
};

inline std::ostream& operator<<(std::ostream& os, const Parameters& p){
    return os << "Parameters: \n lambda: " << p.lambda << ", theta: " << p.theta << ", beta: " << p.beta
              << ", alpha: " << p.alpha << ", \n tau_u: " << p.tau_u << ", tau_eta: " << p.tau_eta
              << ", tau_chi: " << p.tau_chi << ", mu: " << p.mu <<  "\n";

}

struct OpticalFlowData{
    /* data */
    //TODO: This should be outside of this structure
    float * __restrict u1;
    float * __restrict u2;
    float * __restrict u1_ba;
    float * __restrict u2_ba;
    float * __restrict u1_filter;
    float * __restrict u2_filter;
    float * __restrict chi;
    int   * __restrict fixed_points;
    int   * __restrict trust_points;
    float * __restrict saliency; //It stores the saliency value for each pixel.

    Parameters params;
};

struct BilateralFilterData{
    Weights_Bilateral *weights_filtering;
    std::vector<PatchIndexes> indexes_filtering;
};

struct DualVariables{
    float sc[NL_DUAL_VAR]; // value of p(x,y)
    float wp[NL_DUAL_VAR]; // weight of non local
    int   apj[NL_DUAL_VAR]; //absolute positon of p(y,x) (row)
    int   api[NL_DUAL_VAR]; //absolute position of p(y,x) (colum)
    int   rp[NL_DUAL_VAR]; //relative position of p(y,x) in the structure
    float wt = 0.0;
};


//Struct
struct PosNei{
    int   api[DT_NEI]; //absolute positon of Intensity (row)
    int   apj[DT_NEI]; //absolute position of intensity (colum)
    float b[DT_NEI];
    std::vector<float>  ba;
    int n;
};

////Specific struct for the different functionals
struct  Tvl2CoupledOFStuff{
    //Dual variables
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;

    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;

    float *v1;
    float *v2;

    float *rho_c;
    float *grad;

    float *u1_;
    float *u2_;

    float *u1Aux;
    float *u2Aux;

    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;

    float *div_xi1;
    float *div_xi2;
    float *u_N;
};

struct NonLocalTVL1Stuff{
    DualVariables *p;
    DualVariables *q;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};


struct  TvCsadStuff{

    PosNei *pnei;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
};

struct NonLocalTvCsadStuff{

    DualVariables *p;
    DualVariables *q;
    PosNei *pnei;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};

///////////////////PESOS///////////////////////////////////////////////////////
struct  Tvl2CoupledOFStuff_W{

    int iiw;
    int ijw;
    float *weight;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1Aux;
    float *u2Aux;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
    float *u_N;
};

struct NonLocalTvCsadStuff_W{

    int iiw;
    int ijw;
    float *weight;
    DualVariables *p;
    DualVariables *q;
    PosNei *pnei;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};


struct NonLocalTVL1Stuff_W{

    int iiw;
    int ijw;
    float *weight;
    DualVariables *p;
    DualVariables *q;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};

struct  TvCsadStuff_W{

    int iiw;
    int ijw;
    float *weight;
    PosNei *pnei;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
};

////////////////////////
/////////OCCLUSIONS/////
///////////////////////
struct  Tvl2CoupledOFStuff_occ {
    //Occlusion variable
    float * __restrict chi;
    float * __restrict chix;
    float * __restrict chiy;

    //Weigth
    float * __restrict g;

    float * __restrict diff_u_N;

    //Dual variables
    float * __restrict xi11;
    float * __restrict xi12;
    float * __restrict xi21;
    float * __restrict xi22;

    float * __restrict u1x;
    float * __restrict u1y;
    float * __restrict u2x;
    float * __restrict u2y;
    //Dual variables for chi
    float * __restrict eta1;
    float * __restrict eta2;


    float * __restrict v1;
    float * __restrict v2;

    float * __restrict rho_c1;
    float * __restrict rho_c_1;
    float * __restrict grad_1;
    float * __restrict grad__1;

    float * __restrict I0x;
    float * __restrict I0y;

    float * __restrict I1x;
    float * __restrict I1y;
    float * __restrict I1w;
    float * __restrict I1wx;
    float * __restrict I1wy;

    float * __restrict I_1x;
    float * __restrict I_1y;
    float * __restrict I_1w;
    float * __restrict I_1wx;
    float * __restrict I_1wy;


    float * __restrict vi_div1;
    float * __restrict grad_x1;
    float * __restrict grad_y1;
    float * __restrict vi_div2;
    float * __restrict grad_x2;
    float * __restrict grad_y2;
    float * __restrict g_xi11;
    float * __restrict g_xi12;
    float * __restrict g_xi21;
    float * __restrict g_xi22;
    float * __restrict div_g_xi1;
    float * __restrict div_g_xi2;


    float * __restrict F;
    float * __restrict G;

    float * __restrict div_u;
    float * __restrict g_eta1;
    float * __restrict g_eta2;
    float * __restrict div_g_eta;
};



/////////////////////////////////////


//General Struct for the auxiliar stuff.
//Each pointer contains the auxiliar necessary information to estimate the OF.
struct  SpecificOFStuff{
    //TODO: Should think of a better option. To link the things
    //Creo que el problema viene por como declaramos los punteros y todo eso.

    Tvl2CoupledOFStuff  tvl2;
    NonLocalTVL1Stuff   nltvl1;
    TvCsadStuff         tvcsad;
    NonLocalTvCsadStuff nltvcsad;

    Tvl2CoupledOFStuff_W  tvl2w;
    NonLocalTvCsadStuff_W nltvcsadw;
    NonLocalTVL1Stuff_W   nltvl1w;
    TvCsadStuff_W         tvcsadw;

    Tvl2CoupledOFStuff_occ  tvl2_occ;

};


// Struct to contain all the partition-specific variables (generalisation of the 1 partition (whole image) case...)
struct PartitionData {
        int idx;                            // Partition idx: from 0 to NUM_PART - 1
        int width;                          // Partition idx's width
        int height;                         //     "       "   height
        int off_x;                          //     "       "   width offset
        int off_y;                          //     "       "   height  "
        float *i0;
        float *i1;
        float *i_1;
        float *i2;
        float *i0n;                         // Partition nº idx of normalized source image (time step 't')
        float *i1n;                         //     "        "    "     "      second   "   (  "    "  t+1)
        float *i_1n;                        //     "        "    "     "      previous "   (  "    "  t-1)
        float *i2n;                         //     "        "    "     "      third    "   (  "    "  t+2)
        float *oft0;                        // Forward 'trusted' optical flow
        float *oft1;                        // Backward 'trusted' optical flow
        float *ene_Go;                      // Forward energy
        float *ene_Ba;                      // Backward energy
        float *occ_Go;                      // Forward occlusion mask
        float *occ_Ba;                      // Backward occlusion mask
        float *sal_go;                      // Forward saliency (not used ATM but just in case)
        float *sal_ba;                      // Backward saliency ( "   "   "   "   "    "   ")
        pq_cand queue_Go;                  // Forward candidates' queue
        pq_cand queue_Ba;                  // Backward candidates' queue
        int queue_Go_size;
	    int queue_Ba_size;
	    SpecificOFStuff stuffGo;           // Specific stuff for each functional (forward)
        SpecificOFStuff stuffBa;           //     "      "    "    "      "     (backward)
        OpticalFlowData ofGo;               // Common OF data (forward)
        OpticalFlowData ofBa;               // Common OF data (backward)
        BilateralFilterData* BiFilt_Go;      // Bilateral filter data (forward)
        BilateralFilterData* BiFilt_Ba;      // Bilateral filter data (backward)
};

#endif// ENERGY_STRUCTURES_H
