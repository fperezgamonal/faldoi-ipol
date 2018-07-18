// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license athis program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// All rights reserved.


#ifndef TVL2_MODEL_OCC
#define TVL2_MODEL_OCC



#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "parameters.h"

extern "C" {
#include "bicubic_interpolation.h"
}
#include "utils.h"
#include <iostream>

////INITIALIZATION OF EACH METHOD
void  initialize_stuff_tvl2coupled_occ(
        SpecificOFStuff& ofStuff,
        const OpticalFlowData& ofCore,
        const int w,
        const int h)
{
    // Occlusion variable
    ofStuff.tvl2_occ.chix = new float[w*h];
    ofStuff.tvl2_occ.chiy = new float[w*h];

    ofStuff.tvl2_occ.diff_u_N = new float[w*h];
    // Weight
    ofStuff.tvl2_occ.g = new float[w*h];

    // Dual variables
    ofStuff.tvl2_occ.xi11 = new float[w*h];
    ofStuff.tvl2_occ.xi12 = new float[w*h];
    ofStuff.tvl2_occ.xi21 = new float[w*h];
    ofStuff.tvl2_occ.xi22 = new float[w*h];

    ofStuff.tvl2_occ.u1x = new float[w*h];
    ofStuff.tvl2_occ.u1y = new float[w*h];
    ofStuff.tvl2_occ.u2x = new float[w*h];
    ofStuff.tvl2_occ.u2y = new float[w*h];

    ofStuff.tvl2_occ.v1 = new float[w*h];
    ofStuff.tvl2_occ.v2 = new float[w*h];

    ofStuff.tvl2_occ.rho_c1 = new float[w*h];
    ofStuff.tvl2_occ.rho_c_1 = new float[w*h];
    ofStuff.tvl2_occ.grad_1 = new float[w*h];
    ofStuff.tvl2_occ.grad__1 = new float[w*h];

    if (ofCore.params.step_algorithm == GLOBAL_STEP) {
        ofStuff.tvl2_occ.I0x = new float[w*h];
        ofStuff.tvl2_occ.I0y = new float[w*h];
    } else {
        ofStuff.tvl2_occ.I0x = nullptr;
        ofStuff.tvl2_occ.I0y = nullptr;
    }

    ofStuff.tvl2_occ.I1x = new float[w*h];
    ofStuff.tvl2_occ.I1y = new float[w*h];
    ofStuff.tvl2_occ.I1w = new float[w*h];
    ofStuff.tvl2_occ.I1wx = new float[w*h];
    ofStuff.tvl2_occ.I1wy = new float[w*h];

    ofStuff.tvl2_occ.I_1x = new float[w*h];
    ofStuff.tvl2_occ.I_1y = new float[w*h];
    ofStuff.tvl2_occ.I_1w = new float[w*h];
    ofStuff.tvl2_occ.I_1wx = new float[w*h];
    ofStuff.tvl2_occ.I_1wy = new float[w*h];

    ofStuff.tvl2_occ.vi_div1 = new float[w*h];
    ofStuff.tvl2_occ.grad_x1 = new float[w*h];
    ofStuff.tvl2_occ.grad_y1 = new float[w*h];
    ofStuff.tvl2_occ.vi_div2 = new float[w*h];
    ofStuff.tvl2_occ.grad_x2 = new float[w*h];
    ofStuff.tvl2_occ.grad_y2 = new float[w*h];
    ofStuff.tvl2_occ.g_xi11 = new float[w*h];
    ofStuff.tvl2_occ.g_xi12 = new float[w*h];
    ofStuff.tvl2_occ.g_xi21 = new float[w*h];
    ofStuff.tvl2_occ.g_xi22 = new float[w*h];
    ofStuff.tvl2_occ.div_g_xi1 = new float[w*h];
    ofStuff.tvl2_occ.div_g_xi2 = new float[w*h];
    ofStuff.tvl2_occ.eta1 = new float[w*h];
    ofStuff.tvl2_occ.eta2 = new float[w*h];
    ofStuff.tvl2_occ.F = new float[w*h];
    ofStuff.tvl2_occ.G = new float[w*h];

    ofStuff.tvl2_occ.div_u = new float[w*h];
    ofStuff.tvl2_occ.g_eta1 = new float[w*h];
    ofStuff.tvl2_occ.g_eta2 = new float[w*h];
    ofStuff.tvl2_occ.div_g_eta = new float[w*h];
}

void  free_stuff_tvl2coupled_occ(SpecificOFStuff *ofStuff){

    delete [] ofStuff->tvl2_occ.chix;
    delete [] ofStuff->tvl2_occ.chiy;

    delete [] ofStuff->tvl2_occ.g;

    delete [] ofStuff->tvl2_occ.diff_u_N;

    delete [] ofStuff->tvl2_occ.xi11;
    delete [] ofStuff->tvl2_occ.xi12;
    delete [] ofStuff->tvl2_occ.xi21;
    delete [] ofStuff->tvl2_occ.xi22;

    delete [] ofStuff->tvl2_occ.v1;
    delete [] ofStuff->tvl2_occ.v2;

    delete [] ofStuff->tvl2_occ.rho_c1;
    delete [] ofStuff->tvl2_occ.rho_c_1;
    delete [] ofStuff->tvl2_occ.grad_1;
    delete [] ofStuff->tvl2_occ.grad__1;

    delete [] ofStuff->tvl2_occ.vi_div1;
    delete [] ofStuff->tvl2_occ.grad_x1;
    delete [] ofStuff->tvl2_occ.grad_y1;
    delete [] ofStuff->tvl2_occ.vi_div2;
    delete [] ofStuff->tvl2_occ.grad_x2;
    delete [] ofStuff->tvl2_occ.grad_y2;

    delete [] ofStuff->tvl2_occ.g_xi11;
    delete [] ofStuff->tvl2_occ.g_xi12;
    delete [] ofStuff->tvl2_occ.g_xi21;
    delete [] ofStuff->tvl2_occ.g_xi22;
    delete [] ofStuff->tvl2_occ.div_g_xi1;
    delete [] ofStuff->tvl2_occ.div_g_xi2;

    delete [] ofStuff->tvl2_occ.eta1;
    delete [] ofStuff->tvl2_occ.eta2;
    delete [] ofStuff->tvl2_occ.F;
    delete [] ofStuff->tvl2_occ.G;

    delete [] ofStuff->tvl2_occ.div_u;
    delete [] ofStuff->tvl2_occ.g_eta1;
    delete [] ofStuff->tvl2_occ.g_eta2;
    delete [] ofStuff->tvl2_occ.div_g_eta;

    delete [] ofStuff->tvl2_occ.I1x;
    delete [] ofStuff->tvl2_occ.I1y;
    delete [] ofStuff->tvl2_occ.I1w;
    delete [] ofStuff->tvl2_occ.I1wx;
    delete [] ofStuff->tvl2_occ.I1wy;

    if (ofStuff->tvl2_occ.I0x){
        delete [] ofStuff->tvl2_occ.I0x;
        delete [] ofStuff->tvl2_occ.I0y;
    }

    delete [] ofStuff->tvl2_occ.I_1x;
    delete [] ofStuff->tvl2_occ.I_1y;
    delete [] ofStuff->tvl2_occ.I_1w;
    delete [] ofStuff->tvl2_occ.I_1wx;
    delete [] ofStuff->tvl2_occ.I_1wy;

}

//////////////////////////////////////////////////////////////
////TV-l2 COUPLED OPTICAL FLOW PROBLEM WITH OCCLUSIONS////////
/////////////////////////////////////////////////////////////



float eval_tvl2coupled_occ(
        const float *I0,           // source image
        const float *I1,           // forward image
        const float *I_1,           // backward image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
        const PatchIndexes index,
        Parameters params,
        const int nx,
        const int ny
        ){

    const float *u1 = ofD->u1;
    const float *u2 = ofD->u2;
    float *u1_ba = ofD->u1_ba;
    float *u2_ba = ofD->u2_ba;

    //Occlusion variables
    const float *chi = ofD->chi;
    float *chix = tvl2_occ->chix;
    float *chiy = tvl2_occ->chiy;

    // Added changes for subimages
    //Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;

    //Optical flow derivatives
    const float *v1   = tvl2_occ->v1;
    const float *v2   = tvl2_occ->v2;
    float *u1x  = tvl2_occ->u1x;
    float *u1y  = tvl2_occ->u1y;
    float *u2x  = tvl2_occ->u2x;
    float *u2y  = tvl2_occ->u2y;

    //Weigth
    const float *g = tvl2_occ->g;

    float *div_u = tvl2_occ->div_u;

    float *I_1w = tvl2_occ->I_1w;
    float *I1w = tvl2_occ->I1w;

    //Derivatives and warping of I2
    const float *I1x = tvl2_occ->I1x;
    const float *I1y = tvl2_occ->I1y;
    float *I1wx = tvl2_occ->I1wx;
    float *I1wy = tvl2_occ->I1wy;

    //Derivatives and warping of I0
    const float *I_1x = tvl2_occ->I_1x;
    const float *I_1y = tvl2_occ->I_1y;
    float *I_1wx = tvl2_occ->I_1wx;
    float *I_1wy = tvl2_occ->I_1wy;


    float ener = 0.0;

    forward_gradient_patch(u1, u1x, u1y, index.ii, index.ij, index.ei, index.ej, nx);
    forward_gradient_patch(u2, u2x, u2y, index.ii, index.ij, index.ei, index.ej, nx);
    forward_gradient_patch(chi, chix, chiy, index.ii, index.ij, index.ei, index.ej, nx);
    divergence_patch(u1, u2, div_u, index.ii, index.ij, index.ei, index.ej, nx);

    //#pragma omp simd collapse(2)
    //#pragma omp for schedule(dynamic, 1) collapse(2)
    for (int l = index.ij; l < index.ej; l++){
        for (int k = index.ii; k < index.ei; k++){
            const int  i = l*nx + k;
            u1_ba[i] = -u1[i];
            u2_ba[i] = -u2[i];
        }
    }
    // Compute the warping of I1 and its derivatives I1(x + u1o), I1x (x + u1o) and I1y (x + u2o)
    bicubic_interpolation_warp_patch(I1,  u1, u2, I1w, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
    bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
    bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy, index.ii, index.ij, index.ei, index.ej, nx, ny, false);

    // Compute the warping of I0 and its derivatives I0(x - u1o), I0x (x - u1o) and I0y (x - u2o)
    bicubic_interpolation_warp_patch(I_1,  u1_ba, u2_ba, I_1w, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
    bicubic_interpolation_warp_patch(I_1x, u1_ba, u2_ba, I_1wx, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
    bicubic_interpolation_warp_patch(I_1y, u1_ba, u2_ba, I_1wy, index.ii, index.ij, index.ei, index.ej, nx, ny, false);

    //Energy for all the patch. Maybe it would be useful only the 8 pixels around the seed.
    int m  = 0;

    //#pragma omp simd collapse(2) reduction(+:ener)
    //#pragma omp for schedule(dynamic, 1) collapse(2)
    for (int l = index.ij; l < index.ej; l++){
        for (int k = index.ii; k < index.ei; k++){
            const int i = l*nx + k;

            const float diff_uv_term = (1/(2*params.theta))*
                    ((u1[i] - v1[i])*(u1[i]- v1[i]) + (u2[i] - v2[i])*(u2[i] - v2[i]));
            const float norm_v_term = (params.alpha/2)*chi[i]*(v1[i]*v1[i] + v2[i]*v2[i]);

            const float div_u_term = params.beta*chi[i]*div_u[i];

            const float rho_1 = fabs(I1w[i] - I1wx[i] * u1[i] - I1wy[i] * u2[i] - I0[i]
                                     + I1wx[i] * v1[i] + I1wy[i] * v2[i]);
            const float rho__1 = fabs(I_1w[i] - I_1wx[i] * u1[i] - I_1wy[i] * u2[i] - I0[i]
                                      + I_1wx[i] * v1[i] + I_1wy[i] * v2[i]);

            const float data_term = params.lambda * ((1 - chi[i])*rho_1 + chi[i]*rho__1);


            const float grad_u1 = sqrt(u1x[i] * u1x[i] + u1y[i] * u1y[i]);
            const float grad_u2 = sqrt(u2x[i] * u2x[i] + u2y[i] * u2y[i]);
            const float grad_chi = sqrt(chix[i] * chix[i] + chiy[i] * chiy[i]);


            const float smooth_term = g[i]*(grad_u1 + grad_u2 + params.mu*grad_chi);


            if (!std::isfinite(data_term)){
                std::printf("Corrupt data\n");
            }
            if (!std::isfinite(smooth_term)){
                std::printf("Corrupt regularization\n");
            }

            ener += data_term + smooth_term + div_u_term + norm_v_term + diff_uv_term;
            m++;
        }
    }

    ener /= (m*1.0);
    return ener;
}

/////////////////////////////////////
//////////// Minimization ///////////
/////////////////////////////////////


//Dual variables for u
static void tvl2coupled_get_xi_patch(
        float *xi11, //Dual variable
        float *xi12, //Dual variable
        float *xi21, //Dual variable
        float *xi22, //Dual variable
        const float *g,
        const float *v1,
        const float *v2,
        const float *chix,
        const float *chiy,
        float *vi_div1,
        float *grad_x1,
        float *grad_y1,
        float *vi_div2,
        float *grad_x2,
        float *grad_y2,
        float *g_xi11,
        float *g_xi12,
        float *g_xi21,
        float *g_xi22,
        float *div_g_xi1,
        float *div_g_xi2,
        const PatchIndexes index,
        const Parameters params,
        const int nx
        ){
    float tau_theta = params.tau_u/params.theta;
    //int nx = params.w;
    for (int k = 1; k < ITER_XI; k++){
        //What goes inside gradient
        for (int l = index.ij; l < index.ej; l++){
            for (int j = index.ii; j < index.ei; j++){
                const int i = l*nx + j;
                g_xi11[i] = g[i]*xi11[i];
                g_xi12[i] = g[i]*xi12[i];

                g_xi21[i] = g[i]*xi21[i];
                g_xi22[i] = g[i]*xi22[i];
            }
        }

        divergence_patch(g_xi11, g_xi12, div_g_xi1, index.ii, index.ij, index.ei, index.ej, nx);
        divergence_patch(g_xi21, g_xi22, div_g_xi2, index.ii, index.ij, index.ei, index.ej, nx);


        for (int l = index.ij; l < index.ej; l++){
            for (int j = index.ii; j < index.ei; j++){
                const int i = l*nx + j;
                vi_div1[i] = v1[i] + params.theta*div_g_xi1[i] + params.theta*params.beta*chix[i];
                vi_div2[i] = v2[i] + params.theta*div_g_xi2[i] + params.theta*params.beta*chiy[i];
            }
        }


        forward_gradient_patch(vi_div1, grad_x1, grad_y1, index.ii, index.ij, index.ei, index.ej, nx);
        forward_gradient_patch(vi_div2, grad_x2, grad_y2, index.ii, index.ij, index.ei, index.ej, nx);

#ifdef _OPENMP
#pragma omp parallel for
        for (int l = index.ij; l < index.ej; l++){
            for (int j = index.ii; j < index.ei; j++){
                const int i = l*nx + j;

                //Dual variables x11 and x12
                const float vec11 = g[i]*grad_x1[i];
                const float vec12 = g[i]*grad_y1[i];
                const float norm_vec1 = sqrt(vec11 * vec11 + vec12 * vec12);
                xi11[i] = (xi11[i] + tau_theta*vec11)/(1 + tau_theta*norm_vec1);
                xi12[i] = (xi12[i] + tau_theta*vec12)/(1 + tau_theta*norm_vec1);


                //Dual variables x21 and x22
                const float vec21 = g[i]*grad_x2[i];
                const float vec22 = g[i]*grad_y2[i];
                const float norm_vec2 = sqrt(vec21 * vec21 + vec22 * vec22);
                xi21[i] = (xi21[i] + tau_theta*vec21)/(1 + tau_theta*norm_vec2);
                xi22[i] = (xi22[i] + tau_theta*vec22)/(1 + tau_theta*norm_vec2);
            }
        }
#endif
    }
    //Compute divergence for last time
    for (int l = index.ij; l < index.ej; l++){
        for (int j = index.ii; j < index.ei; j++){
            const int i = l*nx + j;
            g_xi11[i] = g[i]*xi11[i];
            g_xi12[i] = g[i]*xi12[i];

            g_xi21[i] = g[i]*xi21[i];
            g_xi22[i] = g[i]*xi22[i];
        }
    }

    divergence_patch(g_xi11, g_xi12, div_g_xi1, index.ii, index.ij, index.ei, index.ej, nx);
    divergence_patch(g_xi21, g_xi22, div_g_xi2, index.ii, index.ij, index.ei, index.ej, nx);
}


//Minimization of occlusion variable
static void tvl2coupled_get_chi_patch(
        float *chi,
        float *chix,
        float *chiy,
        const float *F,
        const float *G,
        const float *g,
        float *eta1,
        float *eta2,
        const float *div_u,
        float *g_eta1,
        float *g_eta2,
        float *div_g_eta,
        const PatchIndexes index,
        Parameters params,
        const int nx
){

    // Added changes for subimages
    //int nx = params.w;
    for (int k = 1; k < ITER_CHI; k++){

        //Compute dual variable eta
        //#pragma omp simd collapse(2)
#ifdef _OPENMP
#pragma omp parallel for
        for (int l = index.ij; l < index.ej; l++){
            for (int j = index.ii; j < index.ei; j++){
                const int i = l*nx + j;

                //Compute new values of eta
                const float eta_new1 = eta1[i] + params.mu*params.tau_eta * g[i] * chix[i];
                const float eta_new2 = eta2[i] + params.mu*params.tau_eta * g[i] * chiy[i];
                const float norm_eta = sqrt(eta_new1*eta_new1 + eta_new2*eta_new2);
                //Put eta in the unit ball, [0, 1]
                if (norm_eta <= 1){
                    eta1[i] = eta_new1;
                    eta2[i] = eta_new2;
                }else{
                    eta1[i] = eta_new1/norm_eta;
                    eta2[i] = eta_new2/norm_eta;
                }

                //Compute values for the needed divergence
                g_eta1[i] = g[i]*eta1[i];
                g_eta2[i] = g[i]*eta2[i];
            }
        }
#endif

        divergence_patch(g_eta1, g_eta2, div_g_eta, index.ii, index.ij, index.ei, index.ej, nx);

        //Compute chi
        //#pragma omp simd collapse(2)
        //#pragma omp for schedule(dynamic, 1) collapse(2)
        for (int l = index.ij; l < index.ej; l++){
            for (int j = index.ii; j < index.ei; j++){
                const int i = l*nx + j;
                const float chi_new = chi[i] + params.tau_chi*(params.mu*div_g_eta[i] - params.beta*div_u[i] - F[i] - G[i]);
                chi[i] = max(min(chi_new, 1), 0);
            }
        }
        forward_gradient_patch(chi, chix, chiy, index.ii, index.ij, index.ei, index.ej, nx);
    }
    //Make thresholding in chi
    //#pragma omp for schedule(dynamic,1) collapse(2)
    for (int l = index.ij; l < index.ej; l++){
        for (int j = index.ii; j < index.ei; j++){
            const int i = l*nx + j;
            chi[i] = (chi[i] > THRESHOLD_DELTA) ? 1 : 0;
            //std::cout << chi[i] << '\n';
        }
    }
}



// Variational Optical flow method based on initial fixed values
// It minimizes the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)|
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)
void guided_tvl2coupled_occ(
        const float *I0,           // source image
        const float *I1,           // forward image
        const float *I_1,           // backward image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
        float *ener_N,
        const PatchIndexes index,
        const int nx,
        const int ny
        ) {


    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float *u1_ba = ofD->u1_ba;
    float *u2_ba = ofD->u2_ba;
    // w, h as params in the function call

    //Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;
    const int size = nx*ny;

    float *diff_u_N = tvl2_occ->diff_u_N;

    //Occlusion variables
    float *chi = ofD->chi;
    float *chix = tvl2_occ->chix;
    float *chiy = tvl2_occ->chiy;

    //Weigth
    float *g = tvl2_occ->g;

    //Dual variables
    float *xi11 = tvl2_occ->xi11;
    float *xi12 = tvl2_occ->xi12;
    float *xi21 = tvl2_occ->xi21;
    float *xi22 = tvl2_occ->xi22;

    //Dual variables of chi
    float *eta1 = tvl2_occ->eta1;
    float *eta2 = tvl2_occ->eta2;

    float *v1 = tvl2_occ->v1;
    float *v2 = tvl2_occ->v2;

    float *rho_c1 = tvl2_occ->rho_c1;
    float *rho_c_1 = tvl2_occ->rho_c_1;
    float *grad_1 = tvl2_occ->grad_1;
    float *grad__1 = tvl2_occ->grad__1;

    float *I0x = tvl2_occ->I0x;
    float *I0y = tvl2_occ->I0y;

    //Derivatives and warping of I2
    float *I1x = tvl2_occ->I1x;
    float *I1y = tvl2_occ->I1y;
    float *I1w = tvl2_occ->I1w;
    float *I1wx = tvl2_occ->I1wx;
    float *I1wy = tvl2_occ->I1wy;

    //Derivatives and warping of I0
    float *I_1x = tvl2_occ->I_1x;
    float *I_1y = tvl2_occ->I_1y;
    float *I_1w = tvl2_occ->I_1w;
    float *I_1wx = tvl2_occ->I_1wx;
    float *I_1wy = tvl2_occ->I_1wy;


    float *vi_div1 = tvl2_occ->vi_div1;
    float *grad_x1 = tvl2_occ->grad_x1;
    float *grad_y1 = tvl2_occ->grad_y1;
    float *vi_div2 = tvl2_occ->vi_div2;
    float *grad_x2 = tvl2_occ->grad_x2;
    float *grad_y2 = tvl2_occ->grad_y2;
    float *g_xi11 = tvl2_occ->g_xi11;
    float *g_xi12 = tvl2_occ->g_xi12;
    float *g_xi21 = tvl2_occ->g_xi21;
    float *g_xi22 = tvl2_occ->g_xi22;
    float *div_g_xi1 = tvl2_occ->div_g_xi1;
    float *div_g_xi2 = tvl2_occ->div_g_xi2;

    float *F = tvl2_occ->F;
    float *G = tvl2_occ->G;

    float *div_u = tvl2_occ->div_u;
    float *g_eta1 = tvl2_occ->g_eta1;
    float *g_eta2 = tvl2_occ->g_eta2;
    float *div_g_eta = tvl2_occ->div_g_eta;



    const float alpha = ofD->params.alpha;
    const float theta = ofD->params.theta;
    const float lambda = ofD->params.lambda;

    const float l_t = lambda * theta;

    //Initialization of dual variables and updating backward flow
    for (int l = index.ij; l < index.ej; l++){
        for (int k = index.ii; k < index.ei; k++){
            const int  i = l*nx + k;
            xi11[i] = xi12[i] = xi21[i] = xi22[i] = 0.0;

            u1_ba[i] = -u1[i];
            u2_ba[i] = -u2[i];

        }
    }

    //Initialize gradients from images and weight for global step
    if (ofD->params.step_algorithm == GLOBAL_STEP){

        //Compute derivatives of all images
        centered_gradient(I1, I1x, I1y, nx, ny);
        centered_gradient(I_1, I_1x, I_1y, nx, ny);
        centered_gradient(I0, I0x, I0y, nx, ny);
        //Initialize weight
        init_weight(g, I0x, I0y, size);
    }

    for (int warpings = 0; warpings < ofD->params.warps; warpings++) {
        // Compute the warping of I1 and its derivatives I1(x + u1o), I1x (x + u1o) and I1y (x + u2o)
        bicubic_interpolation_warp_patch(I1,  u1, u2, I1w, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy, index.ii, index.ij, index.ei, index.ej, nx, ny, false);

        // Compute the warping of I0 and its derivatives I0(x - u1o), I0x (x - u1o) and I0y (x - u2o)
        bicubic_interpolation_warp_patch(I_1,  u1_ba, u2_ba, I_1w, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I_1x, u1_ba, u2_ba, I_1wx, index.ii, index.ij, index.ei, index.ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I_1y, u1_ba, u2_ba, I_1wy, index.ii, index.ij, index.ei, index.ej, nx, ny, false);


        //Compute values that will not change during the whole wraping
        //#pragma omp simd collapse(2)
        for (int l = index.ij; l < index.ej; l++){
            for (int k = index.ii; k < index.ei; k++){

                const int i = l*nx + k;

                const float I1_x2 = I1wx[i] * I1wx[i];
                const float I1_y2 = I1wy[i] * I1wy[i];
                const float I_1_x2 = I_1wx[i] * I_1wx[i];
                const float I_1_y2 = I_1wy[i] * I_1wy[i];

                // store the |Grad(I2)|^2
                grad_1[i] = (I1_x2 + I1_y2);
                grad__1[i] = (I_1_x2 + I_1_y2);

                // Compute the constant part of the rho function
                rho_c1[i] = I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i];
                rho_c_1[i] = I_1w[i] - I_1wx[i] * u1[i]
                        - I_1wy[i] * u2[i] - I0[i];
            }
        }

        //Minimization of the functional
        int n = 0;
        float err_D = INFINITY;
        while (err_D > ofD->params.tol_OF*ofD->params.tol_OF && n < ofD->params.iterations_of){

            n++;
            // estimate the values of the variable (v1, v2)
#ifdef _OPENMP
#pragma omp parallel for
            for (int l = index.ij; l < index.ej; l++){
                for (int k = index.ii; k < index.ei; k++){
                    const int i = l*nx + k;
                    // rho function forward and backward
                    const float rho_1 = rho_c1[i]
                            + I1wx[i] * u1[i] + I1wy[i] * u2[i];
                    const float rho__1 = rho_c_1[i]
                            + I_1wx[i] * u1[i] + I_1wy[i] * u2[i];

                    //Stuff depending if pixel is occluded or not
                    int eps;
                    float alpha_i, mu, Lambda, grad, Iwx, Iwy, rho;
                    if (chi[i] == 0){
                        eps = 1;
                        alpha_i = 1;
                        mu = l_t;
                        Lambda = rho_1;
                        grad = grad_1[i];
                        Iwx = I1wx[i];
                        Iwy = I1wy[i];
                        rho = rho_1;
                    }else{
                        eps = -1;
                        alpha_i = 1/(1 + alpha*theta);
                        mu = l_t/(1 + alpha*theta);
                        Lambda = rho__1 +
                                alpha*theta/(1 + alpha*theta) * (u1[i]*I_1wx[i] + u2[i]*I_1wy[i]);
                        grad = grad__1[i];
                        Iwx = I_1wx[i];
                        Iwy = I_1wy[i];
                        rho = rho__1;
                    }
                    //Decide what to assign to v
                    if (Lambda > mu * grad){
                        v1[i] = alpha_i * u1[i] - mu * eps * Iwx;
                        v2[i] = alpha_i * u2[i] - mu * eps * Iwy;
                    }else{
                        if (Lambda < - mu * grad){
                            v1[i] = alpha_i * u1[i] + mu * eps * Iwx;
                            v2[i] = alpha_i * u2[i] + mu * eps * Iwy;
                        }else{
                            // if gradient is too small, we treat it as zero
                            if (grad < GRAD_IS_ZERO){
                                v1[i] = u1[i];
                                v2[i] = u2[i];
                            }else{

                                v1[i] = u1[i] - eps * rho * Iwx/grad;
                                v2[i] = u2[i] - eps * rho * Iwy/grad;
                            }
                        }
                    }
                }
            }
#endif
            //Estimate the values of the variable (u1, u2)
            //Compute derivatives of chi
            forward_gradient_patch(chi, chix, chiy, index.ii, index.ij, index.ei, index.ej, nx);

            //Compute dual variables
            tvl2coupled_get_xi_patch(xi11, xi12, xi21, xi22, g, v1, v2,
                                     chix, chiy, vi_div1, grad_x1, grad_y1, vi_div2, grad_x2, grad_y2,
                                     g_xi11, g_xi12, g_xi21, g_xi22, div_g_xi1, div_g_xi2,
                                     index, ofD->params, nx);

            //Compute several stuff
            //#pragma omp simd collapse(2)
            for (int l = index.ij; l < index.ej; l++){
                for (int k = index.ii; k < index.ei; k++){
                    const int i = l*nx + k;

                    //Previous value for u
                    const float u1k = u1[i];
                    const float u2k = u2[i];

                    //New value for (u1, u2)
                    u1[i] = v1[i] + theta*div_g_xi1[i] + theta * ofD->params.beta * chix[i];
                    u2[i] = v2[i] + theta*div_g_xi2[i] + theta * ofD->params.beta * chiy[i];

                    //Difference between previous and new value of u
                    diff_u_N[i] = (u1[i] - u1k) * (u1[i] - u1k) +
                            (u2[i] - u2k) * (u2[i] - u2k);


                    const float rho__1 = rho_c_1[i]
                            + I_1wx[i] * v1[i] + I_1wy[i] * v2[i];
                    const float rho_1 = rho_c1[i]
                            + I1wx[i] * v1[i] + I1wy[i] * v2[i];

                    F[i] = lambda*(std::abs(rho__1) - std::abs(rho_1));
                    G[i] = alpha/2*(v1[i]*v1[i] + v2[i]*v2[i]);
                }
            }

            //Compute chi
            tvl2coupled_get_chi_patch(chi, chix, chiy, F, G,
                                      g, eta1, eta2, div_u, g_eta1, g_eta2,
                                      div_g_eta, index, ofD->params, nx);


            //Get the max val
            float err_u = 0;
            int i = 0;
            for (int l = index.ij; l < index.ej; l++){
                for (int k = index.ii; k < index.ei; k++){
                    i = l*nx + k;
                    if (err_u < diff_u_N[i]){
                        err_u = diff_u_N[i];
                    }
                }
            }
            err_D = err_u;
        }
        if (ofD->params.step_algorithm == GLOBAL_STEP && ofD->params.verbose)
            std::printf("Warping: %d, Iter: %d "
                        "Error: %f\n", warpings, n, err_D);
    }
    if (ofD->params.step_algorithm == LOCAL_STEP){
        *ener_N = eval_tvl2coupled_occ(I0, I1, I_1, ofD, tvl2_occ, index, ofD->params, nx, ny);
    }
}



#endif //TVL2-L1 functional
