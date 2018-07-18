#ifndef TVL2W_MODEL
#define TVL2W_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "utils.h"
extern "C" {
#include "bicubic_interpolation.h"
}

////INITIALIZATION OF EACH METHOD
void  initialize_stuff_tvl2coupled_w(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        const int w,
        const int h)

{
    //fprintf(stderr, "W x H :%d x %d\n", w, h);
    ofStuff->tvl2w.weight = new float[ofCore->params.w_radio*2 + 1];
    ofStuff->tvl2w.xi11 = new float[w*h];
    ofStuff->tvl2w.xi12 = new float[w*h];
    ofStuff->tvl2w.xi21 = new float[w*h];
    ofStuff->tvl2w.xi22 = new float[w*h];
    ofStuff->tvl2w.u1x = new float[w*h];
    ofStuff->tvl2w.u1y = new float[w*h];
    ofStuff->tvl2w.u2x = new float[w*h];
    ofStuff->tvl2w.u2y = new float[w*h];
    ofStuff->tvl2w.v1 =  new float[w*h];
    ofStuff->tvl2w.v2 =  new float[w*h];
    ofStuff->tvl2w.rho_c =  new float[w*h];
    ofStuff->tvl2w.grad =  new float[w*h];
    ofStuff->tvl2w.u1_ =  new float[w*h];
    ofStuff->tvl2w.u2_ =  new float[w*h];
    ofStuff->tvl2w.u1Aux = new float[w*h];
    ofStuff->tvl2w.u2Aux = new float[w*h];
    ofStuff->tvl2w.I1x = new float[w*h];
    ofStuff->tvl2w.I1y = new float[w*h];
    ofStuff->tvl2w.I1w = new float[w*h];
    ofStuff->tvl2w.I1wx = new float[w*h];
    ofStuff->tvl2w.I1wy = new float[w*h];
    ofStuff->tvl2w.div_xi1 = new float[w*h];
    ofStuff->tvl2w.div_xi2 = new float[w*h];
    ofStuff->tvl2w.u_N = new float[w*h];
}

void  free_stuff_tvl2coupled_w(SpecificOFStuff *ofStuff){
    delete [] ofStuff->tvl2w.weight;
    delete [] ofStuff->tvl2w.xi11;
    delete [] ofStuff->tvl2w.xi12;
    delete [] ofStuff->tvl2w.xi21;
    delete [] ofStuff->tvl2w.xi22;
    delete [] ofStuff->tvl2w.u1x;
    delete [] ofStuff->tvl2w.u1y;
    delete [] ofStuff->tvl2w.u2x;
    delete [] ofStuff->tvl2w.u2y;
    delete [] ofStuff->tvl2w.v1;
    delete [] ofStuff->tvl2w.v2;
    delete [] ofStuff->tvl2w.rho_c;
    delete [] ofStuff->tvl2w.grad;
    delete [] ofStuff->tvl2w.u1_;
    delete [] ofStuff->tvl2w.u2_;
    delete [] ofStuff->tvl2w.u1Aux;
    delete [] ofStuff->tvl2w.u2Aux;
    delete [] ofStuff->tvl2w.I1x;
    delete [] ofStuff->tvl2w.I1y;
    delete [] ofStuff->tvl2w.I1w;
    delete [] ofStuff->tvl2w.I1wx;
    delete [] ofStuff->tvl2w.I1wy;
    delete [] ofStuff->tvl2w.div_xi1;
    delete [] ofStuff->tvl2w.div_xi2;
    delete [] ofStuff->tvl2w.u_N;
}


//////////////////////////////////////////
////TV-l2 COUPLED OPTICAL FLOW PROBLEM////////
//Dual variable
inline void tvl2coupled_w_getD(
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        float *u1x,
        float *u1y,
        float *u2x,
        float *u2y,
        float tau,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int nx
        ){
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            const float g11 = xi11[i]*xi11[i];
            const float g12 = xi12[i]*xi12[i];
            const float g21 = xi21[i]*xi21[i];
            const float g22 = xi22[i]*xi22[i];

            float xi_N = sqrt(g11 + g12 + g21 + g22);

            xi_N = MAX(1,xi_N);

            xi11[i] = (xi11[i] + tau*u1x[i])/xi_N;
            xi12[i] = (xi12[i] + tau*u1y[i])/xi_N;
            xi21[i] = (xi21[i] + tau*u2x[i])/xi_N;
            xi22[i] = (xi22[i] + tau*u2y[i])/xi_N;
        }
    }
}


//Primal variable
inline void tvl2coupled_w_getP(
        float *u1,
        float *u2,
        float *v1,
        float *v2,
        float *div_xi1,
        float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int nx,
        float *err
        ){
    float err_D = 0.0;
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            const float u1k = u1[i];
            const float u2k = u2[i];

            //Only modify the inpainting domain
            // if (mask[i]==0){
            u1[i] = u1k -tau*(-div_xi1[i] + (u1k - v1[i])/theta);
            u2[i] = u2k -tau*(-div_xi2[i] + (u2k - v2[i])/theta);

            u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
                    (u2[i] - u2k) * (u2[i] - u2k);
            // }else {
            //   u_N[i] = 0.0;
            // }
        }


    //Get the max val
    err_D = 0;
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            if (err_D < u_N[l*nx + k]){
                err_D = u_N[l*nx + k];
            }
        }
    }

    (*err) = err_D;

}

void eval_tvl2coupled_w(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_W *tvl2w,
        float *ener_N,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const float lambda,  // weight of the data term
        const float theta,
        const int nx,
        const int ny
        )
{

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;

    // Added changes for subimages

    //Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;

    //Optical flow derivatives
    float *v1   = tvl2w->v1;
    float *v2   = tvl2w->v2;
    float *u1x  = tvl2w->u1x;
    float *u1y  = tvl2w->u1y;
    float *u2x  = tvl2w->u2x;
    float *u2y  = tvl2w->u2y;

    float *I1w = tvl2w->I1w;

    float ener = 0.0;

    //TODO:Pesos
    const int iiw = tvl2w->iiw;
    const int ijw = tvl2w->ijw;
    float *weight = tvl2w->weight;

    //forward_gradient_mixed_bound_patch(u1,u1x,u1y,ii,ij,ei,ej,nx,ny);
    //forward_gradient_mixed_bound_patch(u2,u2x,u2y,ii,ij,ei,ej,nx,ny);
    forward_gradient_patch(u1,u1x,u1y,ii,ij,ei,ej,nx);
    forward_gradient_patch(u2,u2x,u2y,ii,ij,ei,ej,nx);
    bicubic_interpolation_warp_patch(I1,  u1, u2, I1w,
                                     ii, ij, ei, ej, nx, ny, false);
    //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
    int m  = 0;
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            float dt = lambda*fabs(I1w[i]-I0[i])*weight[l-ij + ijw]*weight[k-ii + iiw];
            float dc = (1/(2*theta))*
                    ((u1[i]-v1[i])*(u1[i]-v1[i]) + (u2[i] - v2[i])*(u2[i] - v2[i]));
            float g1  = u1x[i]*u1x[i];
            float g12 = u1y[i]*u1y[i];
            float g21 = u2x[i]*u2x[i];
            float g2  = u2y[i]*u2y[i];
            float g  = sqrt(g1 + g12 + g21 + g2);
            if (!std::isfinite(dt)){
                std::printf("Corrupt data\n");
            }
            if (!std::isfinite(g)){
                std::printf("Corrupt regularization\n");
            }
            ener += dc + dt + g;
            m++;
        }
    }
    ener /=(m*1.0);
    (*ener_N) = ener;
    assert(ener >= 0.0);
}

// Variational Optical flow method based on initial fixed values
// It minimize the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)| 
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)
void guided_tvl2coupled_w(
        const float *I0,                // source image
        const float *I1,                // target image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_W *tvl2w,
        float *ener_N,
        const int ii,                   // initial column
        const int ij,                   // initial row
        const int ei,                   // end column
        const int ej,                   // end row
        const float lambda,             // weight of the data term
        const float theta,              // weight of the data term
        const float tau,                // time step
        const float tol_OF,             // tol max allowed
        const int   warps,              // number of warpings per scale
        const bool  verbose,            // enable/disable the verbose mode
        const int nx,
        const int ny){

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    // w, h as params in the function call


    //Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;


    float *u1_  = tvl2w->u1_;
    float *u2_  = tvl2w->u2_;
    float *u_N  = tvl2w->u_N;;
    //Optical flow derivatives
    float *u1x  = tvl2w->u1x;
    float *u1y  = tvl2w->u1y;
    float *u2x  = tvl2w->u2x;
    float *u2y  = tvl2w->u2y;
    //Dual variables
    float *xi11 = tvl2w->xi11;
    float *xi12 = tvl2w->xi12;
    float *xi21 = tvl2w->xi21;
    float *xi22 = tvl2w->xi22;

    float *v1 = tvl2w->v1;
    float *v2 = tvl2w->v2;

    float *rho_c = tvl2w->rho_c;
    float *grad  = tvl2w->grad;

    float *u1Aux = tvl2w->u1Aux;
    float *u2Aux = tvl2w->u2Aux;

    float *I1x = tvl2w->I1x;
    float *I1y = tvl2w->I1y;

    float *I1w = tvl2w->I1w;
    float *I1wx = tvl2w->I1wx;
    float *I1wy = tvl2w->I1wy;

    //Divergence
    float *div_xi1 = tvl2w->div_xi1;
    float *div_xi2 = tvl2w->div_xi2;


    //TODO: Weights
    const int iiw = tvl2w->iiw;
    const int ijw = tvl2w->ijw;
    float *weight = tvl2w->weight;

    const float l_t = lambda * theta;

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int  i = l*nx + k;
            //Inizialization dual variables
            xi11[i] = xi12[i] = xi21[i] = xi22[i] = 0.0;
        }
    }

    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp_patch(I1,  u1, u2, I1w,
                                         ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx,
                                         ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy,
                                         ii, ij, ei, ej, nx, ny, false);

        for (int l = ij; l < ej; l++)
            for (int k = ii; k < ei; k++)
            {

                const int i = l*nx + k;
                const float Ix2 = I1wx[i] * I1wx[i];
                const float Iy2 = I1wy[i] * I1wy[i];

                // store the |Grad(I1)|^2
                grad[i] = (Ix2 + Iy2);

                // compute the constant part of the rho function
                rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                            - I1wy[i] * u2[i] - I0[i]);
            }


        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
                const int i = l*nx + k;
                u1_[i] = u1[i];
                u2_[i] = u2[i];
            }
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF*tol_OF && n < ofD->params.max_iter_patch)
        {

            n++;
            // estimate the values of the variable (v1, v2)
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const float l_t_w = l_t * weight[l-ij + ijw]*weight[k-ii + iiw];
                    const int i = l*nx + k;
                    const float rho = rho_c[i]
                            + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                    float d1, d2;

                    if (rho < - l_t_w * grad[i])
                    {
                        d1 = l_t_w * I1wx[i];
                        d2 = l_t_w * I1wy[i];
                    }
                    else
                    {
                        if (rho > l_t_w * grad[i])
                        {
                            d1 = -l_t_w * I1wx[i];
                            d2 = -l_t_w * I1wy[i];
                        }
                        else
                        {
                            if (grad[i] < GRAD_IS_ZERO)
                                d1 = d2 = 0;
                            else
                            {
                                float fi = -rho/grad[i];
                                d1 = fi * I1wx[i];
                                d2 = fi * I1wy[i];
                            }
                        }
                    }

                    v1[i] = u1[i] + d1;
                    v2[i] = u2[i] + d2;
                }
            }

            //Dual variables
            forward_gradient_patch(u1_,u1x,u1y,ii,ij,ei,ej,nx);
            forward_gradient_patch(u2_,u2x,u2y,ii,ij,ei,ej,nx);
            tvl2coupled_w_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y,
                               tau, ii, ij, ei, ej, nx);
            //Primal variables
            divergence_patch(xi11,xi12,div_xi1,ii,ij,ei,ej,nx);
            divergence_patch(xi21,xi22,div_xi2,ii,ij,ei,ej,nx);

            //Almacenamos la iteracion anterior
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*nx + k;
                    u1Aux[i] = u1[i];
                    u2Aux[i] = u2[i];
                }
            }

            tvl2coupled_w_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N,
                                theta, tau, ii, ij, ei, ej, nx, &err_D);
            //(aceleration = 1);

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*nx + k;
                    u1_[i] = 2*u1[i] - u1Aux[i];
                    u2_[i] = 2*u2[i] - u2Aux[i];
                }
            }


        }
        if (verbose)
            std::printf("Warping: %d,Iter: %d "
                        "Error: %f\n", warpings,n, err_D);
    }
    eval_tvl2coupled_w(I0, I1, ofD, tvl2w, ener_N, ii, ij, ei, ej, lambda, theta, nx, ny);
}

#endif //TVL2-L1 functional
