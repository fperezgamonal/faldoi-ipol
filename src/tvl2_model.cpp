#ifndef TVL2_MODEL
#define TVL2_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "utils.h"

extern "C" {
#include "bicubic_interpolation.h"
}

#include <omp.h>

//// INITIALIZATION OF EACH METHOD
void intialize_stuff_tvl2coupled(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore, const int w, const int h) {

    // Added w, h in as function params
    //const int w = ofCore->params.w;
    //const int h = ofCore->params.h;
    //fprintf(stderr, "W x H :%d x %d\n", w, h);
    // Dual variables
    ofStuff->tvl2.xi11 = new float[w * h];
    ofStuff->tvl2.xi12 = new float[w * h];
    ofStuff->tvl2.xi21 = new float[w * h];
    ofStuff->tvl2.xi22 = new float[w * h];

    ofStuff->tvl2.u1x = new float[w * h];
    ofStuff->tvl2.u1y = new float[w * h];
    ofStuff->tvl2.u2x = new float[w * h];
    ofStuff->tvl2.u2y = new float[w * h];

    ofStuff->tvl2.v1 = new float[w * h];
    ofStuff->tvl2.v2 = new float[w * h];

    ofStuff->tvl2.rho_c = new float[w * h];
    ofStuff->tvl2.grad = new float[w * h];
    ofStuff->tvl2.u1_ = new float[w * h];
    ofStuff->tvl2.u2_ = new float[w * h];
    ofStuff->tvl2.u1Aux = new float[w * h];
    ofStuff->tvl2.u2Aux = new float[w * h];
    ofStuff->tvl2.I1x = new float[w * h];
    ofStuff->tvl2.I1y = new float[w * h];
    ofStuff->tvl2.I1w = new float[w * h];
    ofStuff->tvl2.I1wx = new float[w * h];
    ofStuff->tvl2.I1wy = new float[w * h];
    ofStuff->tvl2.div_xi1 = new float[w * h];
    ofStuff->tvl2.div_xi2 = new float[w * h];
    ofStuff->tvl2.u_N = new float[w * h];
}

void free_stuff_tvl2coupled(SpecificOFStuff *ofStuff) {
    delete[] ofStuff->tvl2.xi11;
    delete[] ofStuff->tvl2.xi12;
    delete[] ofStuff->tvl2.xi21;
    delete[] ofStuff->tvl2.xi22;
    delete[] ofStuff->tvl2.u1x;
    delete[] ofStuff->tvl2.u1y;
    delete[] ofStuff->tvl2.u2x;
    delete[] ofStuff->tvl2.u2y;
    delete[] ofStuff->tvl2.v1;
    delete[] ofStuff->tvl2.v2;
    delete[] ofStuff->tvl2.rho_c;
    delete[] ofStuff->tvl2.grad;
    delete[] ofStuff->tvl2.u1_;
    delete[] ofStuff->tvl2.u2_;
    delete[] ofStuff->tvl2.u1Aux;
    delete[] ofStuff->tvl2.u2Aux;
    delete[] ofStuff->tvl2.I1x;
    delete[] ofStuff->tvl2.I1y;
    delete[] ofStuff->tvl2.I1w;
    delete[] ofStuff->tvl2.I1wx;
    delete[] ofStuff->tvl2.I1wy;
    delete[] ofStuff->tvl2.div_xi1;
    delete[] ofStuff->tvl2.div_xi2;
    delete[] ofStuff->tvl2.u_N;
}

//////////////////////////////////////////////
////TV-l2 COUPLED OPTICAL FLOW PROBLEM////////
// Dual variable
static void tvl2coupled_getD(
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float *u1x,
        const float *u1y,
        const float *u2x,
        const float *u2y,
        float tau,
        const int ii,           // initial column
        const int ij,           // initial row
        const int ei,           // end column
        const int ej,           // end row
        const int nx
) {
    // Compute the value of xi
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int l = ij; l < ej; l++) {
        for (int k = ii; k < ei; k++) {
            const int i = l * nx + k;
            const float g11 = xi11[i] * xi11[i];
            const float g12 = xi12[i] * xi12[i];
            const float g21 = xi21[i] * xi21[i];
            const float g22 = xi22[i] * xi22[i];

            float xi_N = sqrt(g11 + g12 + g21 + g22);

            xi_N = MAX(1, xi_N);

            xi11[i] = (xi11[i] + tau * u1x[i]) / xi_N;
            xi12[i] = (xi12[i] + tau * u1y[i]) / xi_N;
            xi21[i] = (xi21[i] + tau * u2x[i]) / xi_N;
            xi22[i] = (xi22[i] + tau * u2y[i]) / xi_N;
        }
    }
}


//Primal variable
static void tvl2coupled_getP(
        float *u1,
        float *u2,
        const float *v1,
        const float *v2,
        const float *div_xi1,
        const float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        const int ii,           // initial column
        const int ij,           // initial row
        const int ei,           // end column
        const int ej,           // end row
        const int nx,
        float *err
) {
    float err_D = 0.0;

//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++) {
            const int i = l * nx + k;
            const float u1k = u1[i];
            const float u2k = u2[i];

            // Only modify the inpainting domain
            // if (mask[i]==0){
            u1[i] = u1k - tau * (-div_xi1[i] + (u1k - v1[i]) / theta);
            u2[i] = u2k - tau * (-div_xi2[i] + (u2k - v2[i]) / theta);

            u_N[i] = (u1[i] - u1k) * (u1[i] - u1k) +
                     (u2[i] - u2k) * (u2[i] - u2k);
            // }else {
            //   u_N[i] = 0.0;
            // }
        }


    //Get the max val
    err_D = 0;
    for (int l = ij; l < ej; l++) {
        for (int k = ii; k < ei; k++) {
            if (err_D < u_N[l * nx + k]) {
                err_D = u_N[l * nx + k];
            }
        }
    }
    (*err) = err_D;

}

void eval_tvl2coupled(
        const float *I0,            // source image
        const float *I1,            // target image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff *tvl2,
        float *ener_N,
        const int ii,               // initial column
        const int ij,               // initial row
        const int ei,               // end column
        const int ej,               // end row
        const float lambda,         // weight of the data term
        const float theta,
        const int nx,
        const int ny
) {

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;


    // Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;

    // Optical flow derivatives
    float *v1 = tvl2->v1;
    float *v2 = tvl2->v2;
    float *u1x = tvl2->u1x;
    float *u1y = tvl2->u1y;
    float *u2x = tvl2->u2x;
    float *u2y = tvl2->u2y;

    float *I1w = tvl2->I1w;

    float ener = 0.0;

    //forward_gradient_mixed_bound_patch(u1,u1x,u1y,ii,ij,ei,ej,nx,ny);
    //forward_gradient_mixed_bound_patch(u2,u2x,u2y,ii,ij,ei,ej,nx,ny);
    forward_gradient_patch(u1, u1x, u1y, ii, ij, ei, ej, nx);
    forward_gradient_patch(u2, u2x, u2y, ii, ij, ei, ej, nx);
    bicubic_interpolation_warp_patch(I1, u1, u2, I1w,
                                     ii, ij, ei, ej, nx, ny, false);

    // Energy for all the patch. Maybe it is useful only the 8 pixel around the seed.
    int m = 0;
    for (int l = ij; l < ej; l++) {
        for (int k = ii; k < ei; k++) {
            const int i = l * nx + k;
            float dt = lambda * fabs(I1w[i] - I0[i]);
            float dc = (1 / (2 * theta)) *
                       ((u1[i] - v1[i]) * (u1[i] - v1[i]) + (u2[i] - v2[i]) * (u2[i] - v2[i]));
            float g1 = u1x[i] * u1x[i];
            float g12 = u1y[i] * u1y[i];
            float g21 = u2x[i] * u2x[i];
            float g2 = u2y[i] * u2y[i];
            float g = sqrt(g1 + g12 + g21 + g2);
            if (!std::isfinite(dt)) {
                std::printf("Corrupt data\n");
            }
            if (!std::isfinite(g)) {
                std::printf("Corrupt regularization\n");
            }
            ener += dc + dt + g;
            m++;
        }
    }
    ener /= (m * 1.0);
    (*ener_N) = ener;
    assert(ener >= 0.0);
}

// Variational Optical flow method based on initial fixed values
// It minimizes the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)|
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)
void guided_tvl2coupled(
        const float *I0,            // source image
        const float *I1,            // target image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff *tvl2,
        float *ener_N,
        const int ii,               // initial column
        const int ij,               // initial row
        const int ei,               // end column
        const int ej,               // end row
        const float lambda,         // weight of the data term
        const float theta,          // weight of the data term
        const float tau,            // time step
        const float tol_OF,         // tol max allowed
        const int warps,            // number of warpings per scale
        const bool verbose,         // enable/disable the verbose mode
        const int nx,               // width of I0 (I1 may have diff. size if we use partitions)
        const int ny                // height of I0 ( " " " " " " " " ")
) {

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;

    // Columns and Rows
    //const int nx = ofD->params.w;
    //const int ny = ofD->params.h;


    float *u1_ = tvl2->u1_;
    float *u2_ = tvl2->u2_;
    float *u_N = tvl2->u_N;

    // Optical flow derivatives
    float *u1x = tvl2->u1x;
    float *u1y = tvl2->u1y;
    float *u2x = tvl2->u2x;
    float *u2y = tvl2->u2y;

    // Dual variables
    float *xi11 = tvl2->xi11;
    float *xi12 = tvl2->xi12;
    float *xi21 = tvl2->xi21;
    float *xi22 = tvl2->xi22;

    float *v1 = tvl2->v1;
    float *v2 = tvl2->v2;

    float *rho_c = tvl2->rho_c;
    float *grad = tvl2->grad;

    float *u1Aux = tvl2->u1Aux;
    float *u2Aux = tvl2->u2Aux;

    float *I1x = tvl2->I1x;
    float *I1y = tvl2->I1y;

    float *I1w = tvl2->I1w;
    float *I1wx = tvl2->I1wx;
    float *I1wy = tvl2->I1wy;

    // Divergence
    float *div_xi1 = tvl2->div_xi1;
    float *div_xi2 = tvl2->div_xi2;

    const float l_t = lambda * theta;

    for (int l = ij; l < ej; l++) {
        for (int k = ii; k < ei; k++) {
            const int i = l * nx + k;
            //Initialization dual variables
            xi11[i] = xi12[i] = xi21[i] = xi22[i] = 0.0;
        }
    }
    //int tmp_count = 0;
    omp_set_dynamic(0); // Consider changing for optimal cluster performance (maybe)
//omp_set_num_threads(2);
    for (int warpings = 0; warpings < warps; warpings++) {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp_patch(I1, u1, u2, I1w,
                                         ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx,
                                         ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy,
                                         ii, ij, ei, ej, nx, ny, false);

        //Compute values that will not change during the whole wraping
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
        for (int l = ij; l < ej; l++)
            for (int k = ii; k < ei; k++) {
                /*
                ++tmp_count;
                printf("Combination of warping, l and k nº= %d\n", tmp_count);
                */
                const int i = l * nx + k;
                const float Ix2 = I1wx[i] * I1wx[i];
                const float Iy2 = I1wy[i] * I1wy[i];

                // store the |Grad(I1)|^2
                grad[i] = (Ix2 + Iy2);

                // compute the constant part of the rho function
                rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                            - I1wy[i] * u2[i] - I0[i]);
            }

//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
        for (int l = ij; l < ej; l++) {
            for (int k = ii; k < ei; k++) {
                const int i = l * nx + k;
                u1_[i] = u1[i];
                u2_[i] = u2[i];
            }
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF * tol_OF && n < ofD->params.max_iter_patch) {

            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++) {
                for (int k = ii; k < ei; k++) {
                    const int i = l * nx + k;
                    const float rho = rho_c[i]
                                      + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                    float d1, d2;

                    if (rho < -l_t * grad[i]) {
                        d1 = l_t * I1wx[i];
                        d2 = l_t * I1wy[i];
                    } else {
                        if (rho > l_t * grad[i]) {
                            d1 = -l_t * I1wx[i];
                            d2 = -l_t * I1wy[i];
                        } else {
                            // if gradient is too small, we treat it as zero
                            if (grad[i] < GRAD_IS_ZERO)
                                d1 = d2 = 0;
                            else {
                                float fi = -rho / grad[i];
                                d1 = fi * I1wx[i];
                                d2 = fi * I1wy[i];
                            }
                        }
                    }
                    v1[i] = u1[i] + d1;
                    v2[i] = u2[i] + d2;
                }
            }
            // Estimate the values of the variable (u1, u2)

            // Compute dual variables
            forward_gradient_patch(u1_, u1x, u1y, ii, ij, ei, ej, nx);
            forward_gradient_patch(u2_, u2x, u2y, ii, ij, ei, ej, nx);
            tvl2coupled_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y,
                             tau, ii, ij, ei, ej, nx);

            // Primal variables
            divergence_patch(xi11, xi12, div_xi1, ii, ij, ei, ej, nx);
            divergence_patch(xi21, xi22, div_xi2, ii, ij, ei, ej, nx);

            // Save previous iteration
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++) {
                for (int k = ii; k < ei; k++) {
                    const int i = l * nx + k;
                    u1Aux[i] = u1[i];
                    u2Aux[i] = u2[i];
                }
            }

            tvl2coupled_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N,
                             theta, tau, ii, ij, ei, ej, nx, &err_D);

            // (acceleration = 1);
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++) {
                for (int k = ii; k < ei; k++) {
                    const int i = l * nx + k;
                    u1_[i] = 2 * u1[i] - u1Aux[i];
                    u2_[i] = 2 * u2[i] - u2Aux[i];
                }
            }


        }
        if (ofD->params.step_algorithm == GLOBAL_STEP && ofD->params.verbose)
            std::printf("Warping: %d,Iter: %d "
                                "Error: %f\n", warpings, n, err_D);
    }
    eval_tvl2coupled(I0, I1, ofD, tvl2, ener_N, ii, ij, ei, ej, lambda, theta, nx, ny);
}

#endif //TVL2-L1 functional