// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license athis program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2013, Roberto P.Palomares <roberto.palomares@upf.edu>
// Copyright (C) 2018, Ferran Pérez <fperez.gamonal@gmail.com>
// All rights reserved.

#ifndef GLOBAL_FALDOI
#define GLOBAL_FALDOI

#include <cmath>
#include <cstdio>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <vector>
#include <algorithm>

extern "C" {
#include "bicubic_interpolation.h"
#include "iio.h"
}

#include "tvl2_model_occ.h"
#include "utils.h"
#include "parameters.h"
#include "utils_preprocess.h"
#include "energy_model.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

using namespace std;


#define MAX(x, y) ((x)>(y)?(x):(y))


//////////
///////////WARNING
/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(
        const float *I0,        // source image
        float *I1,              // target image
        float *u1,              // x component of the optical flow
        float *u2,              // y component of the optical flow
        const int nx,           // image width
        const int ny,           // image height
        const float tau,        // time step
        const float lambda,     // weight parameter for the data term
        const float theta,      // weight parameter for (u - v)²
        const int warps,        // number of warpings per scale
        const float epsilon,    // tolerance for numerical convergence
        const bool verbose      // enable/disable the verbose mode
) {
    const int size = nx * ny;
    const float l_t = lambda * theta;

    auto *I1x = new float[size];
    auto *I1y = new float[size];
    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];
    auto *rho_c = new float[size];
    auto *v1 = new float[size];
    auto *v2 = new float[size];
    auto *p11 = new float[size];
    auto *p12 = new float[size];
    auto *p21 = new float[size];
    auto *p22 = new float[size];
    auto *div = new float[size];
    auto *grad = new float[size];
    auto *div_p1 = new float[size];
    auto *div_p2 = new float[size];
    auto *u1x = new float[size];
    auto *u1y = new float[size];
    auto *u2x = new float[size];
    auto *u2y = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    // Initialization of p
    for (int i = 0; i < size; i++) {
        p11[i] = p12[i] = 0.0;
        p21[i] = p22[i] = 0.0;
    }

    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of the target image and its derivatives
        bicubic_interpolation_warp(I1, u1, u2, I1w, nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

        //#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // Compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        int n = 0;
        float error = INFINITY;
        while (error > epsilon * epsilon && n < MAX_ITERATIONS_GLOBAL) {
            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++) {
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

            // Compute the divergence of the dual variable (p1, p2)
            divergence(p11, p12, div_p1, nx, ny);
            divergence(p21, p22, div_p2, nx, ny);

            // Estimate the values of the optical flow (u1, u2)
            error = 0.0;
            //#pragma omp parallel for reduction(+:error)
            for (int i = 0; i < size; i++) {
                const float u1k = u1[i];
                const float u2k = u2[i];

                u1[i] = v1[i] + theta * div_p1[i];
                u2[i] = v2[i] + theta * div_p2[i];

                error += (u1[i] - u1k) * (u1[i] - u1k) +
                         (u2[i] - u2k) * (u2[i] - u2k);
            }
            error /= size;

            // Compute the gradient of the optical flow (Du1, Du2)
            forward_gradient(u1, u1x, u1y, nx, ny);
            forward_gradient(u2, u2x, u2y, nx, ny);

            // Estimate the values of the dual variable (p1, p2)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                const float taut = tau / theta;
                const float g1 = hypotf(u1x[i], u1y[i]);
                const float g2 = hypotf(u2x[i], u2y[i]);
                const float ng1 = 1.0 + taut * g1;
                const float ng2 = 1.0 + taut * g2;

                p11[i] = (p11[i] + taut * u1x[i]) / ng1;
                p12[i] = (p12[i] + taut * u1y[i]) / ng1;
                p21[i] = (p21[i] + taut * u2x[i]) / ng2;
                p22[i] = (p22[i] + taut * u2y[i]) / ng2;
            }
        }

        if (verbose)
            fprintf(stderr, "Warping: %d, "
                    "Iterations: %d, "
                    "Error: %f\n", warpings, n, error);
    }

    // Delete allocated memory
    free(I1x);
    free(I1y);
    free(I1w);
    free(I1wx);
    free(I1wy);
    free(rho_c);
    free(v1);
    free(v2);
    free(p11);
    free(p12);
    free(p21);
    free(p22);
    free(div);
    free(grad);
    free(div_p1);
    free(div_p2);
    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);
}

// What is this warning for?
////////////WARNING////////
///////////////////////////////////
/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void ofDu_getP(
        float *u1,
        float *u2,
        const float *v1,
        const float *v2,
        const float *div_xi1,
        const float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        int size,
        float *err
) {
    float err_D = 0.0;
    float min, max;

    //#pragma omp parallel for reduction(+:err_D)
    for (int i = 0; i < size; i++) {

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k - tau * (-div_xi1[i] + (u1k - v1[i]) / theta);
        u2[i] = u2k - tau * (-div_xi2[i] + (u2k - v2[i]) / theta);

        u_N[i] = (u1[i] - u1k) * (u1[i] - u1k) +
                 (u2[i] - u2k) * (u2[i] - u2k);
    }

    getminmax(&min, &max, u_N, size);

    err_D = max;
    (*err) = err_D;
}


/*
 * - Name: getD_Du

 *
*/
void ofDu_getD(
        float *xi11,
        float *xi12,
        float *xi22,
        const float *u1x,
        const float *u1y,
        const float *u2x,
        const float *u2y,
        const float tau,
        int size
) {
    //#pragma omp parallel for

    for (int i = 0; i < size; i++) {

        const float g11 = xi11[i] * xi11[i];
        const float g12 = xi12[i] * xi12[i];
        const float g22 = xi22[i] * xi22[i];

        float xi_N = sqrt(g11 + g22 + 2 * g12);

        xi_N = MAX(1, xi_N);

        xi11[i] = (xi11[i] + tau * u1x[i]) / xi_N;
        xi12[i] = (xi12[i] + 0.5 * tau * (u1y[i] + u2x[i])) / xi_N;
        xi22[i] = (xi22[i] + tau * u2y[i]) / xi_N;
    }
}


/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void ofTVl2_getP(
        float *u1,
        float *u2,
        const float *v1,
        const float *v2,
        const float *div_xi1,
        const float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        int size,
        float *err
) {
    float err_D = 0.0;
    float min, max;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++) {

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k - tau * (-div_xi1[i] + (u1k - v1[i]) / theta);
        u2[i] = u2k - tau * (-div_xi2[i] + (u2k - v2[i]) / theta);

        u_N[i] = (u1[i] - u1k) * (u1[i] - u1k) +
                 (u2[i] - u2k) * (u2[i] - u2k);
    }
#endif

    getminmax(&min, &max, u_N, size);

    err_D = max;
    (*err) = err_D;
}

/*
 * - Name: ofTVl2_getD

 *
*/
void ofTVl2_getD(
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float *u1x,
        const float *u1y,
        const float *u2x,
        const float *u2y,
        float tau,
        int size
) {

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++) {

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
#endif
}


void duOF(
        const float *I0,        // source image
        float *I1,              // target image
        float *u1,              // x component of the optical flow
        float *u2,              // y component of the optical flow
        float *xi11,
        float *xi12,
        float *xi22,
        const float lambda,     // weight of the data term
        const float theta,      // weight of the data term
        const float tau,        // time step
        const float tol_OF,     // tol max allowed
        const int nx,           // image width
        const int ny,           // image height
        const int warps,        // number of warpings per scale
        const bool verbose      // enable/disable the verbose mode
) {

    const float l_t = lambda * theta;
    const int size = nx * ny;


    auto *u1x = new float[size];
    auto *u1y = new float[size];
    auto *u2x = new float[size];
    auto *u2y = new float[size];

    auto *v1 = new float[size];
    auto *v2 = new float[size];

    auto *rho_c = new float[size];
    auto *grad = new float[size];

    auto *u1_ = new float[size];
    auto *u2_ = new float[size];

    auto *u1Aux = new float[size];
    auto *u2Aux = new float[size];

    auto *I1x = new float[size];
    auto *I1y = new float[size];

    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];

    // Divergence
    auto *div_xi1 = new float[size];
    auto *div_xi2 = new float[size];

    auto *u_N = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1, u1, u2, I1w, nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // Compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        for (int i = 0; i < nx * ny; i++) {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF * tol_OF && n < MAX_ITERATIONS_GLOBAL) {

            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++) {
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

            // Dual variables
            forward_gradient(u1_, u1x, u1y, nx, ny);
            forward_gradient(u2_, u2x, u2y, nx, ny);
            ofDu_getD(xi11, xi12, xi22, u1x, u1y, u2x, u2y, tau, size);

            // Primal variables
            divergence(xi11, xi12, div_xi1, nx, ny);
            divergence(xi12, xi22, div_xi2, nx, ny);

            // Store previous iteration
            for (int i = 0; i < size; i++) {
                u1Aux[i] = u1[i];
                u2Aux[i] = u2[i];
            }

            ofDu_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N, theta, tau, size, &err_D);

            // (acceleration = 1);
            for (int i = 0; i < size; i++) {
                u1_[i] = 2 * u1[i] - u1Aux[i];
                u2_[i] = 2 * u2[i] - u2Aux[i];
            }


        }
        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                    "Error: %f\n", warpings, n, err_D);
    }

    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);

    free(v1);
    free(v2);

    free(rho_c);
    free(grad);

    free(u1_);
    free(u2_);

    free(u1Aux);
    free(u2Aux);

    free(I1x);
    free(I1y);

    free(I1w);
    free(I1wx);
    free(I1wy);

    free(div_xi1);
    free(div_xi2);

    free(u_N);
}

void tvl2OF(
        const float *I0,        // source image
        float *I1,              // target image
        float *u1,              // x component of the optical flow
        float *u2,              // y component of the optical flow
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float lambda,     // weight of the data term
        const float theta,      // weight of the data term
        const float tau,        // time step
        const float tol_OF,     // tol max allowed
        const int nx,           // image width
        const int ny,           // image height
        const int warps,        // number of warpings per scale
        const bool verbose      // enable/disable the verbose mode
) {
    using namespace std::chrono;
    auto clk_tvl2OF = system_clock::now();

    const float l_t = lambda * theta;
    const int size = nx * ny;


    auto *u1x = new float[size];
    auto *u1y = new float[size];
    auto *u2x = new float[size];
    auto *u2y = new float[size];

    auto *v1 = new float[size];
    auto *v2 = new float[size];

    auto *rho_c = new float[size];
    auto *grad = new float[size];

    auto *u1_ = new float[size];
    auto *u2_ = new float[size];

    auto *u1Aux = new float[size];
    auto *u2Aux = new float[size];

    auto *I1x = new float[size];
    auto *I1y = new float[size];

    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];

    // Divergence
    auto *div_xi1 = new float[size];
    auto *div_xi2 = new float[size];

    auto *u_N = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    auto clk_init_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_init = clk_init_end - clk_tvl2OF; // PROFILING
    cout << "(tvl2OF) initialising everything took "
         << elapsed_secs_init.count() << endl;

    double total_v1v2 = 0.0;
    double total_fwd_grad = 0.0;
    double total_getD = 0.0;
    double total_divergence = 0.0;
    double total_memcpy = 0.0;
    double total_getP = 0.0;
    double total_copy_u1u2 = 0.0;



    for (int warpings = 0; warpings < warps; warpings++) {
        //printf("warpings:%d\n", warpings);
	auto clk_warp_start = system_clock::now();
        // Compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1, u1, u2, I1w, nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

	auto clk_bicubic_end = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_bicubic = clk_bicubic_end - clk_warp_start; // PROFILING
        cout << "(tvl2OF) Bicubic interpolation took "
             << elapsed_secs_bicubic.count() << endl;

        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // Compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        memcpy(u1_, u1, size * sizeof(float));
        memcpy(u2_, u2, size * sizeof(float));
        // for (int i = 0; i < nx*ny; i++){
        //   u1_[i] = u1[i];
        //   u2_[i] = u2[i];
        // }

	auto clk_constants = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_constants = clk_constants - clk_bicubic_end; // PROFILING
        cout << "(tvl2OF) Computing constant part of functions & auxiliar variables (grad, rho_c, u1_, u2_) took "
             << elapsed_secs_constants.count() << endl;
/*
	double total_v1v2 = 0.0;
	double total_fwd_grad = 0.0;
	double total_getD = 0.0;
	double total_divergence = 0.0;
	double total_memcpy = 0.0;
	double total_getP = 0.0;
	double total_copy_u1u2 = 0.0;
*/
        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF * tol_OF && n < MAX_ITERATIONS_GLOBAL) {

            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
	    auto clk_v1v2 = system_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
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
#endif

	    auto clk_v1v2_end = system_clock::now(); // PROFILING
	    duration<double> elapsed_secs_v1v2 = clk_v1v2_end - clk_v1v2;  // PROFILING
	    total_v1v2 += elapsed_secs_v1v2.count();

            // Dual variables
            forward_gradient(u1_, u1x, u1y, nx, ny);
            forward_gradient(u2_, u2x, u2y, nx, ny);

            auto clk_fwd_grad_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_fwd_grad = clk_fwd_grad_end - clk_v1v2_end;  // PROFILING
            total_fwd_grad += elapsed_secs_fwd_grad.count();

            ofTVl2_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y, tau, size);

            auto clk_getD_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_getD = clk_getD_end - clk_fwd_grad_end;  // PROFILING
            total_getD += elapsed_secs_getD.count();

            // Primal variables
            divergence(xi11, xi12, div_xi1, nx, ny);
            divergence(xi21, xi22, div_xi2, nx, ny);

            auto clk_divergence_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_divergence = clk_divergence_end - clk_getD_end;  // PROFILING
            total_divergence += elapsed_secs_divergence.count();

            // Store previous iteration
            memcpy(u1Aux, u1, size * sizeof(float));
            memcpy(u2Aux, u2, size * sizeof(float));
            // for (int i = 0; i < size; i++){
            //   u1Aux[i] = u1[i];
            //   u2Aux[i] = u2[i];
            // }

            auto clk_memcpy_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_memcpy = clk_memcpy_end - clk_divergence_end;  // PROFILING
            total_memcpy += elapsed_secs_memcpy.count();

            ofTVl2_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N, theta, tau, size, &err_D);

            auto clk_getP_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_getP = clk_getP_end - clk_memcpy_end;  // PROFILING
            total_getP += elapsed_secs_getP.count();

            // (aceleration = 1);
            for (int i = 0; i < size; i++) {
                u1_[i] = 2 * u1[i] - u1Aux[i];
                u2_[i] = 2 * u2[i] - u2Aux[i];
            }

            auto clk_copy_u1u2_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_copy_u1u2 = clk_copy_u1u2_end - clk_getP_end;  // PROFILING
            total_copy_u1u2 += elapsed_secs_copy_u1u2.count();

        }

	auto clk_while_end = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_while = clk_while_end - clk_constants; // PROFILING
        cout << "(tvl2OF) While loop (with 400 it) took "
             << elapsed_secs_while.count() << endl;

        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                    "Error: %f\n", warpings, n, err_D);
	
	auto clk_warp_end = system_clock::now(); // PROFILING
   	duration<double> elapsed_secs_warp = clk_warp_end - clk_warp_start; // PROFILING
    	cout << "(tvl2OF) Warping num. " << warpings;
	cout << " took " << elapsed_secs_warp.count() << endl;

    }
    
    auto clk_all_warps_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_all_warps = clk_all_warps_end - clk_init_end; // PROFILING
    cout << "(tvl2OF) All warpings took "
         << elapsed_secs_all_warps.count() << endl;

    // PROFILING
    cout << "\n\n" << std::endl;
    cout << "Warpings loop profiling (total and %)" << endl;
    cout << "\t(v1-v2 loop) total: " << total_v1v2 << ", perc.: " <<
	    100 * (total_v1v2 / elapsed_secs_all_warps.count()) << "%" << endl;
    cout << "\t(fwd_grad) total: " << total_fwd_grad << ", perc.: " <<
            100 * (total_fwd_grad / elapsed_secs_all_warps.count()) << "%" << endl;
    cout << "\t(getD) total: " << total_getD << ", perc.: " <<
            100 * (total_getD / elapsed_secs_all_warps.count()) << "%" << endl; 
    cout << "\t(divergence) total: " << total_divergence << ", perc.: " <<
            100 * (total_divergence / elapsed_secs_all_warps.count()) << "%" << endl;
    cout << "\t(memcpy) total: " << total_memcpy << ", perc.: " <<
            100 * (total_memcpy / elapsed_secs_all_warps.count()) << "%" << endl;
    cout << "\t(getP) total: " << total_getP << ", perc.: " <<
            100 * (total_getP / elapsed_secs_all_warps.count()) << "%" << endl;
    cout << "\t(copy_u1u2) total: " << total_copy_u1u2 << ", perc.: " <<
            100 * (total_copy_u1u2 / elapsed_secs_all_warps.count()) << "%" << endl;

    cout << "\n\n" << std::endl;



    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);

    free(v1);
    free(v2);

    free(rho_c);
    free(grad);

    free(u1_);
    free(u2_);

    free(u1Aux);
    free(u2Aux);

    free(I1x);
    free(I1y);

    free(I1w);
    free(I1wx);
    free(I1wy);

    free(div_xi1);
    free(div_xi2);

    free(u_N);

    auto clk_tvl2OF_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_tvl2OF = clk_tvl2OF_end - clk_tvl2OF; // PROFILING
    cout << "(tvl2OF) All tasks took "
         << elapsed_secs_tvl2OF.count() << endl;

}

////////////////////////////////////NLTVL1//////////////////////////////////////
#define MAX_SPATIAL 2
#define MAX_INTENSITY 5
#define MAX_BETA  2 // Neighbour
#define MAX_DUAL_VAR (2*MAX_BETA + 1)*(2*MAX_BETA + 1)-1 // 5x5

// Define Dual Variables' struct
struct DualVariables_global {
    float sc[MAX_DUAL_VAR];     // value of p(x,y)
    float wp[MAX_DUAL_VAR];     // weight of non local
    int ap[MAX_DUAL_VAR];       // absolute position of p(y,x)
    int rp[MAX_DUAL_VAR];       // relative position of p(y,x) in the structure
    float wt = 0.0;
};

inline bool positive(int val) {
    return val >= 0;
}

float aux_pow2(float f) { return f * f; }

// Non-normalized images are assumed
void image_to_lab(float *in, int size, float *out) {
    const float T = 0.008856;
    const float color_attenuation = 1.5f;
    for (int i = 0; i < size; i++) {
        const float r = in[i] / 255.f;
        const float g = in[i + size] / 255.f;
        const float b = in[i + 2 * size] / 255.f;
        float X = 0.412453 * r + 0.357580 * g + 0.180423 * b;
        float Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
        float Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;
        X /= 0.950456;
        Z /= 1.088754;
        float Y3 = pow(Y, 1. / 3);
        float fX = X > T ? pow(X, 1. / 3) : 7.787 * X + 16 / 116.;
        float fY = Y > T ? Y3 : 7.787 * Y + 16 / 116.;
        float fZ = Z > T ? pow(Z, 1. / 3) : 7.787 * Z + 16 / 116.;
        float L = Y > T ? 116 * Y3 - 16.0 : 903.3 * Y;
        float A = 500 * (fX - fY);
        float B = 200 * (fY - fZ);

        // Correct L*a*b*: dark area or light area have less reliable colors
        float correct_lab = exp(-color_attenuation * aux_pow2(aux_pow2(L / 100) - 0.6));
        out[i] = L;
        out[i + size] = A * correct_lab;
        out[i + 2 * size] = B * correct_lab;
    }
}


static int validate_ap_2(int w, int h, int i, int j, int di, int dj) {
    const int r = j + dj;  // Row
    const int c = i + di;  // Column
    if (c < 0 || c >= w || r < 0 || r >= h)
        return -1;
    return r * w + c;
}

static float get_wspatial_2(int l, int k) {
    float ws = MAX_BETA;
    float w_tmp;
    float difS = 0.0;
    difS = hypot(l, k);
    w_tmp = exp(-difS / ws);
    //std::printf("Dif_S: %f W_S: %f\n",difS, w_tmp);

    return w_tmp;
}

static float get_wcolor_2(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
) {
    float wi = MAX_INTENSITY;
    float w_tmp;
    float difI = 0.0;
    for (int m = 0; m < pd; m++) {
        float aux = getsample_0(a, w, h, pd, i, j, m)
                    - getsample_0(a, w, h, pd, i + l, j + k, m);
        difI += aux * aux;
        //std::printf("valI_%d:%f ",m, difI);
    }
    //std::printf("\n");
    difI = sqrt(difI);
    w_tmp = exp(-difI / wi);

    return w_tmp;
}

static float get_weight_2(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
) {
    float wc = get_wcolor_2(a, w, h, i, j, l, k, pd);
    float ws = get_wspatial_2(l, k);

    return wc * ws;
}

void initialize_dual_variables(
        float *a,
        const int pd,
        const int w,
        const int h,
        const int n_d,
        const int radius,
        DualVariables_global *p,
        DualVariables_global *q
) {
    int size = w * h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++) {
            p[i].sc[j] = -2.0;
            p[i].wp[j] = -2.0;
            p[i].ap[j] = -1;  // Indicates that is out
            p[i].rp[j] = -1;  // Indicates that it is out.

            q[i].sc[j] = -2.0;
            q[i].wp[j] = -2.0;
            q[i].ap[j] = -1;  // Indicates that is out
            q[i].rp[j] = -1;  // Indicates that it is out.
        }


    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int it = 0;
            float ne = 0.0;
            const int pos = j * w + i;
            for (int k = -radius; k < (radius + 1); k++)
                for (int l = -radius; l < (radius + 1); l++) {
                    //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
                    if (!(k == 0 && k == l)) {

                        int ap = validate_ap_2(w, h, i, j, l, k);
                        if (positive(ap)) {
                            p[pos].sc[it] = q[pos].sc[it] = 0.0;
                            p[pos].ap[it] = q[pos].ap[it] = ap;
                            p[pos].rp[it] = q[pos].rp[it] = n_d - (it + 1);
                            assert(p[pos].rp[it] >= 0);
                            assert(q[pos].rp[it] >= 0);

                            // Compute the weight
                            float wp = sqrt(get_weight_2(a, w, h, i, j, l, k, pd));
                            p[pos].wp[it] = q[pos].wp[it] = wp;

                            ne += wp;
                        }
                        it++;
                        //std::printf("Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d Ne:%d\n",radius, j,i,k,l,it,ne);
                    }
                }
            // TODO: later used to normalize
            p[pos].wt = ne;
            q[pos].wt = ne;
        }
    // std::printf(" Ends\n");
}

void non_local_divergence(
        DualVariables_global *p,
        int size,
        int n_d,
        float *div_p
) {

    //#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        div_p[i] = 0.0;
        for (int j = 0; j < n_d; j++) {
            const int ap = p[i].ap[j];
            const int rp = p[i].rp[j];
            if (positive(ap)) {
                assert (p[i].rp[j] >= 0);
                const float pxy = p[i].sc[j];
                const float pyx = p[ap].sc[rp];
                const float w = p[i].wp[j];
                div_p[i] += w * (pxy - pyx);
            }
        }
        div_p[i] /= p[i].wt;
    }
}



// Auxiliar Chambolle Scheme functions

/*
 * - Name: getP

 *
*/
void ofnltv_getP(
        const float *v1,
        const float *v2,
        const float *div_p1,
        const float *div_p2,
        float theta,
        float tau,
        int size,
        float *u1,
        float *u2,
        float *err
) {
    float err_D = 0.0;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++) {

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k - tau * (div_p1[i] + (u1k - v1[i]) / theta);
        u2[i] = u2k - tau * (div_p2[i] + (u2k - v2[i]) / theta);

        err_D += (u1[i] - u1k) * (u1[i] - u1k) +
                 (u2[i] - u2k) * (u2[i] - u2k);
    }
    err_D /= size;
    (*err) = err_D;
#endif
}

/*
 * - Name: getD

 *
*/
void ofnltv_getD(
        const float *u1,
        const float *u2,
        int size,
        int n_d,
        float tau,
        DualVariables_global *p1,
        DualVariables_global *p2
) {
#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++) {
            const int ap1 = p1[i].ap[j];
            // const int rp1 = p1[i].rp[j];
            const float wt1 = p1[i].wt;

            const int ap2 = p2[i].ap[j];
            // const int rp2 = p2[i].rp[j];
            const float wt2 = p2[i].wt;

            if (positive(ap1) && positive(ap2)) {
                // assert(rp1 >=0);
                const float w1 = p1[i].wp[j];
                const float u1x = u1[i];
                const float u1y = u1[ap1];
                const float nlgr1 = w1 * (u1x - u1y) / wt1;
                const float nl1 = sqrt(nlgr1 * nlgr1);
                const float nl1g = 1 + tau * nl1;

                p1[i].sc[j] = (p1[i].sc[j] + tau * nlgr1) / nl1g;
            }

            if (positive(ap1) && positive(ap2)) {
                // assert(rp2 >=0);
                const float w2 = p2[i].wp[j];
                const float u2x = u2[i];
                const float u2y = u2[ap2];
                const float nlgr2 = w2 * (u2x - u2y) / wt2;
                const float nl2 = sqrt(nlgr2 * nlgr2);
                const float nl2g = 1 + tau * nl2;

                p2[i].sc[j] = (p2[i].sc[j] + tau * nlgr2) / nl2g;

            }
        }
#endif
}


void nltvl1_PD(
        const float *I0,              // source image
        float *I1,              // target image
        float *a,               // source image (color)
        int pd,                 // number of channels
        const float lambda,     // weight of the data term
        const float theta,      // weight of the data term
        const float tau,        // time step
        const int w,            // image width
        const int h,            // image height
        const int warps,        // number of warpings per scale
        const bool verbose,     // enable/disable the verbose mode
        float *u1,              // x component of the optical flow
        float *u2               // y component of the optical flow
) {
    const int size = w * h;
    const float l_t = lambda * theta;

    auto *p = new DualVariables_global[size];
    auto *q = new DualVariables_global[size];
    auto *v1 = new float[size];
    auto *v2 = new float[size];
    auto *rho_c = new float[size];
    auto *grad = new float[size];
    auto *u1_ = new float[size];
    auto *u2_ = new float[size];
    auto *u1_tmp = new float[size];
    auto *u2_tmp = new float[size];
    auto *I1x = new float[size];
    auto *I1y = new float[size];
    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];
    auto *div_p = new float[size];
    auto *div_q = new float[size];

    int radius = MAX_BETA;
    int n_d = MAX_DUAL_VAR;


    std::printf("Before\n");
    // Initialization of the Dual variables.
    initialize_dual_variables(a, pd, w, h, n_d, radius, p, q);
    centered_gradient(I1, I1x, I1y, w, h);

    std::printf("Initialization\n");
    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1, u1, u2, I1w, w, h, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, w, h, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, w, h, true);
        //#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // Compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        for (int i = 0; i < size; i++) {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        // while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        while (n < MAX_ITERATIONS_GLOBAL) {

            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
#ifdef _OPENMP
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
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
#endif
            // Dual variables
            ofnltv_getD(u1_, u2_, size, n_d, tau, p, q);
            // Store the previous iteration
            for (int i = 0; i < size; i++) {
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            // Primal variables
            non_local_divergence(p, size, n_d, div_p);
            non_local_divergence(q, size, n_d, div_q);
            ofnltv_getP(v1, v2, div_p, div_q, theta, tau, size, u1, u2, &err_D);

            // (acceleration = 1);
            for (int i = 0; i < size; i++) {
                u1_[i] = 2 * u1[i] - u1_tmp[i];
                u2_[i] = 2 * u2[i] - u2_tmp[i];
            }

        }
        if (verbose)
            std::printf("Warping: %d,Iter: %d Error: %f\n", warpings, n, err_D);
    }

    delete[] v1;
    delete[] v2;

    delete[] rho_c;
    delete[] grad;

    delete[] u1_;
    delete[] u2_;

    delete[] u1_tmp;
    delete[] u2_tmp;

    delete[] I1x;
    delete[] I1y;

    delete[] I1w;
    delete[] I1wx;
    delete[] I1wy;

    delete[] div_p;
    delete[] div_q;
}
//////////////////////////////TV-CSAD///////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void initialize_pos_nei(
        const int w,
        const int h,
        const int n_d,
        const int radius,
        PosNei *p
) {
    int size = w * h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++) {
            p[i].api[j] = -1;  // Indicates that is out
            p[i].apj[j] = -1;
            p[i].b[j] = 0.0;
        }


    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int it = 0;
            float ne = 0.0;
            const int pos = j * w + i;
            for (int k = -radius; k < (radius + 1); k++)
                for (int l = -radius; l < (radius + 1); l++) {
                    //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
                    if (!(k == 0 && k == l)) {

                        int ap = validate_ap_2(w, h, i, j, l, k);
                        if (positive(ap)) {
                            p[pos].api[it] = i + l;
                            p[pos].apj[it] = j + k;
                            assert(p[pos].api[it] >= 0 && p[pos].apj[it] >= 0);
                            ne++;
                        }
                        it++;
                        //std::printf("Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d Ne:%d\n",radius, j,i,k,l,it,ne);
                    }
                }
            p[pos].n = ne;
            for (int k = 0; k < (2 * p[pos].n + 1); k++) {
                p[pos].ba.push_back(0.0);
            }
        }
}


//Auxiliar Chambolle Scheme functions

//Chambolle functions
/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void tvcsad_getP(float *u1,
                 float *u2,
                 const float *v1,
                 const float *v2,
                 const float *div_xi1,
                 const float *div_xi2,
                 float theta,
                 float tau,
                 int size,
                 float *err) {
    float err_D = 0.0;

#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++) {

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k - tau * (-div_xi1[i] + (u1k - v1[i]) / theta);
        u2[i] = u2k - tau * (-div_xi2[i] + (u2k - v2[i]) / theta);

        err_D += (u1[i] - u1k) * (u1[i] - u1k) +
                 (u2[i] - u2k) * (u2[i] - u2k);

        // u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
        //   (u2[i] - u2k) * (u2[i] - u2k);
    }
#endif

    // getminmax(&min,&max,u_N,size);

    // err_D =max;
    err_D /= size;
    (*err) = err_D;
}

/*
 * - Name: getD_Du

 *
*/
void tvcsad_getD(float *xi11, float *xi12, float *xi21, float *xi22, float *u1x, float *u1y, float *u2x, float *u2y,
                 float tau, int size) {
#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        float xi1_N = hypot(xi11[i], xi12[i]);
        float xi2_N = hypot(xi21[i], xi22[i]);

        xi1_N = MAX(1, xi1_N);
        xi2_N = MAX(1, xi2_N);

        xi11[i] = (xi11[i] + tau * u1x[i]) / xi1_N;
        xi12[i] = (xi12[i] + tau * u1y[i]) / xi1_N;

        xi21[i] = (xi21[i] + tau * u2x[i]) / xi2_N;
        xi22[i] = (xi22[i] + tau * u2y[i]) / xi2_N;
    }
#endif
}


void tvcsad_PD(
        const float *I0,            // source image
        float *I1,                  // target image
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float lambda,         // weight of the data term
        const float theta,          // weight of the data term
        const float tau,            // time step
        const float tol_OF,         // tol max allowed
        const int nx,               // image width
        const int ny,               // image height
        const int warps,            // number of warpings per scale
        const bool verbose,         // enable/disable the verbose mode
        float *u1,                  // x component of the optical flow
        float *u2                   // y component of the optical flow
) {

    const float l_t = lambda * theta;
    const int size = nx * ny;
    const int n_d = DT_NEI;
    const int r = DT_R;
    auto *p = new PosNei[size];

    auto *u1x = new float[size];
    auto *u1y = new float[size];
    auto *u2x = new float[size];
    auto *u2y = new float[size];

    auto *v1 = new float[size];
    auto *v2 = new float[size];

    auto *rho_c = new float[size];
    auto *grad = new float[size];

    auto *u1_ = new float[size];
    auto *u2_ = new float[size];

    auto *u1_tmp = new float[size];
    auto *u2_tmp = new float[size];

    auto *I1x = new float[size];
    auto *I1y = new float[size];

    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];

    // Divergence
    auto *div_xi1 = new float[size];
    auto *div_xi2 = new float[size];

    auto *u_N = new float[size];

    // Five point gradient of the right,left view. (1/12)*[-1 8 0 -8 1]
    centered_gradient(I1, I1x, I1y, nx, ny);
    initialize_pos_nei(nx, ny, n_d, r, p);

    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1, u1, u2, I1w, nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);
        // #pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1(p + u))| (Warping image)
            grad[i] = hypot(Ix2 + Iy2, 0.01);

            for (int j = 0; j < n_d; j++) {
                // std::printf("I:%d Iter:%d J:%d I:%d \n",i, j,  p[i].apj[j], p[i].api[j]);
                if (positive(p[i].api[j]) && positive(p[i].apj[j])) {
                    // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, p[i].apj[j]*nx + p[i].api[j], p[i].apj[j], p[i].api[j]);
                    assert(p[i].api[j] >= 0);
                    assert(p[i].apj[j] >= 0);
                    assert(p[i].apj[j] * nx + p[i].api[j] < nx * ny);
                    const int pos = p[i].apj[j] * nx + p[i].api[j];

                    p[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                                 + I1wy[i] * u2[i]) / grad[i];
                }
            }
        }

        for (int i = 0; i < nx * ny; i++) {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF * tol_OF && n < MAX_ITERATIONS_GLOBAL) {
            n++;
            // Estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
#ifdef _OPENMP
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                int it = 0;
                for (int j = 0; j < n_d; j++) {
                    if (positive(p[i].api[j]) && positive(p[i].apj[j])) {
                        // std::printf("J:%d I:%d\n",p[i].api[j],p[i].apj[j]);
                        p[i].ba[it] = -(p[i].b[j] - (I1wx[i] * u1[i]
                                                     + I1wy[i] * u2[i]) / grad[i]);
                        it++;
                    }
                }
                for (int j = 0; j < (p[i].n + 1); j++) {
                    p[i].ba[it] = (p[i].n - 2 * j) * l_t * grad[i];
                    it++;
                }

                std::sort(p[i].ba.begin(), p[i].ba.begin() + it);
                // v1[i] = u1[i] - l_t*I1wx[i]*p[i].ba[it/2+1]/grad[i];
                // v2[i] = u2[i] - l_t*I1wy[i]*p[i].ba[it/2+1]/grad[i];
                //TODO: possible error in the minimization
                v1[i] = u1[i] - I1wx[i] * p[i].ba[it / 2 + 1] / grad[i];
                v2[i] = u2[i] - I1wy[i] * p[i].ba[it / 2 + 1] / grad[i];
            }
#endif
            // Data term

            // Dual variables
            forward_gradient(u1_, u1x, u1y, nx, ny);
            forward_gradient(u2_, u2x, u2y, nx, ny);
            tvcsad_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y, tau, size);

            // Primal variables
            divergence(xi11, xi12, div_xi1, nx, ny);
            divergence(xi21, xi22, div_xi2, nx, ny);

            // Store previous iteration
            for (int i = 0; i < size; i++) {
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            tvcsad_getP(u1, u2, v1, v2, div_xi1, div_xi2, theta, tau, size, &err_D);

            // (acceleration = 1);
            for (int i = 0; i < size; i++) {
                u1_[i] = 2 * u1[i] - u1_tmp[i];
                u2_[i] = 2 * u2[i] - u2_tmp[i];

            }


        }
        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                    "Error: %f\n", warpings, n, err_D);
    }

    delete[] p;

    delete[] u1x;
    delete[] u1y;
    delete[] u2x;
    delete[] u2y;

    delete[] v1;
    delete[] v2;

    delete[] rho_c;
    delete[] grad;

    delete[] u1_;
    delete[] u2_;

    delete[] u1_tmp;
    delete[] u2_tmp;

    delete[] I1x;
    delete[] I1y;

    delete[] I1w;
    delete[] I1wx;
    delete[] I1wy;

    delete[] div_xi1;
    delete[] div_xi2;

    delete[] u_N;
    std::printf("Exits current level\n");

}

///////////////////////NLTV-CSAD///////////////////////


void nltvcsad_PD(
        const float *I0,              // source image
        float *I1,              // target image
        float *a,               // source image (color)
        int pd,                 // number of channels
        const float lambda,     // weight of the data term
        const float theta,      // weight of the data term
        const float tau,        // time step
        const int w,            // image width
        const int h,            // image height
        const int warps,        // number of warpings per scale
        const bool verbose,     // enable/disable the verbose mode
        float *u1,              // x component of the optical flow
        float *u2               // y component of the optical flow
) {

    const int size = w * h;
    const float l_t = lambda * theta;

    auto *p = new DualVariables_global[size];
    auto *q = new DualVariables_global[size];
    auto *pnei = new PosNei[size];
    auto *v1 = new float[size];
    auto *v2 = new float[size];
    auto *rho_c = new float[size];
    auto *grad = new float[size];
    auto *u1_ = new float[size];
    auto *u2_ = new float[size];
    auto *u1_tmp = new float[size];
    auto *u2_tmp = new float[size];
    auto *I1x = new float[size];
    auto *I1y = new float[size];
    auto *I1w = new float[size];
    auto *I1wx = new float[size];
    auto *I1wy = new float[size];
    auto *div_p = new float[size];
    auto *div_q = new float[size];

    int radius = MAX_BETA;
    int n_d = MAX_DUAL_VAR;

    const int ndt = DT_NEI;
    const int rdt = DT_R;


    // Initialization of the Dual variables.
    initialize_dual_variables(a, pd, w, h, n_d, radius, p, q);
    initialize_pos_nei(w, h, ndt, rdt, pnei);
    centered_gradient(I1, I1x, I1y, w, h);

    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1, u1, u2, I1w, w, h, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, w, h, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, w, h, true);
        //#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // Store the |Grad(I1(p + u))| (Warping image)
            grad[i] = Ix2 + Iy2;
            if (grad[i] > GRAD_IS_ZERO) {
                for (int j = 0; j < ndt; j++) {
                    // std::printf("I:%d Iter:%d J:%d I:%d \n",i, j,  p[i].apj[j], p[i].api[j]);
                    if (positive(pnei[i].api[j]) && positive(pnei[i].apj[j])) {
                        // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, p[i].apj[j]*nx + p[i].api[j], p[i].apj[j], p[i].api[j]);
                        assert(pnei[i].api[j] >= 0);
                        assert(pnei[i].apj[j] >= 0);
                        assert(pnei[i].apj[j] * w + pnei[i].api[j] < w * h);
                        const int pos = pnei[i].apj[j] * w + pnei[i].api[j];

                        pnei[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                                        + I1wy[i] * u2[i]) / sqrt(grad[i]);
                    }
                }
            }
        }

        for (int i = 0; i < size; i++) {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        // while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        while (n < MAX_ITERATIONS_GLOBAL) {
            n++;
            // Estimate the values of the variable (v1, v2)
#ifdef _OPENMP
#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                v1[i] = u1[i];
                v2[i] = u2[i];
                if (grad[i] > GRAD_IS_ZERO) {
                    int it = 0;
                    for (int j = 0; j < ndt; j++) {
                        if (positive(pnei[i].api[j]) && positive(pnei[i].apj[j])) {
                            pnei[i].ba[it] = -(pnei[i].b[j] - (I1wx[i] * u1[i]
                                                               + I1wy[i] * u2[i]) / sqrt(grad[i]));
                            it++;
                        }
                    }
                    for (int j = 0; j < (pnei[i].n + 1); j++) {
                        pnei[i].ba[it] = (pnei[i].n - 2 * j) * l_t * sqrt(grad[i]);
                        it++;
                    }

                    std::sort(pnei[i].ba.begin(), pnei[i].ba.begin() + it);
                    // v1[i] = u1[i] - l_t*I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
                    // v2[i] = u2[i] - l_t*I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
                    // TODO: possible error in the minimization
                    v1[i] = u1[i] - I1wx[i] * pnei[i].ba[it / 2 + 1] / sqrt(grad[i]);
                    v2[i] = u2[i] - I1wy[i] * pnei[i].ba[it / 2 + 1] / sqrt(grad[i]);
                }
            }
#endif
            // Dual variables
            ofnltv_getD(u1_, u2_, size, n_d, tau, p, q);
            // Store previous iteration
            for (int i = 0; i < size; i++) {
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            // Primal variables
            non_local_divergence(p, size, n_d, div_p);
            non_local_divergence(q, size, n_d, div_q);
            ofnltv_getP(v1, v2, div_p, div_q, theta, tau, size, u1, u2, &err_D);
            // (acceleration = 1);
            for (int i = 0; i < size; i++) {
                u1_[i] = 2 * u1[i] - u1_tmp[i];
                u2_[i] = 2 * u2[i] - u2_tmp[i];
            }

        }
        if (verbose)
            std::printf("Warping: %d,Iter: %d Error: %f\n", warpings, n, err_D);
    }

    delete[] p;
    delete[] q;
    delete[] pnei;

    delete[] v1;
    delete[] v2;

    delete[] rho_c;
    delete[] grad;

    delete[] u1_;
    delete[] u2_;

    delete[] u1_tmp;
    delete[] u2_tmp;

    delete[] I1x;
    delete[] I1y;

    delete[] I1w;
    delete[] I1wx;
    delete[] I1wy;

    delete[] div_p;
    delete[] div_q;
}








/////////////////////////MAIN/////////////////////


void rgb2gray(float *in, int w, int h, float *out) {
    int size = w * h;
    for (int i = 0; i < size; i++) {
        out[i] = .299 * in[i] + .587 * in[size + i] + .114 * in[2 * size + i];

    }

}


/*
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nprocs      number of threads to use (OpenMP library)
 *   -I0          first image
 *   -I1          second image
 *   -tau         time step in the numerical scheme
 *   -theta       attachment parameter between E_Data and E_Smooth
 *   -nwarps      number of warps per scales
 *   -out         name of the output flow field
 *   -verbose     switch on/off messages
 *
 */

int main(int argc, char *argv[]) {

    using namespace std::chrono;
    system_clock::time_point today = system_clock::now();
    time_t tt;

    tt = system_clock::to_time_t(today);
    std::cerr << "today is: " << ctime(&tt);

    // process input
    std::vector<std::string> args(argv, argv + argc);
    auto warps_val = pick_option(args, "w",
                                 to_string(PAR_DEFAULT_NWARPS_GLOBAL));     // Warpings
    auto var_reg = pick_option(args, "m", to_string(M_TVL1));               // Methods
    auto file_params = pick_option(args, "p", "");                          // Params' file
    auto global_iters = pick_option(args, "glb_iters",
                                    to_string(MAX_ITERATIONS_GLOBAL));      // Faldoi global iterations

    if (args.size() != 6 && args.size() != 4) {
        fprintf(stderr, "Without occlusions:\n");
        fprintf(stderr, "Usage: %lu  ims.txt in_flow.flo  out.flo "
                "[-m method_val] [-w num_warps] [-p file of parameters] val [-glb_iters global_iters] \n", args.size());
        fprintf(stderr, "With occlusions:\n");
        fprintf(stderr, "Usage: %lu  ims.txt in_flow.flo  out.flo occl_input.png occl_out.png"
                " [-m method_val] [-w num_warps] [-p file of parameters] val [-glb_iters global_iters] \n", args.size());

        return EXIT_FAILURE;
    }

    // O TV-l2 coupled 1 - ||Du + Du'||_{F}
    // Optional input arguments
    int val_method = stoi(var_reg);
    int nwarps = stoi(warps_val);
    int glb_it = stoi(global_iters);


    // Read the parameters
    // Filename that contains the paths to the images to use
    const std::string &filename_images = args[1];
    const std::string &image_flow_name = args[2];
    const std::string &outfile = args[3];
    // Initialize occlusion-specific params
    std::string occ_input;
    std::string occ_output;

    //TODO: works with occlusions but we lost 'const' advantages
    // being constant references this may only be possible if the const is applied
    // after somehow...(??)
    if (args.size() == 6) { // if occlusions
        occ_input = args[4];
        occ_output = argv[5];
    }


    // Filename of images
    std::string filename_i_1, filename_i0, filename_i1;

    // Read txt file of images

    ifstream infile(filename_images);
    int num_files = 0;
    string line;
    while (getline(infile, line)) {

        ++num_files;
        if (num_files == 3) {
            filename_i_1 = line;
        } else {
            if (num_files == 1) {
                filename_i0 = line;
            } else {
                if (num_files == 2) {
                    filename_i1 = line;
                }
            }
        }
    }
    infile.close();


    // Open input images
    int w[5], h[5], pd[5];
    float *i_1;
    if (num_files == 4) {
        i_1 = iio_read_image_float_split(filename_i_1.c_str(), w + 3, h + 3, pd + 3);
    } else {
        i_1 = iio_read_image_float_split(filename_i1.c_str(), w + 3, h + 3, pd + 3);
    }

    float *i0 = iio_read_image_float_split(filename_i0.c_str(), w + 0, h + 0, pd + 0);
    float *i1 = iio_read_image_float_split(filename_i1.c_str(), w + 1, h + 1, pd + 1);
    float *flow = iio_read_image_float_split(image_flow_name.c_str(), w + 2, h + 2, pd + 2);

    float *occ = nullptr;
    if (val_method >= 8) {
        occ = iio_read_image_float_split(occ_input.c_str(), w + 4, h + 4, pd + 4);
    }

    // Ensure that dimensions match
    if (num_files == 3) {
        if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
        if (w[0] != w[3] || h[0] != h[3] || pd[0] != pd[3])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
        if (w[1] != w[3] || h[1] != h[3] || pd[1] != pd[3])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
    } else {
        if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
    }
    // Ensure dimensions match between flow and images
    if (w[0] != w[2] || h[0] != h[2] || pd[2] != 2)
        return fprintf(stderr, "ERROR: input flow field size mismatch\n");

    // Print method used
    if (num_files == 2 && val_method == M_TVL1_OCC) {
        // If only two images given for occ, revert to tvl1 without occlusions
        // TODO: when new methods with occlusions implemented, add here
        switch (val_method) {
            case M_TVL1_OCC:
                fprintf(stderr, "Since only two images given, method is changed to TV-l2 coupled\n");
                fprintf(stderr, "Occlusion estimation requires 4 frames: i_1 ==> i0 ==> i1 ==> i2\n");
                val_method = M_TVL1;
                break;
            default:
                fprintf(stderr, "Method unknown\n");
                break;
        }
    } else {
        // If four images given for without occ, two not needed
        if (num_files == 4 && val_method >= 0 && val_method <= 7) {
            fprintf(stderr, "Only two of the four images given will be used, according to the method selected\n");
            fprintf(stderr, "Method: ");
            switch (val_method) {
                case M_NLTVL1:      // NLTVL1
                    fprintf(stderr, "NLTV-L1\n");
                    break;
                case M_TVCSAD:      // TV-CSAD
                    fprintf(stderr, "TV-CSAD\n");
                    break;
                case M_NLTVCSAD:    // NLTV-CSAD
                    fprintf(stderr, "NLTV-CSAD\n");
                    break;
                case M_TVL1_W:      // TV-l2 with weights
                    fprintf(stderr, "TV-l2 coupled Weights\n");
                    break;
                case M_NLTVCSAD_W:  // NLTV-CSAD with weights
                    fprintf(stderr, "NLTV-CSAD Weights\n");
                    break;
                case M_NLTVL1_W:    // NLTV-L1 with weights
                    fprintf(stderr, " NLTV-L1 Weights\n");
                    break;
                case M_TVCSAD_W:    // TV-CSAD with weights
                    fprintf(stderr, "TV-CSAD Weights\n");
                    break;
                default:            // TV-l2 coupled
                    fprintf(stderr, "TV-l2 coupled\n");
            }
        } else {
            fprintf(stderr, "Method: ");
            switch (val_method) {
                case M_TVL1_OCC:    // TV-l2 with occlusion
                    fprintf(stderr, "TV-l2 occlusions\n");
                    break;
                default:
                    break;
            }
        }
    }

    // Initialize parameters
    int step_alg = GLOBAL_STEP;
    Parameters params = init_params(file_params, step_alg);
    params.w = w[0];
    params.h = h[0];
    params.warps = nwarps;
    params.val_method = val_method;
    params.iterations_of = glb_it;
    if (params.verbose)
        cerr << params;

    OpticalFlowData ofD = init_Optical_Flow_Data(params);

    float *a = nullptr;
    float *xi11 = nullptr;
    float *xi12 = nullptr;
    float *xi21 = nullptr;
    float *xi22 = nullptr;

    int size = w[0] * h[0];
    // 0 - TVl2 coupled, otherwise Du
    if (val_method == M_NLTVL1 || val_method == M_NLTVL1_W || val_method == M_NLTVCSAD || val_method == M_NLTVCSAD_W) {
        //printf("NL-TVL1 or NLTV-CSAD\n");
        std::printf("W:%d H:%d Pd:%d\n", w[0], h[0], pd[0]);
        a = new float[size * pd[0]];
        image_to_lab(i0, size, a);
    }

    auto *i0n = new float[size];
    auto *i1n = new float[size];
    auto *i_1n = new float[size];

    if (pd[0] != 1) {

        rgb2gray(i0, w[0], h[0], i0n);
        rgb2gray(i1, w[0], h[0], i1n);
        rgb2gray(i_1, w[0], h[0], i_1n);

    } else {

        memcpy(i0n, i0, size * sizeof(float));
        memcpy(i1n, i1, size * sizeof(float));
        memcpy(i_1n, i_1, size * sizeof(float));
    }
    image_normalization_3(i0n, i1n, i_1n, i0n, i1n, i_1n, size);
    gaussian(i0n, w[0], h[0], PRESMOOTHING_SIGMA);
    gaussian(i1n, w[0], h[0], PRESMOOTHING_SIGMA);
    gaussian(i_1n, w[0], h[0], PRESMOOTHING_SIGMA);



    // Allocate memory for the flow
    float *u = ofD.u1;
    float *v = ofD.u2;
    float *chi = ofD.chi;

    auto clk_init_start = system_clock::now();
    // Initialize flow with flow from local faldoi
    for (int i = 0; i < size; i++) {
        u[i] = flow[i];
        v[i] = flow[size + i];
        if (val_method >= 8) {
            chi[i] = occ[i];
        }
    }

    SpecificOFStuff stuffOF{};
    PatchIndexes index{};
    // Initialize dual variables if necessary (TV) and other stuff for TVL1_occ
    if (val_method == M_TVL1 || val_method == M_TVL1_W || val_method == M_TVCSAD || val_method == M_TVCSAD_W
        || val_method == M_TVL1_OCC) {

        if (val_method == M_TVL1_OCC) {

            index.ii = 0;
            index.ij = 0;
            index.ei = w[0];
            index.ej = h[0];


            initialize_auxiliar_stuff(stuffOF, ofD, params.w, params.h);
            // Derivatives of I0 to compute weight g

            xi11 = stuffOF.tvl2_occ.xi11;
            xi12 = stuffOF.tvl2_occ.xi12;
            xi21 = stuffOF.tvl2_occ.xi21;
            xi22 = stuffOF.tvl2_occ.xi22;
        } else {
            xi11 = new float[size];
            xi12 = new float[size];
            xi21 = new float[size];
            xi22 = new float[size];
        }


        for (int i = 0; i < size; i++) {
            xi11[i] = 0.0;
            xi12[i] = 0.0;
            xi21[i] = 0.0;
            xi22[i] = 0.0;
        }
    }

    auto clk_init_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_init = clk_init_end - clk_init_start; // PROFILING
    cout << "(global_faldoi.cpp) initialising everything took "
         << elapsed_secs_init.count() << endl;

    // 0 - TVl2 coupled, otherwise Du
    if (val_method == M_TVL1 || val_method == M_TVL1_W) {
        //printf("TV-l2 coupled\n");
        tvl2OF(i0n, i1n, u, v, xi11, xi12, xi21, xi22, params.lambda, params.theta, params.tau, params.tol_OF, w[0],
               h[0], params.warps, params.verbose);

    } else if (val_method == M_NLTVCSAD || val_method == M_NLTVCSAD_W) {
        params.lambda = 0.85;
        params.theta = 0.3;
        params.tau = 0.1;
        //printf("NLTV-CSAD\n");
        nltvcsad_PD(i0n, i1n, a, pd[0], params.lambda, params.theta, params.tau, w[0], h[0], params.warps,
                    params.verbose, u, v);

    } else if (val_method == M_NLTVL1 || val_method == M_NLTVL1_W) {
        params.lambda = 2.0;
        params.theta = 0.3;
        params.tau = 0.1;
        //printf("NLTV-L1\n");
        nltvl1_PD(i0n, i1n, a, pd[0], params.lambda, params.theta, params.tau, w[0], h[0], params.warps,
                  params.verbose, u, v);

    } else if (val_method == M_TVCSAD || val_method == M_TVCSAD_W) {
        params.lambda = 0.85;
        params.theta = 0.3;
        params.tau = 0.125;
        //printf("TV-CSAD\n");
        tvcsad_PD(i0n, i1n, xi11, xi12, xi21, xi22, params.lambda, params.theta, params.tau, params.tol_OF, w[0], h[0],
                  params.warps, params.verbose, u, v);

    } else if (val_method == M_TVL1_OCC) {
        //fprintf(stderr, "TV-l2 occlusions\n");
        float ener_N;

        guided_tvl2coupled_occ(i0n, i1n, i_1n, &ofD, &(stuffOF.tvl2_occ), &ener_N, index, params.w, params.h);

    }

    iio_save_image_float_split(outfile.c_str(), u, w[0], h[0], 2);


    if (val_method == M_TVL1_OCC) {
        auto *out_occ_int = new int[w[0] * h[0]];
        for (int i = 0; i < w[0] * h[0]; i++) {

            out_occ_int[i] = chi[i];
        }
        iio_save_image_int(occ_output.c_str(), out_occ_int, w[0], h[0]);
        delete[] out_occ_int;
        free_auxiliar_stuff(&stuffOF, &ofD);
    }

    // Delete allocated memory

    delete[] u;
    delete[] chi;
    if (val_method == M_TVL1 || val_method == M_TVL1_W || val_method == M_TVCSAD || val_method == M_TVCSAD_W) {
        delete[] xi11;
        delete[] xi12;
        delete[] xi21;
        delete[] xi22;
    } else {
        delete[] a;
    }

    delete[] i0n;
    delete[] i1n;
    delete[] i_1n;
    today = system_clock::now();

    tt = system_clock::to_time_t(today);
    std::cerr << "today is: " << ctime(&tt);
    return EXIT_SUCCESS;

    auto clk_global_min_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_global_min = clk_global_min_end - clk_init_end; // PROFILING
    cout << "(global_faldoi.cpp) global minimisation (functional-specific) took "
         << elapsed_secs_global_min.count() << endl;

}

#endif//GLOBAL_FALDOI
