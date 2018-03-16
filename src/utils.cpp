// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING


#include <vector>
#include <cmath>
#include <cassert>
#include <cstdio>
#include "energy_structures.h"

extern "C" {
#include "bicubic_interpolation.h"
#include "xmalloc.h"
}


#define PAR_DEFAULT_GAMMA 0.05
PatchIndexes get_index_patch(
        const int wr,
        const int w,
        const int h,
        const int i,
        const int j,
        const int factor
        ) {
    PatchIndexes index;
    //Points to begin and end. End is the previous value
    index.i = i;
    index.j = j;
    index.ii = ((i - factor * wr) < 0)? 0 : (i - factor * wr);
    index.ij = ((j - factor * wr) < 0)? 0 : (j - factor * wr);
    index.ei = ((i + 1 + factor * wr) > w)? w : (i + 1 + factor * wr);
    index.ej = ((j + 1 + factor * wr) > h)? h : (j + 1 + factor * wr);
    return index;

}


///////////////////////////////////////////
/////////DERIVATIVES PATCH/////////////////
///////////////////////////////////////////


//Compute the divergence (backward differences) over a patch
void divergence_patch(
        const float *v1, // x component of the vector field
        const float *v2, // y component of the vector field
        float *div,      // output divergence
        const int ii,     // initial column
        const int ij,     // initial row
        const int ei,     // end column
        const int ej,     // end row
        const int nx    // image width
        ){

    // compute the divergence on the central body of the image
//#pragma omp simd collapse(2)
//#pragma omp for schedule(dynamic, 1) collapse(2)
    for (int j = ij + 1; j < ej - 1; j++){
        for (int i = ii + 1; i < ei - 1; i++){
            const int p  = j * nx + i;
            const int p1 = p - 1;
            const int p2 = p - nx;

            const float v1x = v1[p] - v1[p1];
            const float v2y = v2[p] - v2[p2];

            div[p] = v1x + v2y;
        }
    }
    // compute the divergence on the first and last rows
    for (int i = ii + 1; i < ei - 1; i++) {
        const int p = (ej - 1) * nx + i;

        div[i] = v1[i] - v1[i-1] + v2[i];
        div[p] = v1[p] - v1[p-1] - v2[p-nx];
    }

    // compute the divergence on the first and last columns
    for (int j = ij + 1; j < ej - 1; j++) {
        const int p1 = j * nx;
        const int p2 = (j + 1) * nx - 1;

        div[p1] =  v1[p1]     + v2[p1] - v2[p1 - nx];
        div[p2] = -v1[p2 - 1] + v2[p2] - v2[p2 - nx];

    }
    div[ij*nx + ii]           =  v1[ij*nx + ii]          + v2[ij*nx + ii];
    //div[ei - 1 + ii]          = -v1[ei + ii  - 2]        + v2[ei + ii - 1];
    div[ij*nx + ei - 1]       = -v1[ij*nx + ei - 2]      + v2[ij*nx + ei - 1];
    div[(ej - 1)*nx + ii]     =  v1[(ej - 1)*nx + ii]    - v2[(ej - 2)*nx + ii];
    div[(ej - 1)*nx + ei - 1] = -v1[(ej - 1)*nx + ei -2] - v2[(ej - 2)*nx + ei -1];

}


//Compute the forward gradient (forward difference) over a patch
void forward_gradient_mixed_bound_patch(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int ii,     // initial column
        const int ij,     // initial row
        const int ei,     // end column
        const int ej,     // end row
        const int nx,   //image width
        const int ny   //image height
        ){

    // compute the divergence on the central body of the image
//#pragma omp simd collapse(2)
//#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int j = ij; j < ej - 1; j++){
        for (int i = ii; i < ei - 1; i++){
            const int p = j*nx + i;
            const int p1 = p + 1;
            const int p2 = p + nx;

            fx[p] = f[p1] - f[p];
            fy[p] = f[p2] - f[p];
        }
    }

    // compute the gradient on the last row
    for (int i = ii; i < ei - 1; i++) {
        const int p = (ej-1) * nx + i;


        fx[p] = f[p + 1] - f[p];
        if (ej == ny){
            fy[p] = 0;
        }else{
            fy[p] = f[p + nx] - f[p];
        }

    }

    // compute the gradient on the last column
    for (int j = ij; j < ej-1; j++) {
        const int p = j*nx + ei -1;

        if (ei == nx){
            fx[p] = 0;
        }else{
            fx[p] = f[p + 1] - f[p];
        }
        fy[p] = f[p+nx] - f[p];
    }

    fx[(ej - 1)*nx + ei - 1] = 0;
    fy[(ej - 1)*nx + ei - 1] = 0;
}



//Compute the forward gradient (forward difference) over a patch
void forward_gradient_patch(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int ii,     // initial column
        const int ij,     // initial row
        const int ei,     // end column
        const int ej,     // end row
        const int nx   //image width
        ){

    // compute the divergence on the central body of the image
//#pragma omp simd collapse(2)
//#pragma omp for schedule(dynamic, 1) collapse(2)
    for (int j = ij; j < ej-1; j++){
        for (int i = ii; i < ei-1; i++){
            const int p = j*nx + i;
            const int p1 = p + 1;
            const int p2 = p + nx;

            fx[p] = f[p1] - f[p];
            fy[p] = f[p2] - f[p];
        }
    }

    // compute the gradient on the last row
    for (int i = ii; i < ei - 1; i++){
        const int p = (ej - 1) * nx + i;

        fx[p] = f[p + 1] - f[p];
        fy[p] = 0;
    }

    // compute the gradient on the last column
    for (int j = ij; j < ej-1; j++){

        const int p = j*nx + ei -1;

        fx[p] = 0;
        fy[p] = f[p + nx] - f[p];
    }

    fx[(ej - 1)*nx + ei - 1] = 0;
    fy[(ej - 1)*nx + ei - 1] = 0;
}






/////////////////////////////////////////
///////////DERIVATIVES IMAGE////////////
////////////////////////////////////////

/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/

void divergence(
        const float *v1, // x component of the vector field
        const float *v2, // y component of the vector field
        float *div,      // output divergence
        const int nx,    // image width
        const int ny     // image height
        ) {
    // compute the divergence on the central body of the image
//#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < ny - 1; i++) {
        for(int j = 1; j < nx - 1; j++){
            const int p  = i * nx + j;
            const int p1 = p - 1;
            const int p2 = p - nx;

            const float v1x = v1[p] - v1[p1];
            const float v2y = v2[p] - v2[p2];

            div[p] = v1x + v2y;
        }
    }

    // compute the divergence on the first and last rows
    for (int j = 1; j < nx - 1; j++){
        const int p = (ny - 1) * nx + j;

        div[j] = v1[j] - v1[j - 1] + v2[j];
        div[p] = v1[p] - v1[p - 1] - v2[p - nx];
    }

    // compute the divergence on the first and last columns
    for (int i = 1; i < ny - 1; i++){
        const int p1 = i * nx;
        const int p2 = (i+1) * nx - 1;

        div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
        div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];

    }

    div[0]         =  v1[0] + v2[0];
    div[nx-1]      = -v1[nx - 2] + v2[nx - 1];
    div[(ny-1)*nx] =  v1[(ny - 1)*nx] - v2[(ny - 2)*nx];
    div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny - 1)*nx - 1];
}

void forward_gradient(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int nx,   //image width
        const int ny    //image height
        ){
    // compute the gradient on the central body of the image
//#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ny-1; i++){
        for(int j = 0; j < nx-1; j++){
            const int p  = i * nx + j;
            const int p1 = p + 1;
            const int p2 = p + nx;

            fx[p] = f[p1] - f[p];
            fy[p] = f[p2] - f[p];
        }
    }

    // compute the gradient on the last row
    for (int j = 0; j < nx-1; j++){
        const int p = (ny-1) * nx + j;

        fx[p] = f[p+1] - f[p];
        fy[p] = 0;
    }

    // compute the gradient on the last column
    for (int i = 1; i < ny; i++){
        const int p = i * nx-1;

        fx[p] = 0;
        fy[p] = f[p+nx] - f[p];
    }

    fx[ny * nx - 1] = 0;
    fy[ny * nx - 1] = 0;
}


void backward_gradient(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int nx,   //image width
        const int ny    //image height
        ){
    // compute the gradient on the central body of the image
//#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < ny; i++){
        for(int j = 1; j < nx; j++){
            const int p  = i * nx + j;
            const int p1 = p - 1;
            const int p2 = p - nx;

            fx[p] = f[p] - f[p1];
            fy[p] = f[p] - f[p2];
        }
    }

    // compute the gradient on the first row
    for (int j = 1; j < nx; j++){
        const int p = j;

        fx[p] = f[p] - f[p-1];
        fy[p] = 0;
    }

    // compute the gradient on the first column
    for (int i = 1; i < ny; i++){
        const int p = i * nx;

        fx[p] = 0;
        fy[p] = f[p] - f[p-nx];
    }

    fx[0] = 0;
    fy[0] = 0;
}


void centered_gradient(
        const float *input,  //input image
        float *dx,           //computed x derivative
        float *dy,           //computed y derivative
        const int nx,        //image width
        const int ny         //image height
        ) {

    // compute the gradient on the center body of the image
//#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < ny-1; i++){
        for(int j = 1; j < nx-1; j++){

            const int k = i * nx + j;
            dx[k] = 0.5*(input[k + 1] - input[k - 1]);
            dy[k] = 0.5*(input[k + nx] - input[k - nx]);

        }
    }

    // compute the gradient on the first and last rows
    for (int j = 1; j < nx-1; j++){

        dx[j] = 0.5*(input[j + 1] - input[j - 1]);
        dy[j] = 0.5*(input[j + nx] - input[j]);

        const int k = (ny - 1) * nx + j;

        dx[k] = 0.5*(input[k + 1] - input[k - 1]);
        dy[k] = 0.5*(input[k] - input[k - nx]);
    }

    // compute the gradient on the first and last columns
    for(int i = 1; i < ny-1; i++){

        const int p = i * nx;
        dx[p] = 0.5*(input[p + 1] - input[p]);
        dy[p] = 0.5*(input[p + nx] - input[p - nx]);

        const int k = (i + 1) * nx - 1;
        dx[k] = 0.5*(input[k] - input[k - 1]);
        dy[k] = 0.5*(input[k + nx] - input[k - nx]);
    }

    // compute the gradient at the four corners
    dx[0] = 0.5*(input[1] - input[0]);
    dy[0] = 0.5*(input[nx] - input[0]);

    dx[nx - 1] = 0.5*(input[nx-1] - input[nx-2]);
    dy[nx - 1] = 0.5*(input[2*nx-1] - input[nx-1]);

    dx[(ny - 1)*nx] = 0.5*(input[(ny - 1)*nx + 1] - input[(ny - 1)*nx]);
    dy[(ny - 1)*nx] = 0.5*(input[(ny - 1)*nx] - input[(ny - 2)*nx]);

    dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
    dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}




/**
 *
 * Function to compute the gradient with five-point derivaties (1/12)*[-1 8 0 -8 1]
 */
void five_point_gradient(
        const float *input,  //input image
        float *dx,           //computed x derivative
        float *dy,           //computed y derivative
        const int nx,        //image width
        const int ny         //image height
        ){


    // compute the gradient on the center body of the image
//#pragma omp parallel for schedule(dynamic)
    for (int i = 2; i < ny-2; i++){
        for(int j = 2; j < nx-2; j++)
        {
            const int k = i * nx + j;
            dx[k] = (1.0/12)*(input[k-2] - 8*input[k-1] + 8*input[k+1] - input[k+2]);
            dy[k] = (1.0/12)*(input[k-2*nx] - 8*input[k-nx] + 8*input[k+ nx] - input[k+2*nx]);
        }
    }


    //Centered gradient for the second column an the penultimate column
    for (int j = 1; j < ny-1; j++){
        dx[j*nx + 1] = 0.5*(input[j*nx + 2] - input[j*nx -1]);
        dy[j*nx + 1] = 0.5*(input[j*nx + nx +1] - input[j*nx - nx +1]);

        const int k = j*nx + (ny-2);

        dx[k] = 0.5*(input[k+1] - input[k-1]);
        dy[k] = 0.5*(input[k + nx] - input[k-nx]);
    }

    //Centered gradient for the second row an the penultimate row
    for (int j = 1; j < ny-1; j++){
        dx[nx + j] = 0.5*(input[nx + j + 1] - input[nx + j -1]);
        dy[nx + j] = 0.5*(input[2*nx + j] - input[j]);

        const int k = nx*(ny-2);

        dx[k+j] = 0.5*(input[k+j+1] - input[k+j-1]);
        dy[k+j] = 0.5*(input[k + nx + j] - input[k - nx +j]);
    }

    // compute the gradient on the first and last rows
    for (int j = 1; j < nx-1; j++){
        dx[j] = 0.5*(input[j+1] - input[j-1]);
        dy[j] = 0.5*(input[j+nx] - input[j]);

        const int k = (ny - 1) * nx + j;

        dx[k] = 0.5*(input[k+1] - input[k-1]);
        dy[k] = 0.5*(input[k] - input[k-nx]);
    }

    // compute the gradient on the first and last columns
    for(int i = 1; i < ny-1; i++){
        const int p = i * nx;
        dx[p] = 0.5*(input[p+1] - input[p]);
        dy[p] = 0.5*(input[p+nx] - input[p-nx]);

        const int k = (i+1) * nx - 1;

        dx[k] = 0.5*(input[k] - input[k-1]);
        dy[k] = 0.5*(input[k+nx] - input[k-nx]);
    }

    // compute the gradient at the four corners
    dx[0] = 0.5*(input[1] - input[0]);
    dy[0] = 0.5*(input[nx] - input[0]);

    dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
    dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

    dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
    dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

    dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
    dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}


////////////////////////////////
/////////SMOOTH FUNCTIONS/////
//////////////////////////////
/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void gaussian(
        float *I,             // input/output image
        const int xdim,       // image width
        const int ydim,       // image height
        const float sigma    // Gaussian sigma
        ){
    const int boundary_condition = DEFAULT_BOUNDARY_CONDITION;
    const int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

    const float den  = 2*sigma*sigma;
    const int   size = (int) (window_size * sigma) + 1 ;
    const int   bdx  = xdim + size;
    const int   bdy  = ydim + size;

    if (boundary_condition && size > xdim) {
        fprintf(stderr, "GaussianSmooth: sigma too large\n");
        abort();
    }

    // compute the coefficients of the 1D convolution kernel
    float B[size];
    for(int i = 0; i < size; i++)
        B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

    // normalize the 1D convolution kernel
    float norm = 0;
    for(int i = 0; i < size; i++)
        norm += B[i];
    norm *= 2;
    norm -= B[0];
    for(int i = 0; i < size; i++)
        B[i] /= norm;

    // convolution of each line of the input image
    float *R = (float*) xmalloc((size + xdim + size)*sizeof*R);

    for (int k = 0; k < ydim; k++){
        int i, j;
        for (i = size; i < bdx; i++)
            R[i] = I[k * xdim + i - size];

        switch (boundary_condition){
        case BOUNDARY_CONDITION_DIRICHLET:
            for(i = 0, j = bdx; i < size; i++, j++)
                R[i] = R[j] = 0;
            break;

        case BOUNDARY_CONDITION_REFLECTING:
            for(i = 0, j = bdx; i < size; i++, j++) {
                R[i] = I[k * xdim + size-i];
                R[j] = I[k * xdim + xdim-i-1];
            }
            break;

        case BOUNDARY_CONDITION_PERIODIC:
            for(i = 0, j = bdx; i < size; i++, j++) {
                R[i] = I[k * xdim + xdim-size+i];
                R[j] = I[k * xdim + i];
            }
            break;
        }

        for (i = size; i < bdx; i++){
            float sum = B[0] * R[i];
            for (j = 1; j < size; j++ )
                sum += B[j] * ( R[i-j] + R[i+j] );
            I[k * xdim + i - size] = sum;
        }
    }

    // convolution of each column of the input image
    float *T = (float*) xmalloc((size + ydim + size)*sizeof*T);

    for (int k = 0; k < xdim; k++){
        int i, j;
        for (i = size; i < bdy; i++)
            T[i] = I[(i - size) * xdim + k];

        switch (boundary_condition){
        case BOUNDARY_CONDITION_DIRICHLET:
            for (i = 0, j = bdy; i < size; i++, j++)
                T[i] = T[j] = 0;
            break;

        case BOUNDARY_CONDITION_REFLECTING:
            for (i = 0, j = bdy; i < size; i++, j++) {
                T[i] = I[(size-i) * xdim + k];
                T[j] = I[(ydim-i-1) * xdim + k];
            }
            break;

        case BOUNDARY_CONDITION_PERIODIC:
            for( i = 0, j = bdx; i < size; i++, j++) {
                T[i] = I[(ydim-size+i) * xdim + k];
                T[j] = I[i * xdim + k];
            }
            break;
        }

        for (i = size; i < bdy; i++){
            float sum = B[0] * T[i];
            for (j = 1; j < size; j++ )
                sum += B[j] * (T[i-j] + T[i+j]);
            I[(i - size) * xdim + k] = sum;
        }
    }

    free(R);
    free(T);
}


void gaussian1Dweight(
        float *I,             // input/output image
        const int r       // image width
        ){
    const float sigma = r*0.3333;

    const float den  = 2*sigma*sigma;

    for (int i = 0; i < (2*r + 1); i++){
        I[i] = 1.0;
        I[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-(i-r) * (i-r) / den);
        // fprintf(stderr,"%d: %f \n",i, I[i]);
    }

    // // compute the coefficients of the 1D convolution kernel
    // for(int i = -r; i < (2*r + 1); i++)
    //   I[i + r] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);
    //   I[i+r] = 1.0;
}





//////////////////////////////////
////////MAX AND MIN FUNCTIONS////
/////////////////////////////////

float max(float a, float b){
    if (a > b){
        return a;
    }else{
        return b;
    }
}

float min(float a, float b){
    if (a < b){
        return a;
    }else{
        return b;
    }
}
/**
 * Compute the max and min of an array
 **/
void getminmax(
        float *min,     // output min
        float *max,     // output max
        const float *x, // input array
        int n           // array size
        ) {
    *min = *max = x[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < *min)
            *min = x[i];
        if (x[i] > *max)
            *max = x[i];
    }
}

/////////////////////////////////////
/////////IMAGE NORMALIZATION/////////
////////////////////////////////////

/**
 *
 * Function to normalize two images between 0 and 1
 *
 **/
void image_normalization(
        const float *I0,  // input image0
        const float *I1,  // input image1
        float *I0n,       // normalized output image0
        float *I1n,       // normalized output image1
        int size        // size of the image
        ) {
    float max0, max1, min0, min1;

    // obtain the max and min of each image
    getminmax(&min0, &max0, I0, size);
    getminmax(&min1, &max1, I1, size);

    // obtain the max and min of both images
    const float max = (max0 > max1)? max0 : max1;
    const float min = (min0 < min1)? min0 : min1;
    const float den = max - min;

    if (den > 0)
        // normalize both images between [0,1]
        for (int i = 0; i < size; i++){
            I0n[i] = (I0[i] - min) / den;
            I1n[i] = (I1[i] - min) / den;
        }

    else
        // copy the original images
        for (int i = 0; i < size; i++){
            I0n[i] = I0[i];
            I1n[i] = I1[i];
        }
}


/**
 *
 * Function to normalize three images between 0 and 1
 *
 **/

void image_normalization_3(
        const float *I1,  // input image1
        const float *I2,  // input image2
        const float *I0,  // input image0
        float *I1n,       // normalized output image1
        float *I2n,       // normalized output image2
        float *I0n,       // normalized output image0
        int size        // size of the image
        ) {
    float max0, max1, max2, min0, min1, min2;

    // Obtain the max and min of each image
    getminmax(&min0, &max0, I0, size);
    getminmax(&min1, &max1, I1, size);
    getminmax(&min2, &max2, I2, size);

    // Obtain the max and min of all images
    const float max01 = (max0 > max1)? max0 : max1;
    const float max = (max2 > max01)? max2 : max01;
    const float min01 = (min0 < min1)? min0 : min1;
    const float min = (min2 > min01)? min2 : min01;
    const float den = max - min;

    if (den > 0)
        // normalize both images between [0, 1]
        for (int i = 0; i < size; i++){
            I0n[i] = (I0[i] - min) / den;
            I1n[i] = (I1[i] - min) / den;
            I2n[i] = (I2[i] - min) / den;
        }

    else
        // copy the original images
        for (int i = 0; i < size; i++){
            I0n[i] = I0[i];
            I1n[i] = I1[i];
            I2n[i] = I2[i];
        }
}


/**
 *
 * Function to normalize three images between 0 and 1
 *
 **/

void image_normalization_4(
        float *I0, // input image-1
        float *I1,  // input image0
        float *I_1,  // input image1
        float *I2, //Smooth vefsion of I0
        float *I0n,		// normalized output image -1
        float *I1n,       // normalized output image0
        float *I_1n,       // normalized output image1
        float *I2n,   // normalized output image filtI0
        int size          // size of the image
        ) {
    float max_1, max0, max1, max2, min_1, min0, min1, min2;

    // obtain the max and min of each image
    getminmax(&min_1, &max_1, I_1, size);
    getminmax(&min0, &max0, I0, size);
    getminmax(&min1, &max1, I1, size);
    getminmax(&min2, &max2, I2, size);

    // obtain the max and min of images
    double max = (max_1 > max0) ? max_1 : max0;
    max = (max > max1) ? max : max1;
    max = (max > max2) ? max : max2;

    double min = (min_1 < min0) ? min_1 : min0;
    min = (min < min1) ? min : min1;
    min = (min < min2) ? min : min2;

    const double den = max - min;

    if (den > 0)
        // normalize both images
        for (int i = 0; i < size; i++) {
            I_1n[i] = (I_1[i] - min) / den;
            I0n[i] = (I0[i] - min) / den;
            I1n[i] = (I1[i] - min) / den;
            I2n[i] = (I2[i] - min) / den;
        }else{
        // copy the original images
        for (int i = 0; i < size; i++) {
            I_1n[i] = I_1[i];
            I0n[i] = I0[i];
            I1n[i] = I1[i];
            I2n[i] = I2[i];
        }
    }
}

void init_weight(
        float *g,
        float *Ix,
        float *Iy,
        float size){

    const float gamma = PAR_DEFAULT_GAMMA;

    for (int i = 0; i < size; i++) {
        const float grad = sqrt(Ix[i] * Ix[i] + Iy[i] * Iy[i]);
        g[i] = 1/(1 + gamma*grad);

    }

}
