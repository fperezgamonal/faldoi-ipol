#ifndef UTILS_H
#define UTILS_H

// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include <chrono>
#include <iostream>
#include <utility>
#include "energy_structures.h"


PatchIndexes get_index_patch(
        const int wr,
        const int w,
        const int h,
        const int i,
        const int j,
        const int factor
        );


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
        );


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
        );


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
        );


/////////////////////////////////////////
///////////DERIVATIVES IMAGE////////////
////////////////////////////////////////


void divergence(
        const float *v1, // x component of the vector field
        const float *v2, // y component of the vector field
        float *div,      // output divergence
        const int nx,    // image width
        const int ny     // image height
        );
void forward_gradient(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int nx,   //image width
        const int ny    //image height
        );


void backward_gradient(
        const float *f, //input image
        float *fx,      //computed x derivative
        float *fy,      //computed y derivative
        const int nx,   //image width
        const int ny    //image height
        );


void centered_gradient(
        const float *input,  //input image
        float *dx,           //computed x derivative
        float *dy,           //computed y derivative
        const int nx,        //image width
        const int ny         //image height
        );

void five_point_gradient(
        const float *input,  //input image
        float *dx,           //computed x derivative
        float *dy,           //computed y derivative
        const int nx,        //image width
        const int ny         //image height
        );


////////////////////////////////
/////////SMOOTH FUNCTIONS/////
//////////////////////////////

void gaussian(
        float *I,             // input/output image
        const int xdim,       // image width
        const int ydim,       // image height
        const float sigma    // Gaussian sigma
        );


void gaussian1Dweight(
        float *I,             // input/output image
        const int r       // image width
        );





//////////////////////////////////
////////MAX AND MIN FUNCTIONS////
/////////////////////////////////

float max(float a, float b);
float min(float a, float b);
void getminmax(
        float *min,     // output min
        float *max,     // output max
        const float *x, // input array
        int n           // array size
        );

/////////////////////////////////////
/////////IMAGE NORMALIZATION/////////
////////////////////////////////////

void image_normalization(
        const float *I0,  // input image0
        const float *I1,  // input image1
        float *I0n,       // normalized output image0
        float *I1n,       // normalized output image1
        int size        // size of the image
        );

void image_normalization_3(
        const float *I1,  // input image1
        const float *I2,  // input image2
        const float *I0,  // input image0
        float *I1n,       // normalized output image1
        float *I2n,       // normalized output image2
        float *I0n,       // normalized output image0
        int size        // size of the image
        );

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
        );

void init_weight(
        float *g,
        float *Ix,
        float *Iy,
        float size);

template <typename F>
void time_it(const char* name, F&& f)
{
    using namespace std::chrono;

    auto start = steady_clock::now();

    std::forward<F>(f)();

    auto end = steady_clock::now();

    std::cout << "Time elapsed for " << name << ": " << duration_cast<milliseconds>(end - start).count() << " ms\n";
}

#endif // UTILS_H
