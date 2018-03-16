#ifndef MASK_H
#define MASK_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(
		const float *v1, // x component of the vector field
		const float *v2, // y component of the vector field
		float *div,      // output divergence
		const int nx,    // image width
		const int ny     // image height
	       );


/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(
		const float *f, //input image
		float *fx,      //computed x derivative
		float *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
		);
/**
 *
 * Function to compute the gradient with backward differences
 * (see [2] for details)
 *
 **/
void backward_gradient(
		const float *f, //input image
		float *fx,      //computed x derivative
		float *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
		);

/**
 *
 * Function to compute the gradient with centered differences
 *
 **/
void centered_gradient(
		const float *input,  //input image
		float *dx,           //computed x derivative
		float *dy,           //computed y derivative
		const int nx,        //image width
		const int ny         //image height
		);


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
);

void gaussian1Dweight(
  float *I,             // input/output image
  const int r       // image width
);

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
		);


#endif//MASK_C
