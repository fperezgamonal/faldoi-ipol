#ifndef TVCSAD_MODEL
#define TVCSAD_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm> 

#include "energy_structures.h"
#include "aux_energy_model.h"
#include "utils.h"
extern "C" {
#include "bicubic_interpolation.h"
}
//OPTICAL FLOW PARAMETERS
#define TVCSAD_LAMBDA  40//40
#define TVCSAD_THETA   0.3
#define TVCSAD_TAU     0.125 //0.25
#define TVCSAD_NWARPS  1  //5
#define TVCSAD_TOL_D   0.01
#define TVCSAD_VERBOSE 0  //0



void  intialize_stuff_tvcsad(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore)

{
  const int w = ofCore->params.w;
  const int h = ofCore->params.h;
  //fprintf(stderr, "W x H :%d x %d\n", w, h);
  ofStuff->tvcsad.pnei = new PosNei[w*h];
  ofStuff->tvcsad.xi11 = new float[w*h];
  ofStuff->tvcsad.xi12 = new float[w*h];
  ofStuff->tvcsad.xi21 = new float[w*h];
  ofStuff->tvcsad.xi22 = new float[w*h];
  ofStuff->tvcsad.u1x = new float[w*h];
  ofStuff->tvcsad.u1y = new float[w*h];
  ofStuff->tvcsad.u2x = new float[w*h];
  ofStuff->tvcsad.u2y = new float[w*h];
  ofStuff->tvcsad.v1 =  new float[w*h];
  ofStuff->tvcsad.v2 =  new float[w*h];
  ofStuff->tvcsad.rho_c =  new float[w*h];
  ofStuff->tvcsad.grad =  new float[w*h];
  ofStuff->tvcsad.u1_ =  new float[w*h];
  ofStuff->tvcsad.u2_ =  new float[w*h];
  ofStuff->tvcsad.u1_tmp = new float[w*h];
  ofStuff->tvcsad.u2_tmp = new float[w*h];  
  ofStuff->tvcsad.I1x = new float[w*h];
  ofStuff->tvcsad.I1y = new float[w*h]; 
  ofStuff->tvcsad.I1w = new float[w*h]; 
  ofStuff->tvcsad.I1wx = new float[w*h]; 
  ofStuff->tvcsad.I1wy = new float[w*h]; 
  ofStuff->tvcsad.div_xi1 = new float[w*h]; 
  ofStuff->tvcsad.div_xi2 = new float[w*h]; 
}


void  free_stuff_tvcsad(SpecificOFStuff *ofStuff){
  delete [] ofStuff->tvcsad.pnei;
  delete [] ofStuff->tvcsad.xi11;
  delete [] ofStuff->tvcsad.xi12;
  delete [] ofStuff->tvcsad.xi21;
  delete [] ofStuff->tvcsad.xi22;
  delete [] ofStuff->tvcsad.u1x;
  delete [] ofStuff->tvcsad.u1y;
  delete [] ofStuff->tvcsad.u2x;
  delete [] ofStuff->tvcsad.u2y;
  delete [] ofStuff->tvcsad.v1;
  delete [] ofStuff->tvcsad.v2;
  delete [] ofStuff->tvcsad.rho_c;
  delete [] ofStuff->tvcsad.grad;
  delete [] ofStuff->tvcsad.u1_;
  delete [] ofStuff->tvcsad.u2_;
  delete [] ofStuff->tvcsad.u1_tmp;
  delete [] ofStuff->tvcsad.u2_tmp;  
  delete [] ofStuff->tvcsad.I1x;
  delete [] ofStuff->tvcsad.I1y; 
  delete [] ofStuff->tvcsad.I1w; 
  delete [] ofStuff->tvcsad.I1wx; 
  delete [] ofStuff->tvcsad.I1wy; 
  delete [] ofStuff->tvcsad.div_xi1; 
  delete [] ofStuff->tvcsad.div_xi2; 
}

void eval_tvcsad(
    const float *I0,
    const float *I1,
    OpticalFlowData *ofD,
    TvCsadStuff *tvcsad,
    float *ener_N,
    const int ii,
    const int ij,
    const int ei,
    const int ej,
    const float lambda,
    const float theta
    )
{
  float *u1 = ofD->u1;
  float *u2 = ofD->u2;


  //Columns and Rows
  const int nx = ofD->params.w;
  const int ny = ofD->params.h;

  //Optical flow derivatives
  float *v1   = tvcsad->v1;
  float *v2   = tvcsad->v2;
  float *u1x  = tvcsad->u1x;
  float *u1y  = tvcsad->u1y;
  float *u2x  = tvcsad->u2x;
  float *u2y  = tvcsad->u2y;


  float *I1w = tvcsad->I1w;

  PosNei *pnei  = tvcsad->pnei;

  int ndt = DT_NEI;

  bicubic_interpolation_warp_patch(I1,  u1, u2, I1w, 
                              ii, ij, ei, ej, nx, ny, false);
  float ener = 0.0;
  int m = 0;
  for (int l = ij; l < ej; l++){
  for (int k = ii; k < ei; k++){
    const int i = l*nx + k;
    float dc = (1/(2*theta))*
        ((u1[i]-v1[i])*(u1[i]-v1[i]) + (u2[i] - v2[i])*(u2[i] - v2[i]));
    float g1  = u1x[i]*u1x[i];
    float g12 = u1y[i]*u1y[i];
    float g21 = u2x[i]*u2x[i];
    float g2  = u2y[i]*u2y[i];
    float g  = sqrt(g1 + g12 + g21 + g2);

    float dt = 0.0;
    for (int j = 0; j < ndt; j++)
    {
      const int api = pnei[i].api[j];
      const int apj = pnei[i].apj[j];
      const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);
      if (positive(ap))
      {
        // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, pnei[i].apj[j]*nx + pnei[i].api[j], pnei[i].apj[j], pnei[i].api[j]);
        assert(pnei[i].api[j] >= 0);
        assert(pnei[i].apj[j] >= 0);
        assert(pnei[i].apj[j]*nx + pnei[i].api[j] < nx*ny);
        const int pos = apj*nx + api;
        dt +=  fabs(I0[i] - I0[pos] - I1w[i] + I1w[pos]);
      }
    }
    dt *=lambda;

    assert(g>=0);
    assert(dt>=0);
    ener +=dc + dt + g;
    m++;
    if (!std::isfinite(dt))
        std::printf("Datos corruptos\n");
    if (!std::isfinite(g))
        std::printf("Regularizacion corrupta\n");

  }
  }
  assert(ener>=0);
  ener /=(m*1.0);
  (*ener_N) = ener;
}

//Chambolle functions
/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/

void tvcsad_getP(
      float *u1,
      float *u2,
      float *v1,
      float *v2,
      float *div_xi1,
      float *div_xi2, 
      float theta,
      float tau,
      const int ii, // initial column
      const int ij, // initial row
      const int ei, // end column
      const int ej, // end row
      const int nx,
      float *err)
{
  float err_D = 0.0;

//#pragma omp parallel for reduction(+:err_D)
  for (int l = ij; l < ej; l++){
  for (int k = ii; k < ei; k++){
    const int i = l*nx + k;

    const float u1k = u1[i];
    const float u2k = u2[i];
    // if (mask[i]==0){
      u1[i] = u1k  -tau*(-div_xi1[i]  + (u1k - v1[i])/theta);
      u2[i] = u2k  -tau*(-div_xi2[i]  + (u2k - v2[i])/theta);

      err_D += (u1[i] - u1k) * (u1[i] - u1k) +
            (u2[i] - u2k) * (u2[i] - u2k);
    // }
  }
  }
    
  err_D /= (ej-ij)*(ei-ii);
  (*err) = err_D;
}


/*
 * - Name: getD_Du

 *
*/
void tvcsad_getD(
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
//#pragma omp parallel for
#pragma omp parallel for schedule(dynamic,1) collapse(2)
  for (int l = ij; l < ej; l++){
  for (int k = ii; k < ei; k++){
    const int i = l*nx + k;
    float xi1_N = hypot(xi11[i],xi12[i]);
    float xi2_N = hypot(xi21[i],xi22[i]);
    
    xi1_N = MAX(1,xi1_N);
    xi2_N = MAX(1,xi2_N);
    
    xi11[i] = (xi11[i] + tau*u1x[i])/xi1_N;
    xi12[i] = (xi12[i] + tau*u1y[i])/xi1_N;
    
    xi21[i] = (xi21[i] + tau*u2x[i])/xi2_N;
    xi22[i] = (xi22[i] + tau*u2y[i])/xi2_N;
  }
  }
}


void guided_tvcsad(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    TvCsadStuff *tvcsad,
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
  )
{

  const float l_t = lambda * theta;

  const int ndt = DT_NEI;

  float *u1 = ofD->u1;
  float *u2 = ofD->u2;

  //Columns and Rows
  const int nx = ofD->params.w;
  const int ny = ofD->params.h;

  PosNei *pnei  = tvcsad->pnei;
  float *u1_  = tvcsad->u1_;
  float *u2_  = tvcsad->u2_;
  //Optical flow derivatives
  float *u1x  = tvcsad->u1x;
  float *u1y  = tvcsad->u1y;
  float *u2x  = tvcsad->u2x;
  float *u2y  = tvcsad->u2y;
  //Dual variables
  float *xi11 = tvcsad->xi11;
  float *xi12 = tvcsad->xi12;
  float *xi21 = tvcsad->xi21;
  float *xi22 = tvcsad->xi22;
  
  //Decouple variables
  float *v1 = tvcsad->v1;
  float *v2 = tvcsad->v2;

  float *grad  = tvcsad->grad;

  float *u1_tmp = tvcsad->u1_tmp;
  float *u2_tmp = tvcsad->u2_tmp;

  float *I1x = tvcsad->I1x;
  float *I1y = tvcsad->I1y;

  float *I1w = tvcsad->I1w;
  float *I1wx = tvcsad->I1wx;
  float *I1wy = tvcsad->I1wy;

  //Divergence
  float *div_xi1 = tvcsad->div_xi1;
  float *div_xi2 = tvcsad->div_xi2;

#pragma omp parallel for schedule(dynamic,1) collapse(2)
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
// #pragma omp parallel for
#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++){
    for (int k = ii; k < ei; k++){

      const int i = l*nx + k;
      const float Ix2 = I1wx[i] * I1wx[i];
      const float Iy2 = I1wy[i] * I1wy[i];

      // store the |Grad(I1(p + u))| (Warping image)
      grad[i] = hypot(Ix2 + Iy2,0.01);
      int n_tmp = 0;
      for (int j = 0; j < ndt; j++)
      {
        const int api = pnei[i].api[j];
        const int apj = pnei[i].apj[j];
        const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);

        // std::printf("I:%d Iter:%d J:%d I:%d \n",i, j,  pnei[i].apj[j], pnei[i].api[j]);
        if (ap == 0)
        {
          // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, pnei[i].apj[j]*nx + pnei[i].api[j], pnei[i].apj[j], pnei[i].api[j]);
          assert(api >= 0);
          assert(apj >= 0);
          assert(apj*nx + api < nx*ny);
          const int pos = apj*nx + api;
          
          pnei[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                 + I1wy[i] * u2[i])/grad[i];
          n_tmp ++;
        }
      }
      pnei[i].n= n_tmp;
    }
    }

#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++){
    for (int k = ii; k < ei; k++){
      const int i = l*nx + k;
      u1_[i] = u1[i];
      u2_[i] = u2[i];
    }
    }

    int n = 0;
    float err_D = INFINITY;
    while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS_LOCAL)
    {

      n++;
      // estimate the values of the variable (v1, v2)
      // (thresholding opterator TH)
// #pragma omp parallel for
#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*nx + k;
        int it = 0;

        for (int j = 0; j< ndt; j++)
        {
          const int api = pnei[i].api[j];
          const int apj = pnei[i].apj[j];
          const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);

          if (ap == 0)
          {
            pnei[i].ba[it] = -(pnei[i].b[j] -  (I1wx[i] * u1[i]
                   + I1wy[i] * u2[i])/grad[i]);
            it++;
          }
        }
        for (int j = 0; j < (pnei[i].n+1); j++)
        {
           pnei[i].ba[it]= (pnei[i].n - 2*j)*l_t*grad[i];
           it++;
        }
  
        std::sort(pnei[i].ba.begin(), pnei[i].ba.begin() + it);
        // v1[i] = u1[i] - l_t*I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
        // v2[i] = u2[i] - l_t*I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
        //TODO: Posible error en la minimizacion
        v1[i] = u1[i] - I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
        v2[i] = u2[i] - I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
      }
      }

      //Dual variables
      tvcsad_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y, 
          tau, ii, ij, ei, ej, nx);  

      //Primal variables
            //Primal variables
      divergence_patch(xi11,xi12,div_xi1,ii,ij,ei,ej,nx);
      divergence_patch(xi21,xi22,div_xi2,ii,ij,ei,ej,nx);

      //Almacenamos la iteracion anterior
#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*nx + k;
        u1_tmp[i] = u1[i];
        u2_tmp[i] = u2[i];   
      }
      }
      tvcsad_getP(u1, u2, v1, v2, div_xi1, div_xi2, theta, tau,
          ii, ij, ei, ej, nx, &err_D);

      //(aceleration = 1);
#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*nx + k;
        u1_[i] = 2*u1[i] - u1_tmp[i];
        u2_[i] = 2*u2[i] - u2_tmp[i];
      }
      }
    }
    if (verbose)
      fprintf(stderr, "Warping: %d,Iter: %d "
      "Error: %f\n", warpings,n, err_D);
  }
  eval_tvcsad(I0, I1, ofD, tvcsad, ener_N, ii, ij, ei, ej, lambda, theta);

}

#endif //TVCSAD 
