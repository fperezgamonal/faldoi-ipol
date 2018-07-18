#ifndef TVCSAD_MODEL_W
#define TVCSAD_MODEL_W

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
#define TVCSADW_LAMBDA  40//40
#define TVCSADW_THETA   0.3
#define TVCSADW_TAU     0.125 //0.25
#define TVCSADW_NWARPS  1  //5
#define TVCSADW_TOL_D   0.01
#define TVCSADW_VERBOSE 0  //0



void  initialize_stuff_tvcsad_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore,
          const int w,
          const int h)

{
  // w, h as params in the function call
  //const int w = ofCore->params.w;
  //const int h = ofCore->params.h;
  //fprintf(stderr, "W x H :%d x %d\n", w, h);
  ofStuff->tvcsadw.weight = new float[ofCore->params.w_radio*2 + 1];
  ofStuff->tvcsadw.pnei = new PosNei[w*h];
  ofStuff->tvcsadw.xi11 = new float[w*h];
  ofStuff->tvcsadw.xi12 = new float[w*h];
  ofStuff->tvcsadw.xi21 = new float[w*h];
  ofStuff->tvcsadw.xi22 = new float[w*h];
  ofStuff->tvcsadw.u1x = new float[w*h];
  ofStuff->tvcsadw.u1y = new float[w*h];
  ofStuff->tvcsadw.u2x = new float[w*h];
  ofStuff->tvcsadw.u2y = new float[w*h];
  ofStuff->tvcsadw.v1 =  new float[w*h];
  ofStuff->tvcsadw.v2 =  new float[w*h];
  ofStuff->tvcsadw.rho_c =  new float[w*h];
  ofStuff->tvcsadw.grad =  new float[w*h];
  ofStuff->tvcsadw.u1_ =  new float[w*h];
  ofStuff->tvcsadw.u2_ =  new float[w*h];
  ofStuff->tvcsadw.u1_tmp = new float[w*h];
  ofStuff->tvcsadw.u2_tmp = new float[w*h];  
  ofStuff->tvcsadw.I1x = new float[w*h];
  ofStuff->tvcsadw.I1y = new float[w*h]; 
  ofStuff->tvcsadw.I1w = new float[w*h]; 
  ofStuff->tvcsadw.I1wx = new float[w*h]; 
  ofStuff->tvcsadw.I1wy = new float[w*h]; 
  ofStuff->tvcsadw.div_xi1 = new float[w*h]; 
  ofStuff->tvcsadw.div_xi2 = new float[w*h]; 
}


void  free_stuff_tvcsad_w(SpecificOFStuff *ofStuff){
  delete [] ofStuff->tvcsadw.weight;
  delete [] ofStuff->tvcsadw.pnei;
  delete [] ofStuff->tvcsadw.xi11;
  delete [] ofStuff->tvcsadw.xi12;
  delete [] ofStuff->tvcsadw.xi21;
  delete [] ofStuff->tvcsadw.xi22;
  delete [] ofStuff->tvcsadw.u1x;
  delete [] ofStuff->tvcsadw.u1y;
  delete [] ofStuff->tvcsadw.u2x;
  delete [] ofStuff->tvcsadw.u2y;
  delete [] ofStuff->tvcsadw.v1;
  delete [] ofStuff->tvcsadw.v2;
  delete [] ofStuff->tvcsadw.rho_c;
  delete [] ofStuff->tvcsadw.grad;
  delete [] ofStuff->tvcsadw.u1_;
  delete [] ofStuff->tvcsadw.u2_;
  delete [] ofStuff->tvcsadw.u1_tmp;
  delete [] ofStuff->tvcsadw.u2_tmp;  
  delete [] ofStuff->tvcsadw.I1x;
  delete [] ofStuff->tvcsadw.I1y; 
  delete [] ofStuff->tvcsadw.I1w; 
  delete [] ofStuff->tvcsadw.I1wx; 
  delete [] ofStuff->tvcsadw.I1wy; 
  delete [] ofStuff->tvcsadw.div_xi1; 
  delete [] ofStuff->tvcsadw.div_xi2; 
}

void eval_tvcsad_w(
    const float *I0,
    const float *I1,
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
    float *ener_N,
    const int ii,
    const int ij,
    const int ei,
    const int ej,
    const float lambda,
    const float theta,
    const int nx,
    const int ny
    )
{
  float *u1 = ofD->u1;
  float *u2 = ofD->u2;

  // w, h as params in the function call

  //Columns and Rows
  //const int nx = ofD->params.w;
  //const int ny = ofD->params.h;

  //Optical flow derivatives
  float *v1   = tvcsadw->v1;
  float *v2   = tvcsadw->v2;
  float *u1x  = tvcsadw->u1x;
  float *u1y  = tvcsadw->u1y;
  float *u2x  = tvcsadw->u2x;
  float *u2y  = tvcsadw->u2y;


  float *I1w = tvcsadw->I1w;

  PosNei *pnei  = tvcsadw->pnei;

  int ndt = DT_NEI;

  const int iiw = tvcsadw->iiw;
  const int ijw = tvcsadw->ijw;
  float *weight = tvcsadw->weight;


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
    dt *=lambda*weight[l-ij + ijw]*weight[k-ii + iiw];

    assert(g>=0);
    assert(dt>=0);
    ener +=dc + dt + g;
    m++;
    if (!std::isfinite(dt))
        std::printf("Corrupt data\n");
    if (!std::isfinite(g))
        std::printf("Corrupt regularization\n");

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

void tvcsad_w_getP(
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
void tvcsad_w_getD(
      float *xi11, 
      float *xi12, 
      float *xi21, 
      float *xi22, 
      const float *u1x,
      const float *u1y,
      const float *u2x,
      const float *u2y,
      float tau,       
      const int ii, // initial column
      const int ij, // initial row
      const int ei, // end column
      const int ej, // end row
      const int nx
){
//#pragma omp parallel for
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
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


void guided_tvcsad_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    TvCsadStuff_W *tvcsadw,
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
    const bool  verbose, // enable/disable the verbose mode
    const int nx,        // width of I0 (and I1)
    const int ny         // height of I0 (and I1)
  )
{

  const float l_t = lambda * theta;

  const int ndt = DT_NEI;


  // Added changes for subimages

  float *u1 = ofD->u1;
  float *u2 = ofD->u2;
  //Columns and Rows
  //const int nx = ofD->params.w;
  //const int ny = ofD->params.h;

  PosNei *pnei  = tvcsadw->pnei;
  float *u1_  = tvcsadw->u1_;
  float *u2_  = tvcsadw->u2_;
  //Optical flow derivatives
  float *u1x  = tvcsadw->u1x;
  float *u1y  = tvcsadw->u1y;
  float *u2x  = tvcsadw->u2x;
  float *u2y  = tvcsadw->u2y;
  //Dual variables
  float *xi11 = tvcsadw->xi11;
  float *xi12 = tvcsadw->xi12;
  float *xi21 = tvcsadw->xi21;
  float *xi22 = tvcsadw->xi22;
  
  float *v1 = tvcsadw->v1;
  float *v2 = tvcsadw->v2;

  float *grad  = tvcsadw->grad;

  float *u1_tmp = tvcsadw->u1_tmp;
  float *u2_tmp = tvcsadw->u2_tmp;

  float *I1x = tvcsadw->I1x;
  float *I1y = tvcsadw->I1y;

  float *I1w = tvcsadw->I1w;
  float *I1wx = tvcsadw->I1wx;
  float *I1wy = tvcsadw->I1wy;

  //Divergence
  float *div_xi1 = tvcsadw->div_xi1;
  float *div_xi2 = tvcsadw->div_xi2;

  //TODO:Weights
  const int iiw = tvcsadw->iiw;
  const int ijw = tvcsadw->ijw;
  float *weight = tvcsadw->weight;

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
// #pragma omp parallel for
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
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
        if (positive(ap))
        {
          // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, pnei[i].apj[j]*nx + pnei[i].api[j], pnei[i].apj[j], pnei[i].api[j]);
          assert(pnei[i].api[j] >= 0);
          assert(pnei[i].apj[j] >= 0);
          assert(pnei[i].apj[j]*nx + pnei[i].api[j] < nx*ny);
          const int pos = apj*nx + api;
          
          pnei[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                 + I1wy[i] * u2[i])/grad[i];
          n_tmp ++;
        }
      }
      pnei[i].n= n_tmp;
    }
    }

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
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
      // (thresholding opterator TH)
// #pragma omp parallel for
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const float l_t_w = l_t * weight[l-ij + ijw]*weight[k-ii + iiw];
        const int i = l*nx + k;
        int it = 0;
        for (int j = 0; j< ndt; j++)
        {
          const int api = pnei[i].api[j];
          const int apj = pnei[i].apj[j];
          const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);

          if (positive(ap))
          {
            // std::printf("J:%d I:%d\n",pnei[i].api[j],pnei[i].apj[j]);
            // const int pos = apj*nx + api;

            pnei[i].ba[it] = -(pnei[i].b[j] -  (I1wx[i] * u1[i]
                   + I1wy[i] * u2[i])/grad[i]);
            it++;
          }
        }
        for (int j = 0; j < (pnei[i].n+1); j++)
        {
           pnei[i].ba[it]= (pnei[i].n - 2*j)*l_t_w*grad[i];
           it++;
        }
  
        std::sort(pnei[i].ba.begin(), pnei[i].ba.begin() + it);
        // v1[i] = u1[i] - l_t*I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
        // v2[i] = u2[i] - l_t*I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
        //TODO: Possible error in the minimization
        v1[i] = u1[i] - I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
        v2[i] = u2[i] - I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
      }
      }

      //Dual variables
      tvcsad_w_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y, 
          tau, ii, ij, ei, ej, nx);  

      //Primal variables
            //Primal variables
      divergence_patch(xi11,xi12,div_xi1,ii,ij,ei,ej,nx);
      divergence_patch(xi21,xi22,div_xi2,ii,ij,ei,ej,nx);

      //Almacenamos la iteracion anterior
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*nx + k;
        u1_tmp[i] = u1[i];
        u2_tmp[i] = u2[i];   
      }
      }
      tvcsad_w_getP(u1, u2, v1, v2, div_xi1, div_xi2, theta, tau,
          ii, ij, ei, ej, nx, &err_D);

      //(aceleration = 1);
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
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
  eval_tvcsad_w(I0, I1, ofD, tvcsadw, ener_N, ii, ij, ei, ej, lambda, theta, nx, ny);

}

#endif //TVCSAD 
