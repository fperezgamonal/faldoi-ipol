#ifndef NLTVL1_MODEL_W
#define NLTVL1_MODEL_W

#include <cmath>
#include <cstdio>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"
extern "C" {
#include "bicubic_interpolation.h"
}
void  initialize_stuff_nltvl1_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore,
          const int w,
          const int h)

{
  // w, h as params in the function call
  //const int w = ofCore->params.w;
  //const int h = ofCore->params.h;
  ofStuff->nltvl1w.weight = new float[ofCore->params.w_radio*2 + 1];
  ofStuff->nltvl1w.p = new DualVariables[w*h];
  ofStuff->nltvl1w.q = new DualVariables[w*h];
  ofStuff->nltvl1w.v1 =  new float[w*h];
  ofStuff->nltvl1w.v2 =  new float[w*h];
  ofStuff->nltvl1w.rho_c =  new float[w*h];
  ofStuff->nltvl1w.grad =  new float[w*h];
  ofStuff->nltvl1w.u1_ =  new float[w*h];
  ofStuff->nltvl1w.u2_ =  new float[w*h];
  ofStuff->nltvl1w.u1_tmp = new float[w*h];
  ofStuff->nltvl1w.u2_tmp = new float[w*h];  
  ofStuff->nltvl1w.I1x = new float[w*h];
  ofStuff->nltvl1w.I1y = new float[w*h]; 
  ofStuff->nltvl1w.I1w = new float[w*h]; 
  ofStuff->nltvl1w.I1wx = new float[w*h]; 
  ofStuff->nltvl1w.I1wy = new float[w*h]; 
  ofStuff->nltvl1w.div_p = new float[w*h]; 
  ofStuff->nltvl1w.div_q = new float[w*h]; 
}



void  free_stuff_nltvl1_w(SpecificOFStuff *ofStuff)

{
  delete [] ofStuff->nltvl1w.weight;
  delete [] ofStuff->nltvl1w.p;
  delete [] ofStuff->nltvl1w.q;
  delete [] ofStuff->nltvl1w.v1;
  delete [] ofStuff->nltvl1w.v2;
  delete [] ofStuff->nltvl1w.rho_c;
  delete [] ofStuff->nltvl1w.grad;
  delete [] ofStuff->nltvl1w.u1_;
  delete [] ofStuff->nltvl1w.u2_;
  delete [] ofStuff->nltvl1w.u1_tmp;
  delete [] ofStuff->nltvl1w.u2_tmp;  
  delete [] ofStuff->nltvl1w.I1x;
  delete [] ofStuff->nltvl1w.I1y; 
  delete [] ofStuff->nltvl1w.I1w; 
  delete [] ofStuff->nltvl1w.I1wx; 
  delete [] ofStuff->nltvl1w.I1wy; 
  delete [] ofStuff->nltvl1w.div_p; 
  delete [] ofStuff->nltvl1w.div_q; 
}


void eval_nltvl1_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta,
    const int w,
    const int h
){

  float *u1 = ofD->u1;
  float *u2 = ofD->u2;

  // Added changes for subimages

  //Columns and Rows
  //const int w = ofD->params.w;
  //const int h = ofD->params.h;

  float *I1w = nltvl1w->I1w;
  DualVariables *p = nltvl1w->p;
  DualVariables *q = nltvl1w->q;


  float *v1 = nltvl1w->v1;
  float *v2 = nltvl1w->v2;

  float ener = 0.0;
  int n_d = NL_DUAL_VAR;

    //TODO:Pesos
  const int iiw = nltvl1w->iiw;
  const int ijw = nltvl1w->ijw;
  float *weight = nltvl1w->weight;



  bicubic_interpolation_warp_patch(I1,  u1, u2, I1w, 
                              ii, ij, ei, ej, w, h, false);
  //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
  int m  = 0;
  for (int l = ij; l < ej; l++){
  for (int k = ii; k < ei; k++){
    const int i = l*w + k;
    float dt = lambda*fabs(I1w[i]-I0[i])*weight[l-ij + ijw]*weight[k-ii + iiw];
    float dc = (1/(2*theta))*
        ((u1[i]-v1[i])*(u1[i]-v1[i]) + (u2[i] - v2[i])*(u2[i] - v2[i]));
    float g = 0.0;
    for (int j = 0; j < n_d; j++)
    { 
      const int ap1i = p[i].api[j];
      const int ap1j = p[i].apj[j];
      const int ap2i = q[i].api[j];
      const int ap2j = q[i].apj[j];
      
      assert(ap1i == ap2i && ap1j == ap2j);
      //The position should be the same
      const int ap1 = validate_ap_patch(ii, ij, ei, ej, ap1i, ap1j);
      const int ap2 = validate_ap_patch(ii, ij, ei, ej, ap2i, ap2j);
      assert(ap1 == ap2);
      if (positive(ap1) && positive(ap2))
      {
        const float wp1 = p[i].wp[j]; 
        const float wp2 = q[i].wp[j]; 
        assert(wp1 == wp2);
        assert(wp1 >= 0);
        assert(wp2 >= 0);
        g += fabs(u1[i] - u1[ap1j*w + ap1i])*wp1 
            + fabs(u2[i] - u2[ap2j*w + ap2i])*wp2;
      }
    }
    assert(g>=0);
    g /=p[i].wt;
    assert(g>=0);
    assert(dt>=0);
    assert(p[i].wt == q[i].wt);
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




/*
 * - Name: getP

 *
*/

void nltvl1_w_getP(
            float *v1, 
            float *v2,
            float *div_p1, 
            float *div_p2,  
            float theta, 
            float tau, 
            const int ii, // initial column
            const int ij, // initial row
            const int ei, // end column
            const int ej, // end row
            const int w,
            float *u1,
            float *u2,
            float *err
    )
{
  float err_D = 0.0;

//#pragma omp parallel for reduction(+:err_D)
  for (int l = ij; l < ej; l++)
  for (int k = ii; k < ei; k++)
  {
    const int i = l*w + k;      
    //Only modify the inpainting domain
    // if (mask[i]==0){
      const float u1k = u1[i];
      const float u2k = u2[i];
        
      u1[i] = u1k  -tau*(div_p1[i]  + (u1k - v1[i])/theta);
      u2[i] = u2k  -tau*(div_p2[i]  + (u2k - v2[i])/theta);

      err_D += (u1[i] - u1k) * (u1[i] - u1k) +
            (u2[i] - u2k) * (u2[i] - u2k);
      assert(std::isfinite(v1[i]));
      assert(std::isfinite(v2[i]));
      assert(std::isfinite(u1[i]));
      assert(std::isfinite(u2[i]));
    // }
  }
  err_D /= (ej-ij)*(ei-ii);
  (*err) = err_D;
}


 void nltvl1_w_getD(
          float *u1,
          float *u2,
          const int ii, // initial column
          const int ij, // initial row
          const int ei, // end column
          const int ej, // end row
          const int w,
          int n_d,
          float tau,
          DualVariables *p1,
          DualVariables *p2
)
{

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
  for (int l = ij; l < ej; l++)
  for (int k = ii; k < ei; k++)
  {
    const int i = l*w + k;
    const float wt1 = p1[i].wt;
    const float wt2 = p2[i].wt;
    assert (wt1 == wt2 && wt1 >0);
    for (int j = 0; j < n_d; j++)
    {
      const int ap1i   = p1[i].api[j];
      const int ap1j   = p1[i].apj[j];
      const int ap2i   = p2[i].api[j];
      const int ap2j   = p2[i].apj[j];
      const float wp1  = p1[i].wp[j];
      const float wp2  = p2[i].wp[j];
      // std::printf("Ap1i:%d Ap1j:%d Ap2i:%d Ap2j:%d\n",ap1i, ap1j, ap2i, ap2j);
      assert(ap1i == ap2i && ap1j == ap2j);
      assert(wp1 == wp2);
      //The position should be the same
      const int ap1 = validate_ap_patch(ii, ij, ei, ej, ap1i, ap1j);
      const int ap2 = validate_ap_patch(ii, ij, ei, ej, ap2i, ap2j);
      assert(ap1 == ap2);
      if (positive(ap1) && positive(ap2))
      {
        assert(wp1 >=0);
        assert(wp2 >=0);
        
        const float u1x = u1[i];
        const float u2x = u2[i];
        const float u1y = u1[ap1j*w + ap1i];
        const float u2y = u2[ap2j*w + ap1i];

        const float nlgr1 =  wp1 * (u1x -u1y)/wt1;
        const float nlgr2 =  wp2 * (u2x -u2y)/wt2;
        const float nl1 = sqrt(nlgr1*nlgr1);
        const float nl2 = sqrt(nlgr2*nlgr2);
        const float nl1g = 1 + tau * nl1;
        const float nl2g = 1 + tau * nl2;

        p1[i].sc[j] =  (p1[i].sc[j] + tau *nlgr1)/nl1g;
        p2[i].sc[j] =  (p2[i].sc[j] + tau *nlgr2)/nl2g;
        assert(std::isfinite(p1[j].sc[j]));
        assert(std::isfinite(p2[i].sc[j]));
      }
    }
  }
}

void guided_nltvl1_w(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
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
    const bool  verbose,  // enable/disable the verbose mode
    const int w,
    const int h
) {
  float *u1 = ofD->u1;
  float *u2 = ofD->u2;
  // w, h as params in the function call

  //Columns and Rows
  //const int w = ofD->params.w;
  //const int h = ofD->params.h;

  DualVariables *p = nltvl1w->p;
  DualVariables *q = nltvl1w->q;

  float *u1_  = nltvl1w->u1_;
  float *u2_  = nltvl1w->u2_;
  
  float *v1 = nltvl1w->v1;
  float *v2 = nltvl1w->v2;

  float *rho_c = nltvl1w->rho_c;
  float *grad  = nltvl1w->grad;

  float *u1_tmp = nltvl1w->u1_tmp;
  float *u2_tmp = nltvl1w->u2_tmp;

  float *I1x = nltvl1w->I1x;
  float *I1y = nltvl1w->I1y;

  float *I1w = nltvl1w->I1w;
  float *I1wx = nltvl1w->I1wx;
  float *I1wy = nltvl1w->I1wy;
  
  //Divergence
  float *div_p = nltvl1w->div_p;
  float *div_q = nltvl1w->div_q;

  const int n_d = NL_DUAL_VAR;
  const float l_t = lambda * theta;

    //TODO:Weights
  const int iiw = nltvl1w->iiw;
  const int ijw = nltvl1w->ijw;
  float *weight = nltvl1w->weight;

  for (int warpings = 0; warpings < warps; warpings++)
  {   
    // compute the warping of the Right image and its derivatives Ir(x + u1o), xIrx (x + u1o) and Iry (x + u2o)
    bicubic_interpolation_warp_patch(I1,  u1, u2, I1w, 
                              ii, ij, ei, ej, w, h, false);
    bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx, 
                              ii, ij, ei, ej, w, h, false);
    bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy, 
                                ii, ij, ei, ej, w, h, false);

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++)
    for (int k = ii; k < ei; k++)
    {
      const int i = l*w + k;
      const float Ix2 = I1wx[i] * I1wx[i];
      const float Iy2 = I1wy[i] * I1wy[i];

      // store the |Grad(I1)|^2
      grad[i] = (Ix2 + Iy2);

      // compute the constant part of the rho function
      rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
            - I1wy[i] * u2[i] - I0[i]);
    }

    //Get the correct wt to force than the sum will be 1
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++)
    for (int k = ii; k < ei; k++)
    {
      float wt1_tmp = 0.0;
      float wt2_tmp = 0.0;    
      const int i = l*w + k;
      for (int j = 0; j < n_d; j++)
      {
        const int ap1i = p[i].api[j];
        const int ap1j = p[i].apj[j];
        const int ap2i = q[i].api[j];
        const int ap2j = q[i].apj[j];
        
        assert(ap1i == ap2i && ap1j == ap2j);
        //The position should be the same
        const int ap1 = validate_ap_patch(ii, ij, ei, ej, ap1i, ap1j);
        const int ap2 = validate_ap_patch(ii, ij, ei, ej, ap2i, ap2j);
        assert(ap1 == ap2);
        if ((ap1 == 0) && (ap2 == 0))
        {
          const float wp1 = p[i].wp[j]; 
          const float wp2 = q[i].wp[j]; 
          assert(wp1 == wp2);
          assert(wp1 >= 0);
          assert(wp2 >= 0);
          wt1_tmp +=wp1;
          wt2_tmp +=wp2;
        }
      }
      p[i].wt = wt1_tmp;
      q[i].wt = wt2_tmp;
    }

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++)
    for (int k = ii; k < ei; k++)
    {
      const int i = l*w + k;
      u1_[i] = u1[i];
      u2_[i] = u2[i];
    }

    int n = 0;
    float err_D = INFINITY;
    while (err_D > tol_OF*tol_OF && n < ofD->params.max_iter_patch)
    {
      n++;
      // estimate the values of the variable (v1, v2)
      // (thresholding opterator TH)
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const float l_t_w = l_t * weight[l-ij + ijw]*weight[k-ii + iiw];
        const int i = l*w + k;
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
      nltvl1_w_getD(u1_, u2_, ii, ij, ei, ej, w, n_d, tau, p, q);
      //Almacenamos la iteracion anterior
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*w + k;
        u1_tmp[i] = u1[i];
        u2_tmp[i] = u2[i];    
      }
      }

      //Primal variables
      non_local_divergence(p, ii, ij, ei, ej, w, n_d, div_p);
      non_local_divergence(q, ii, ij, ei, ej, w, n_d, div_q);
      nltvl1_w_getP(v1, v2, div_p, div_q, theta, tau,
                        ii, ij, ei, ej, w, u1, u2, &err_D);

      //aceleration = 1
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
      for (int l = ij; l < ej; l++){
      for (int k = ii; k < ei; k++){
        const int i = l*w + k;
        u1_[i] = 2*u1[i] - u1_tmp[i];
        u2_[i] = 2*u2[i] - u2_tmp[i];
      }
      }
    }
    if (verbose)
      fprintf(stderr, "Warping: %d,Iter: %d "
      "Error: %f\n", warpings,n, err_D);
  }
  eval_nltvl1_w(I0, I1, ofD, nltvl1w, ener_N, ii, ij, ei, ej, lambda, theta, w, h);

}

#endif //NLTVL1_MODEL_W
