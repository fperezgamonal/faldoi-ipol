#ifndef NLTVCSAD_MODEL
#define NLTVCSAD_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm>

#include "energy_structures.h"
#include "aux_energy_model.h"
extern "C" {
#include "bicubic_interpolation.h"
}

void  initialize_stuff_nltvcsad(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        const int w,
        const int h)
{
//    const int w = ofCore->w;
//    const int h = ofCore->h;
    ofStuff->nltvcsad.p    = new DualVariables[w*h];
    ofStuff->nltvcsad.q    = new DualVariables[w*h];
    ofStuff->nltvcsad.pnei = new PosNei[w*h];
    ofStuff->nltvcsad.v1 =  new float[w*h];
    ofStuff->nltvcsad.v2 =  new float[w*h];
    ofStuff->nltvcsad.rho_c =  new float[w*h];
    ofStuff->nltvcsad.grad =  new float[w*h];
    ofStuff->nltvcsad.u1_ =  new float[w*h];
    ofStuff->nltvcsad.u2_ =  new float[w*h];
    ofStuff->nltvcsad.u1_tmp = new float[w*h];
    ofStuff->nltvcsad.u2_tmp = new float[w*h];
    ofStuff->nltvcsad.I1x = new float[w*h];
    ofStuff->nltvcsad.I1y = new float[w*h];
    ofStuff->nltvcsad.I1w = new float[w*h];
    ofStuff->nltvcsad.I1wx = new float[w*h];
    ofStuff->nltvcsad.I1wy = new float[w*h];
    ofStuff->nltvcsad.div_p = new float[w*h];
    ofStuff->nltvcsad.div_q = new float[w*h];
}

void  free_stuff_nltvcsad(SpecificOFStuff *ofStuff)

{

    delete [] ofStuff->nltvcsad.p;
    delete [] ofStuff->nltvcsad.q;
    delete [] ofStuff->nltvcsad.pnei;
    delete [] ofStuff->nltvcsad.v1;
    delete [] ofStuff->nltvcsad.v2;
    delete [] ofStuff->nltvcsad.rho_c;
    delete [] ofStuff->nltvcsad.grad;
    delete [] ofStuff->nltvcsad.u1_;
    delete [] ofStuff->nltvcsad.u2_;
    delete [] ofStuff->nltvcsad.u1_tmp;
    delete [] ofStuff->nltvcsad.u2_tmp;
    delete [] ofStuff->nltvcsad.I1x;
    delete [] ofStuff->nltvcsad.I1y;
    delete [] ofStuff->nltvcsad.I1w;
    delete [] ofStuff->nltvcsad.I1wx;
    delete [] ofStuff->nltvcsad.I1wy;
    delete [] ofStuff->nltvcsad.div_p;
    delete [] ofStuff->nltvcsad.div_q;
}




void eval_nltvcsad(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        NonLocalTvCsadStuff *nltvcsad,
        float *ener_N,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const float lambda,  // weight of the data term
        const float theta,
        const int w,
        const int h
)
{

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;


    //Columns and Rows
//    const int w = ofD->w;
//    const int h = ofD->h;

    float *I1w = nltvcsad->I1w;
    DualVariables *p = nltvcsad->p;
    DualVariables *q = nltvcsad->q;
    PosNei *pnei  = nltvcsad->pnei;


    float *v1 = nltvcsad->v1;
    float *v2 = nltvcsad->v2;

    float ener = 0.0;
    int n_d = NL_DUAL_VAR;


    int ndt = DT_NEI;


    bicubic_interpolation_warp_patch(I1,  u1, u2, I1w,
                                     ii, ij, ei, ej, w, h, false);
    //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
    int m  = 0;
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int i = l*w + k;

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
                if ((ap1==0) && (ap2==0))
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
            assert(p[i].wt == q[i].wt);
            float dt = 0.0;
            for (int j = 0; j < ndt; j++)
            {
                const int api = pnei[i].api[j];
                const int apj = pnei[i].apj[j];
                const int ap = validate_ap_patch(ii, ij, ei, ej, api, apj);
                if (ap == 0)
                {
                    // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, pnei[i].apj[j]*nx + pnei[i].api[j], pnei[i].apj[j], pnei[i].api[j]);
                    assert(pnei[i].api[j] >= 0);
                    assert(pnei[i].apj[j] >= 0);
                    assert(pnei[i].apj[j]*w + pnei[i].api[j] < w*h);
                    const int pos = apj*w + api;
                    dt +=  fabs(I0[i] - I0[pos] - I1w[i] + I1w[pos]);
                }
            }
            dt *=lambda;
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

/*
 * - Name: getP

 *
*/

void nltvcsad_getP(
        float *v1,
        float *v2,
        float *div_p1,
        float *div_p2,
        int *mask,
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


void nltvcsad_getD(
        const float *u1,
        const float *u2,
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
                if ((ap1==0) && (ap2==0))
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

void guided_nltvcsad(
        const float *I0,           // source image
        const float *I1,           // target image
        OpticalFlowData *ofD,
        NonLocalTvCsadStuff *nltvcsad,
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
)
{
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int *mask = ofD->fixed_points;
    //Columns and Rows
//    const int w = ofD->w;
//    const int h = ofD->h;

    DualVariables *p = nltvcsad->p;
    DualVariables *q = nltvcsad->q;
    PosNei     *pnei = nltvcsad->pnei;

    float *u1_  = nltvcsad->u1_;
    float *u2_  = nltvcsad->u2_;

    float *v1 = nltvcsad->v1;
    float *v2 = nltvcsad->v2;

    float *grad  = nltvcsad->grad;

    float *u1_tmp = nltvcsad->u1_tmp;
    float *u2_tmp = nltvcsad->u2_tmp;

    float *I1x = nltvcsad->I1x;
    float *I1y = nltvcsad->I1y;

    float *I1w = nltvcsad->I1w;
    float *I1wx = nltvcsad->I1wx;
    float *I1wy = nltvcsad->I1wy;

    //Divergence
    float *div_p = nltvcsad->div_p;
    float *div_q = nltvcsad->div_q;

    const int n_d = NL_DUAL_VAR;
    const int ndt = DT_NEI;
    const float l_t = lambda * theta;


    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp_patch(I1,  u1, u2, I1w,
                                         ii, ij, ei, ej, w, h, false);
        bicubic_interpolation_warp_patch(I1x, u1, u2, I1wx,
                                         ii, ij, ei, ej, w, h, false);
        bicubic_interpolation_warp_patch(I1y, u1, u2, I1wy,
                                         ii, ij, ei, ej, w, h, false);

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){

                const int i = l*w + k;
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
                        assert(pnei[i].api[j] >= 0);
                        assert(pnei[i].apj[j] >= 0);
                        assert(pnei[i].apj[j]*w + pnei[i].api[j] < w*h);
                        const int pos = apj*w + api;

                        pnei[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                                        + I1wy[i] * u2[i])/grad[i];
                        n_tmp ++;
                    }
                }
                pnei[i].n= n_tmp;
            }
        }

//Get the correct wt to force than the sum will be 1
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
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
        }

//#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
                const int i = l*w + k;
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
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*w + k;
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
                    //TODO: Check minimization's integrity
                    v1[i] = u1[i] - I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
                    v2[i] = u2[i] - I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
                }
            }
            //Dual variables
            nltvcsad_getD(u1_, u2_, ii, ij, ei, ej, w, n_d, tau, p, q);
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
            nltvcsad_getP(v1, v2, div_p, div_q, mask, theta, tau,
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
            std::printf("Warping: %d,Iter: %d Error: %f\n", warpings,n, err_D);
    }
    eval_nltvcsad(I0, I1, ofD, nltvcsad, ener_N, ii, ij, ei, ej, lambda, theta, w, h);
}
#endif
