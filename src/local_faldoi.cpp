// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// All rights reserved.


#ifndef LOCAL_FALDOI
#define LOCAL_FALDOI

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <queue>
#include <random>
#include <future>
#include <algorithm>

#include "energy_structures.h"
#include "energy_model.h"
#include "utils_preprocess.h"

extern "C" {
#include "iio.h"
#include "bicubic_interpolation.h"
#include "elap_recsep.h"
}

#include <omp.h>

#include <iostream>
#include <fstream>
#include <string>

#include "utils.h"
#include "parameters.h"
#include <ctime>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////OUTLIERS FUNCTIONS/////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return INFINITY;
    return x[(i + j * w) * pd + l];
}

static int too_uniform(float *a, float tol, int i, int j, int w, int h, int pd) {
    float difference = 0;
    int neighbours[4][2] = {
            {0,  1},
            {0,  -1},
            {1,  0},
            {-1, 0}};
    for (int l = 0; l < pd; l++) {
        float center = getsample_inf(a, w, h, pd, i, j, l);
        for (int k = 0; k < 4; k++) {
            int px = i + neighbours[k][0];
            int py = j + neighbours[k][1];
            float neighborhood = getsample_inf(a, w, h, pd, px, py, l);
            if (isfinite(center) && isfinite(neighborhood)) {
                float tmp = abs(neighborhood - center);
                //printf("Tmp: %f, Tol: %f neighborhood: %f Center: %f\n", tmp, difference, neighborhood, center);
                if (difference < tmp) {
                    difference = tmp;
                }
            }
        }
    }

    if (difference < tol) {
        return 1;
    }
    return 0;
}

void too_uniform_areas(
        float *a,
        float *b,
        float *in0,
        int *trust_in0,
        int w,
        int h,
        float tol
) {

    auto *bw = new float[w * h];
    int size = w * h;
    int n = 0;

    bicubic_interpolation_warp(b, in0, in0 + size, bw, w, h, true);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            //If both areas present too uniform pixels, remove the flow.
            if ((too_uniform(a, tol, i, j, w, h, 1) == 1) || (too_uniform(bw, tol, i, j, w, h, 1) == 1)) {
                trust_in0[j * w + i] = 0;
            } else {
                trust_in0[j * w + i] = 1;
                n++;
            }
        }
    printf("Too-chosen: %f\n", (n * 1.0) / size);

    delete[] bw;
}

// Check forward-backward consistency check for of |u(x) + v(x+u(x))| < eps.
// Energy map related to that flows are put to INFINITY.
void fb_consistency_check(
        float *in0,
        float *in1,
        int *trust_in0,
        int w,
        int h,
        float epsilon
) {

    auto *u1w = new float[w * h];
    auto *u2w = new float[w * h];
    int size = w * h;
    int n = 0;

    bicubic_interpolation_warp(in1, in0, in0 + size, u1w, w, h, true);
    bicubic_interpolation_warp(in1 + w * h, in0, in0 + size, u2w, w, h, true);

    for (int i = 0; i < size; i++) {
        float tolerance = hypotf(in0[i] + u1w[i], in0[size + i] + u2w[i]);
        if (tolerance > epsilon) {
            //(tolerance > epsilon) means pixel is occluded
            trust_in0[i] = 0;
        } else {
            trust_in0[i] = 1;
            n++;
        }
    }
    printf("FB-Chosen: %f\n", (n * 1.0) / size);
    delete[] u2w;
    delete[] u1w;
}

void pruning_method(
        float *i0,          // I0
        float *i1,          // I1
        int w,              // width image
        int h,              // height image
        float *tol,         // tolerance too_uniform and f-b
        const int *method,  // if method[i]!=0, then
        int *trust_Go,      // energy map of u
        float *go,          // of to t, t+1
        int *trust_Ba,      // energy map of v
        float *ba           // of to t+1, t
) {
    auto *go_fb_check = new int[w * h];
    auto *go_cons_check = new int[w * h];
    auto *ba_fb_check = new int[w * h];
    auto *ba_cons_check = new int[w * h];

    for (int i = 0; i < w * h; i++) {
        // 0 - Invalid pixel 1 - Trustable pixel.
        trust_Go[i] = 1;
        trust_Ba[i] = 1;
    }


    // FB - consistency check
    if (method[0] == 1) {
        printf("FB-Consistency: %f\n", tol[0]);
        fb_consistency_check(go, ba, go_fb_check, w, h, tol[0]);
        fb_consistency_check(ba, go, ba_fb_check, w, h, tol[0]);
    }
    // Too-uniform consistency check
    if (method[1] == 1) {
        printf("Too Uniform -Consistency: %f\n", tol[1]);
        too_uniform_areas(i0, i1, go, go_cons_check, w, h, tol[1]);
        too_uniform_areas(i1, i0, ba, ba_cons_check, w, h, tol[1]);
    }
    for (int i = 0; i < w * h; i++) {
        if (method[0] == 1) {
            // FB-Consistency
            if (go_fb_check[i] == 0) {
                trust_Go[i] = 0;
            }
            if (ba_fb_check[i] == 0) {
                trust_Ba[i] = 0;
            }
        }
        // Too uniform -Consistency
        if (method[1] == 1) {
            if (go_cons_check[i] == 0) {
                trust_Go[i] = 0;
            }
            if (ba_cons_check[i] == 0) {
                trust_Ba[i] = 0;
            }
        }
    }

    delete[] go_fb_check;
    delete[] go_cons_check;
    delete[] ba_fb_check;
    delete[] ba_cons_check;
}

void delete_not_trustable_candidates(
        OpticalFlowData *ofD,
        float *in,
        float *ene_val,
        const int w,
        const int h
) {

    int *trust_points = ofD->trust_points;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float *chi = ofD->chi;
    // w, h as params in
    //int w = ofD->params.w;
    //int h = ofD->params.h;
    int n = 0;
    for (int i = 0; i < w * h; i++) {
        if (trust_points[i] == 0) {
            //printf("%f\n", ene_val[i]);
            if (ene_val[i] == 0.0) {
                n++;
            }
            in[i] = NAN;
            in[i + w * h] = NAN;
            u1[i] = NAN;
            u2[i] = NAN;
            ene_val[i] = INFINITY;
            // If the flow is non trustable, is considered
            // to be an occlusion
            chi[i] = 1;
        }
    }
    printf("Total_seeds: %d\n", n);
}

////////////////////////////////////////////////////////////////////////////////
//////////////////LOCAL PARTITION INTO SUBIMAGES////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void init_subimage_partitions(
        const float *i0,                                // Source image I0 (t)
        const float *i1,                                // Second image I1 (t+1)
        const float *i_1,                               // Previous image I-1 (t-1)
        const float *i2,                                // Third frame I2 (t+2)
        float *sal_go,                                  // Forward saliency map
        float *sal_ba,                                  // Backward saliency map
        const int w_src,                                // Width of the images
        const int h_src,                                // Height of the images
        const int h_parts,                              // Number of horizontal parts for the partition
        const int v_parts,                              // Number of vertical parts for the partition
        std::vector<PartitionData*> *p_data,            // Partition data
        Parameters params
) {
    // Define partition-specific variables
    int num_partitions = h_parts * v_parts;

    auto *sub_w = new int[num_partitions];
    auto *sub_h = new int[num_partitions];
    std::fill_n(sub_w, num_partitions, w_src / h_parts);
    std::fill_n(sub_h, num_partitions, h_src / v_parts);

    int rem_width = w_src % h_parts;
    int rem_height = h_src % v_parts;

    // Add extra pixels (if the partitions cannot be of equal size)

    // Horizontally
    if (rem_width > 0)
    {
        for (int i = 0; i < v_parts; i++)
        {
            sub_w[i * h_parts + h_parts-1] += rem_width;
        }
    }

    // Vertically
    if (rem_height > 0)
    {
        for (int j = 0; j < h_parts; j++)
        {
            sub_h[(v_parts-1) * h_parts + j] += rem_height;
        }
    }

    // Define partitions' offset w.r.t original image
    auto *off_x = new int[num_partitions];
    auto *off_y = new int[num_partitions];
    std::fill_n(off_x, num_partitions, 0);
    std::fill_n(off_y, num_partitions, 0);

    for (int v = 0; v < v_parts; v++)
        for (int h = 0; h < h_parts; h++)
        {
            if (h > 0)
            {
                off_x[v * h_parts + h] = h * sub_w[0];
                //std::cout << offset_x[v * h_parts + h] << std::endl;
            }
            if (v > 0)
            {
                off_y[v * h_parts + h] = v * sub_h[0];
                //std::cout << offset_y[v * h_parts + h] << std::endl;
            }
        }

    // No return value, just update the structs' fields via pointer
    for (int p = 0; p < num_partitions; p++)
    {
        auto *t_pdata = new PartitionData;    // Create new struct of partition data
        t_pdata->idx = p;                     // Assign idx
        t_pdata->width = sub_w[p];            // Assign width
        t_pdata->height = sub_h[p];           // Assign height
        t_pdata->off_x = off_x[p];            // Assign offset in x
        t_pdata->off_y = off_y[p];            // Assign offset in y
        t_pdata->oft0 = new float[sub_w[p] * sub_h[p] * 2];
        t_pdata->oft1 = new float[sub_w[p] * sub_h[p] * 2];
        t_pdata->ene_Go = new float[sub_w[p] * sub_h[p]];
        t_pdata->ene_Ba = new float[sub_w[p] * sub_h[p]];
        t_pdata->occ_Go = new float[sub_w[p] * sub_h[p]];
        t_pdata->occ_Ba = new float[sub_w[p] * sub_h[p]];
        t_pdata->sal_go = sal_go;
        t_pdata->sal_ba = sal_ba;

        p_data->push_back(t_pdata);
    }

    // Added the code from 'fill_subimage_partitions' here to only have two functions:
    // one that 'maps' subpartitions to the whole image and another that does the inverse

    //TODO: add other initializations (and necessary params as input): 'init_optical_flow, auxiliar_stuff,...'

    // Fill images

    const int n_channels = params.pd;

    // Check that dimensions match
    int total_size = 0;
    for (unsigned p = 0; p < num_partitions; p++)
    {
        total_size += p_data->at(p)->width * p_data->at(p)->height;
    }
    assert(total_size == w_src * h_src);

    for (unsigned p = 0; p < num_partitions; p++)
    {
        int size = p_data->at(p)->width * p_data->at(p)->height * n_channels;
//        //std::cout << "p = " << p << std::endl;
        auto *i0_p = new float[size];
        auto *i1_p = new float[size];
        auto *i_1_p = new float[size];
        auto *i2_p = new float[size];
        for (int k = 0; k < n_channels; k++)
        {
            //std::cout << "k = " << k << std::endl;
            for (int j = 0; j < p_data->at(p)->height; j++)
                for (int i = 0; i < p_data->at(p)->height; i++)
                {
                    int offset_p = 0;
                    if (p > 0)
                    {
                        offset_p = p * p_data->at(p-1)->width * p_data->at(p-1)->height;
                    }

                    int m = (j * p_data->at(p)->width + i) * n_channels + k;
                    // idx of the subarray is computed as follows (from image array):
                    // idx = (y + j) * p + x + i
                    // where:   x, y (subimage offsets, top-left corner)
                    //          i, j indices (for sub_width[p] cols, sub_height[p] rows)
                    //          p is the width of the image (q is the height)
                    // int idx = (offset_y[p] + j) * w_src + offset_x[p] + i; // equiv for 2D
                    int idx = ((p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i) * n_channels + k;  // " " 3D

                    i0_p[m] = i0[idx];
                    i1_p[m] = i1[idx];
                    i_1_p[m] = i_1[idx];
                    i2_p[m] = i2[idx];

                    //std::cout << partition[m] << "\t";
                    //if (i == p_data->at(p)->width-1)  std::cout << std::endl;

                }
        }
        //TODO: IT SEEMS TO WORK, DO THE SAME WITH ALL VARIABLES (see if here is the best place)
        p_data->at(p)->i0 = i0_p;
        p_data->at(p)->i1 = i1_p;
        p_data->at(p)->i_1 = i_1_p;
        p_data->at(p)->i2 = i2_p;

        p_data->at(p)->ofGo = init_Optical_Flow_Data(sal_go, params, sub_w[p], sub_h[p]);
        p_data->at(p)->ofBa = init_Optical_Flow_Data(sal_ba, params, sub_w[p], sub_h[p]);

        // Initialise Specific OF stuff
        initialize_auxiliar_stuff(p_data->at(p)->stuffGo, p_data->at(p)->ofGo, sub_w[p], sub_h[p]);
        initialize_auxiliar_stuff(p_data->at(p)->stuffBa, p_data->at(p)->ofBa, sub_w[p], sub_h[p]);

        // Prepare auxiliar stuff
        prepare_stuff(&p_data->at(p)->stuffGo, &p_data->at(p)->ofGo, &p_data->at(p)->stuffBa, &p_data->at(p)->ofBa,
                      i0_p, i1_p, i_1_p, i2_p, params.pd, &p_data->at(p)->i0n, &p_data->at(p)->i1n, &p_data->at(p)->i_1n,
                      &p_data->at(p)->i2n, sub_w[p], sub_h[p]);

        // Bilateral filter
        p_data->at(p)->BiFilt_Go = init_weights_bilateral(p_data->at(p)->i0n, sub_w[p], sub_h[p]);
        p_data->at(p)->BiFilt_Ba = init_weights_bilateral(p_data->at(p)->i1n, sub_w[p], sub_h[p]);


    }

    //TODO: add yet another loop for initialization (see if we can reduce it afterwards)

//    for (unsigned p = 0; p < num_partitions; p++)
//    {
//        // Initialise OF data
//        p_data->at(p)->ofGo = init_Optical_Flow_Data(sal_go, params, sub_w[p], sub_h[p]);
//        p_data->at(p)->ofBa = init_Optical_Flow_Data(sal_ba, params, sub_w[p], sub_h[p]);
//
//        // Initialise Specific OF stuff
//        initialize_auxiliar_stuff(p_data->at(p)->stuffGo, p_data->at(p)->ofGo, sub_w[p], sub_h[p]);
//        initialize_auxiliar_stuff(p_data->at(p)->stuffBa, p_data->at(p)->ofBa, sub_w[p], sub_h[p]);
//
//        // Prepare auxiliar stuff
//        prepare_stuff(&p_data->at(p)->stuffGo, &p_data->at(p)->ofGo, &p_data->at(p)->stuffBa, &p_data->at(p)->ofBa,
//                      i0, i1, i_1, i2, params.pd, &p_data->at(p)->i0n, &p_data->at(p)->i1n, &p_data->at(p)->i_1n,
//                      &p_data->at(p)->i2n, sub_w[p], sub_h[p]);
//    }
}

void update_tvl2_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {
    if (n_ch == 1)
    {
        // Chi
        p_data->ofGo.chi[idx_out] = ofGo->chi[idx_src];
        p_data->ofBa.chi[idx_out] = ofBa->chi[idx_src];

        // Fixed points and trust points
        p_data->ofGo.fixed_points[idx_out] = ofGo->fixed_points[idx_src];
        p_data->ofBa.fixed_points[idx_out] = ofBa->fixed_points[idx_src];
        p_data->ofGo.trust_points[idx_out] = ofGo->trust_points[idx_src];
        p_data->ofBa.trust_points[idx_out] = ofBa->trust_points[idx_src];

        // Saliency
        p_data->ofGo.saliency[idx_out] = ofGo->saliency[idx_src];
        p_data->ofBa.saliency[idx_out] = ofBa->saliency[idx_src];

    }
    else
    {
        // OF fields
        p_data->ofGo.u1[idx_out] = ofGo->u1[idx_src];
        p_data->ofGo.u2[idx_out] = ofGo->u2[idx_src];
        p_data->ofBa.u1[idx_out] = ofBa->u1[idx_src];
        p_data->ofBa.u2[idx_out] = ofBa->u2[idx_src];

        p_data->ofGo.u1_ba[idx_out] = ofGo->u1_ba[idx_src];
        p_data->ofGo.u2_ba[idx_out] = ofGo->u2_ba[idx_src];
        p_data->ofBa.u1_ba[idx_out] = ofBa->u1_ba[idx_src];
        p_data->ofBa.u2_ba[idx_out] = ofBa->u2_ba[idx_src];

        // Filters
        p_data->ofGo.u1_filter[idx_out] = ofGo->u1_filter[idx_src];
        p_data->ofGo.u2_filter[idx_out] = ofGo->u2_filter[idx_src];
        p_data->ofBa.u1_filter[idx_out] = ofBa->u1_filter[idx_src];
        p_data->ofBa.u2_filter[idx_out] = ofBa->u2_filter[idx_src];
    }

}

void update_tvl2_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {
    // Xi
    p_data->stuffGo.tvl2.xi11[idx_out] = stuffGo->tvl2.xi11[idx_src];
    p_data->stuffGo.tvl2.xi12[idx_out] = stuffGo->tvl2.xi11[idx_src];
    p_data->stuffGo.tvl2.xi21[idx_out] = stuffGo->tvl2.xi21[idx_src];
    p_data->stuffGo.tvl2.xi22[idx_out] = stuffGo->tvl2.xi22[idx_src];

    p_data->stuffBa.tvl2.xi11[idx_out] = stuffBa->tvl2.xi11[idx_src];
    p_data->stuffBa.tvl2.xi12[idx_out] = stuffBa->tvl2.xi11[idx_src];
    p_data->stuffBa.tvl2.xi21[idx_out] = stuffBa->tvl2.xi21[idx_src];
    p_data->stuffBa.tvl2.xi22[idx_out] = stuffBa->tvl2.xi22[idx_src];

    // u1, u2
    p_data->stuffGo.tvl2.u1x[idx_out] = stuffGo->tvl2.u1x[idx_src];
    p_data->stuffGo.tvl2.u1y[idx_out] = stuffGo->tvl2.u1y[idx_src];
    p_data->stuffGo.tvl2.u2x[idx_out] = stuffGo->tvl2.u2x[idx_src];
    p_data->stuffGo.tvl2.u2y[idx_out] = stuffGo->tvl2.u2y[idx_src];

    p_data->stuffBa.tvl2.u1x[idx_out] = stuffBa->tvl2.u1x[idx_src];
    p_data->stuffBa.tvl2.u1y[idx_out] = stuffBa->tvl2.u1y[idx_src];
    p_data->stuffBa.tvl2.u2x[idx_out] = stuffBa->tvl2.u2x[idx_src];
    p_data->stuffBa.tvl2.u2y[idx_out] = stuffBa->tvl2.u2y[idx_src];

    // v1, v2 (auxiliar minimization variables)
    p_data->stuffGo.tvl2.v1[idx_out] = stuffGo->tvl2.v1[idx_src];
    p_data->stuffGo.tvl2.v2[idx_out] = stuffGo->tvl2.v2[idx_src];

    p_data->stuffBa.tvl2.v1[idx_out] = stuffBa->tvl2.v1[idx_src];
    p_data->stuffBa.tvl2.v2[idx_out] = stuffBa->tvl2.v2[idx_src];

    // Auxiliary (gradients, weighted gradients, divergence, ...)
    p_data->stuffGo.tvl2.rho_c[idx_out] = stuffGo->tvl2.rho_c[idx_src];
    p_data->stuffBa.tvl2.rho_c[idx_out] = stuffBa->tvl2.rho_c[idx_src];

    p_data->stuffGo.tvl2.grad[idx_out] = stuffGo->tvl2.grad[idx_src];
    p_data->stuffBa.tvl2.grad[idx_out] = stuffBa->tvl2.grad[idx_src];

    p_data->stuffGo.tvl2.u1_[idx_out] = stuffGo->tvl2.u1_[idx_src];
    p_data->stuffGo.tvl2.u2_[idx_out] = stuffGo->tvl2.u2_[idx_src];

    p_data->stuffBa.tvl2.u1_[idx_out] = stuffBa->tvl2.u1_[idx_src];
    p_data->stuffBa.tvl2.u2_[idx_out] = stuffBa->tvl2.u2_[idx_src];

    p_data->stuffGo.tvl2.u1Aux[idx_out] = stuffGo->tvl2.u1Aux[idx_src];
    p_data->stuffGo.tvl2.u2Aux[idx_out] = stuffGo->tvl2.u2Aux[idx_src];

    p_data->stuffBa.tvl2.u1Aux[idx_out] = stuffBa->tvl2.u1Aux[idx_src];
    p_data->stuffBa.tvl2.u2Aux[idx_out] = stuffBa->tvl2.u2Aux[idx_src];

    p_data->stuffGo.tvl2.I1x[idx_out] = stuffGo->tvl2.I1x[idx_src];
    p_data->stuffGo.tvl2.I1y[idx_out] = stuffGo->tvl2.I1y[idx_src];

    p_data->stuffBa.tvl2.I1x[idx_out] = stuffBa->tvl2.I1x[idx_src];
    p_data->stuffBa.tvl2.I1y[idx_out] = stuffBa->tvl2.I1y[idx_src];

    p_data->stuffGo.tvl2.I1wx[idx_out] = stuffGo->tvl2.I1wx[idx_src];
    p_data->stuffGo.tvl2.I1wy[idx_out] = stuffGo->tvl2.I1wy[idx_src];

    p_data->stuffBa.tvl2.I1wx[idx_out] = stuffBa->tvl2.I1wx[idx_src];
    p_data->stuffBa.tvl2.I1wy[idx_out] = stuffBa->tvl2.I1wy[idx_src];

    p_data->stuffGo.tvl2.div_xi1[idx_out] = stuffGo->tvl2.div_xi1[idx_src];
    p_data->stuffGo.tvl2.div_xi2[idx_out] = stuffGo->tvl2.div_xi2[idx_src];

    p_data->stuffBa.tvl2.div_xi1[idx_out] = stuffBa->tvl2.div_xi1[idx_src];
    p_data->stuffBa.tvl2.div_xi2[idx_out] = stuffBa->tvl2.div_xi2[idx_src];

    p_data->stuffGo.tvl2.u_N[idx_out] = stuffGo->tvl2.u_N[idx_src];
    p_data->stuffBa.tvl2.u_N[idx_out] = stuffBa->tvl2.u_N[idx_src];

}

// TODO: implement the rest once it works with the TVL2 functional
void update_tvl2w_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_tvl2w_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}


void update_tvl2occ_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_tvl2occ_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_nltvl1_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_nltvl1_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_nltv1w_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_nltv1w_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_tvcsad_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_tvcsad_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_tvcsadw_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_tvcsadw_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_nltvcsad_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_nltvcsad_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}

void update_nltvcsadw_of_data(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out,
        const int n_ch
) {

}
void update_nltvcsadw_stuffof(
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int idx_src,
        const int idx_out
) {

}


void update_partitions_structures(
        OpticalFlowData* ofGo,
        OpticalFlowData* ofBa,
        SpecificOFStuff* stuffGo,
        SpecificOFStuff* stuffBa,
        PartitionData *p_data,
        const int n_ch,
        const int idx_src,
        const int idx_out
) {
    switch(ofGo->params.val_method)
    {
        case M_NLTVL1:          // NLTV-L1
            update_nltvl1_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_nltvl1_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_TVCSAD:          // TV-CSAD
            update_tvcsad_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_tvcsad_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_NLTVCSAD:        // NLTV-CSAD
            update_nltvcsad_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_nltvcsad_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_TVL1_W:          // TV-l2 coupled with weights
            update_tvl2w_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_tvl2w_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_NLTVCSAD_W:      // NLTV-CSAD with weights
            update_nltvcsadw_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_nltvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_NLTVL1_W:        // NLTV-L1 with weights
            update_nltv1w_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_nltv1w_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_TVCSAD_W:        // TV-CSAD with weights
            update_tvcsadw_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_tvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        case M_TVL1_OCC:        // TV-l2 with occlusion
            update_tvl2occ_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_tvl2occ_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
        default:                // TV-l2 coupled
            update_tvl2_of_data(ofGo, ofBa, p_data, idx_src, idx_out, n_ch);
            update_tvl2_stuffof(stuffGo, stuffBa, p_data, idx_src, idx_out);
            break;
    }
}
void image_to_partitions(
        const float *oft0,
        const float *oft1,
        const float *ene_Go,
        const float *ene_Ba,
        const float *occ_Go,
        const float *occ_Ba,
        OpticalFlowData *ofGo,
        OpticalFlowData *ofBa,
        SpecificOFStuff *stuffGo,
        SpecificOFStuff *stuffBa,
        const int n_partitions,
        const int w_src,
        const int h_src,
        std::vector<PartitionData*> *p_data
) {

    // Update variables with one channel (energy, occlusions, some OF and SpecificOF terms)
    int n_channels = 2;
    for (unsigned p = 0; p < n_partitions; p++)
    {
        int size = p_data->at(p)->width * p_data->at(p)->height;
        int size2 = size * n_channels;

        // Define temporal variables
        auto *ene_Go_p = new float[size];
        auto *ene_Ba_p = new float[size];
        auto *occ_Go_p = new float[size];
        auto *occ_Ba_p = new float[size];
        auto *oft0_p = new float[size2];
        auto *oft1_p = new float[size2];

        for (int k = 0; k < n_channels; k++)
            for (int j = 0; j < p_data->at(p)->height; j++)
                for (int i = 0; i < p_data->at(p)->height; i++)
                {
                    int offset_p = 0;
                    if (p > 0)
                    {
                        offset_p = p * p_data->at(p-1)->width * p_data->at(p-1)->height;
                    }

                    // 'Mapping' indices
                    int m = (j * p_data->at(p)->width + i) * n_channels + k;
                    int idx = ((p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i) * n_channels + k + offset_p;

                    // Update everything here
                    if (k <= n_channels-1) {
                        ene_Go_p[m] = ene_Go[idx];
                        ene_Ba_p[m] = ene_Ba[idx];
                        occ_Go_p[m] = occ_Go[idx];
                        occ_Ba_p[m] = occ_Ba[idx];
                        oft0_p[m] = oft0[idx];
                        oft1_p[m] = oft1[idx];
                        // Params
                        p_data->at(p)->ofGo.params = ofGo->params;
                        p_data->at(p)->ofBa.params = ofBa->params;
                    }
                    else
                    {
                        // Only update variables with more than a channel
                        oft0_p[m] = oft0[idx];
                        oft1_p[m] = oft1[idx];
                    }
                    // Update structures in a separate function (more clean):
                    // Only update the values for the functional used (faster)
                    // Being independent of the nº of channels, we can call it outside the if
                    update_partitions_structures(ofGo, ofBa, stuffGo, stuffBa, p_data->at(p), k, idx, m);


                }
        p_data->at(p)->ene_Go = ene_Go_p;
        p_data->at(p)->ene_Ba = ene_Ba_p;
        p_data->at(p)->occ_Go = occ_Go_p;
        p_data->at(p)->occ_Ba = occ_Ba_p;
        p_data->at(p)->oft0 = oft0_p;
        p_data->at(p)->oft1 = oft1_p;
    }
}


void fill_subimage_partition(
        const float *src_img,       // Pointer to source image (e.g.: i0n, i1n, ...)
        int w_src,                  // Width of the source image
        int n_channels,             // Number of channels (for flow array w x h x 2)
        int h_parts,                // Number of horizontal parts for the partition
        int v_parts,                // Number of vertical parts for the partition
        std::vector<PartitionData*> *p_data
) {
    // Define partition-specific variables
    int num_partitions = h_parts * v_parts;

    // Create matrix to store 'num_partitions' partitions
    int subimgs_size = 0;
    for (unsigned p = 0; p < num_partitions; p++) {
        subimgs_size += p_data->at(p)->width * p_data->at(p)->height;
    }


    //auto *subimages = new float[subimgs_size * n_channels];

    for (unsigned p = 0; p < num_partitions; p++) {
        std::cout << "p = " << p << std::endl;
        auto *partition = new float[p_data->at(p)->width * p_data->at(p)->height * n_channels];
        for (int k = 0; k < n_channels; k++) {
            std::cout << "k = " << k << std::endl;
            for (int j = 0; j < p_data->at(p)->height; j++)
                for (int i = 0; i < p_data->at(p)->height; i++) {
                    //int m = j * sub_width[p] + i;   // equiv for 2D mtx
                    int offset_p = 0;
                    if (p > 0) {
                        offset_p = p * p_data->at(p - 1)->width * p_data->at(p - 1)->height;
                    }

                    //int m = ((j + p_data->at(p)->off_y) * p_data->at(p)->width +
                    //        p_data->at(p)->off_x + i) * n_channels + k + offset_p;  // "" 3D "
                    int m = (j * p_data->at(p)->width + i) * n_channels + k;
                    // idx of the subarray is computed as follows (from image array):
                    // idx = (y + j) * p + x + i
                    // where:   x, y (subimage offsets, top-left corner)
                    //          i, j indices (for sub_width[p] cols, sub_height[p] rows)
                    //          p is the width of the image (q is the height)
                    // int idx = (offset_y[p] + j) * w_src + offset_x[p] + i; // equiv for 2D
                    int idx =
                            ((p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i) * n_channels + k;  // " " 3D

                    //subimages[m] = src_img[idx];
                    //p_data->at(p)->i0n[m] = src_img[idx];
                    partition[m] = src_img[idx];
                    std::cout << partition[m] << "\t";
                    if (i == p_data->at(p)->width - 1) std::cout << std::endl;

                }
        }
        //TODO: IT SEEMS TO WORK, DO THE SAME WITH ALL VARIABLES (see if here is the best place)
        p_data->at(p)->i0n = partition;


    }
}

void get_subimage_partition(
        const float *subimages,
        float **part,
        int part_idx,
        const int *sub_width,
        const int *sub_height,
        const int n_channels,
        const int *offset_x,
        const int *offset_y
) {
    //std::cout << "Partition num = "  << part_idx << std::endl;
    auto *partition = new float[sub_width[part_idx] * sub_height[part_idx] * n_channels];
    for (int k = 0; k < n_channels; k++)
    {
        //std::cout << "k = " << k << std::endl;
        for (int j = 0; j < sub_height[part_idx]; j++)
            for (int i = 0; i < sub_width[part_idx]; i++)
            {
                int offset_p = 0;
                if (part_idx > 0)
                {
                    offset_p = part_idx * sub_width[part_idx-1] * sub_height[part_idx-1];
                }

                int m = ((j + offset_y[part_idx]) * sub_width[part_idx] + offset_x[part_idx] + i) * n_channels + k + offset_p;  // "" 3D "
                int m_p = (j * sub_width[part_idx] + i) * n_channels + k;
                partition[m_p] = subimages[m];
                //std::cout << partition[m_p] << "\t";
                //if (i == sub_width[part_idx]-1)  std::cout << std::endl;

            }
    }

    // Return pointer to partition 'part_idx' array
    *part = partition;
}

//TODO: reduce code of partitions by using solely this function (after the current code works)
/*
void image_variables_to_partitions(
        struct PartitionParams* params,     // Struct with size, offset and number of slices
        const int w_src,                    // Width of the source image
        const int h_src,                    // Height of the source image
        const float *i0,
        const float *i1,
        const float *i_1,
        const float *i2,
        const float *oft0,
        const float *oft1,
        const float *ene_Go,
        const float *ene_Ba,
        const float *occ_Go,
        const float *occ_Ba,
        std::vector<float *> *i0_p,
        std::vector<float *> *i1_p,
        std::vector<float *> *i_1_p,
        std::vector<float *> *i2_p,
        std::vector<float *> *oft0_p,
        std::vector<float *> *oft1_p,
        std::vector<float *> *ene_Go_p,
        std::vector<float *> *ene_Ba_p,
        std::vector<float *> *occ_Go_p,
        std::vector<float *> *occ_Ba_p,
        std::vector<OpticalFlowData> ofGo_p,
        std::vector<OpticalFlowData> ofBa_p,
        std::vector<SpecificOFStuff> stuffGo_p,
        std::vector<SpecificOFStuff> stuffBa_p
) {
    int n_channels = 3;

    // Define partition-specific variables
    int num_partitions = params->h_parts * params->v_parts;

    // Create matrix to store num_partitions partitions
    // mtx has sum_{i in num_partitions} sub_width[i] * sub_height[i] which always should be equal to w_src * h_src!
    int subimgs_size[num_partitions] = {0};
    int total_size = 0;
    for (int p = 0; p < num_partitions; p++)
    {
        subimgs_size[p] = params->sub_width[p] * params->sub_height[p];
        total_size += subimgs_size[p];
    }

    // Sanity check
    assert(total_size == w_src * h_src);

    for (int p = 0; p < num_partitions; p++)
    {
        //std::cout << "p = " << p << std::endl;
        // Define temporal variables
        int size = subimgs_size[p] * n_channels;            // 3 channels
        auto *i0_p_tmp = new float[size];
        auto *i1_p_tmp = new float[size];
        auto *i_1_p_tmp = new float[size];
        auto *i2_p_tmp = new float[size];

        int size2 = subimgs_size[p] * (n_channels - 1);     // 2 channels
        auto *oft0_p_tmp = new float[size2];
        auto *oft1_p_tmp = new float[size2];

        int size3 = subimgs_size[p] * (n_channels - 2);     // 1 channel
        auto *ene_Go_p_tmp = new float[size3];
        auto *ene_Ba_p_tmp = new float[size3];
        auto *occ_Go_p_tmp = new float[size3];
        auto *occ_Ba_p_tmp = new float[size3];



        for (int k = 0; k < n_channels; k++)
        {
            //std::cout << "k = " << k << std::endl;
            for (int j = 0; j < params->sub_height[p]; j++)
                for (int i = 0; i < params->sub_width[p]; i++)
                {
                    int m = (j * params->sub_width[p] + i) * n_channels + k;
                    int idx = ((params->offset_y[p] + j) * w_src + params->offset_x[p] + i) * n_channels + k;  // " " 3D

                    // Assign pointers
                    // Frames
                    i0_p_tmp[m] = i0[idx];
                    i1_p_tmp[m] = i1[idx];
                    i_1_p_tmp[m] = i_1[idx];
                    i2_p_tmp[m] = i2[idx];
                    oft0_p_tmp[m] = oft0[idx];
                    oft1_p_tmp[m] = oft1[idx];
                    ene_Go_p_tmp[m] = ene_Go[idx];
                    ene_Ba_p_tmp[m] = ene_Ba[idx];
                    occ_Go_p_tmp[m] = occ_Go[idx];
                    occ_Ba_p_tmp[m] = occ_Ba[idx];

                    //std::cout << i0_p_tmp[m] << "\t";
                    //if (i == params->sub_width[p]-1)  std::cout << std::endl;

                }
        }
        i0_p->push_back(i0_p_tmp);

        //push back everything else (see auto's)

    }

}
*/

////////////////////////////////////////////////////////////////////////////////
//////////////////LOCAL INITIALIZATION//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Poisson Interpolation
void interpolate_poisson(
        OpticalFlowData *ofD,
        const PatchIndexes &patch
) {
    int w = patch.ei - patch.ii;
    int h = patch.ej - patch.ij;

    int *fixed_points = ofD->fixed_points;
    int wR = ofD->params.w;

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;

    float buf_in[2 * MAX_PATCH * MAX_PATCH];
    float buf_out[2 * MAX_PATCH * MAX_PATCH];
    assert(w * h < MAX_PATCH * MAX_PATCH);

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            int x = i + patch.ii;
            int y = j + patch.ij;
            int xy = y * wR + x;
            // 1 fixed - 0 not
            if (fixed_points[xy] == 1) {
                buf_in[j * w + i] = u1[xy];
                buf_in[j * w + i + w * h] = u2[xy];
            } else { // not fixed
                buf_in[j * w + i] = NAN;
                buf_in[j * w + i + w * h] = NAN;
            }
        }
    }

    elap_recursive_separable(buf_out, buf_in, w, h, 2, 0.4, 3, 7);

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            int x = i + patch.ii;
            int y = j + patch.ij;
            int xy = y * wR + x;
            u1[xy] = buf_out[j * w + i];
            u2[xy] = buf_out[j * w + i + w * h];
        }
    }
}


void bilateral_filter(
        OpticalFlowData *ofD,
        BilateralFilterData *BiFilt,
        const PatchIndexes &patch,
        const int w,
        const int h
) {

    int *trust_points = ofD->trust_points;
    int *fixed_points = ofD->fixed_points;

    // Added changes for subimages
    //const int w = ofD->params.w;
    //const int h = ofD->params.h;

    const int w_patch = patch.ei - patch.ii;
    const int h_patch = patch.ej - patch.ij;

    const int wr_patch = ofD->params.w_radio;
    const int wr_filter = PATCH_BILATERAL_FILTER;
    const int wr = wr_filter + wr_patch;

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;

    float *u1_filter = ofD->u1_filter;
    float *u2_filter = ofD->u2_filter;

    PatchIndexes area_interp = get_index_patch(wr, w, h, patch.i, patch.j, 1);

    // Copy all flow values to minimize and surroundings
    for (int j = 0; j < area_interp.ej - area_interp.ij; j++) {
        for (int i = 0; i < area_interp.ei - area_interp.ii; i++) {


            // Coordinates of pixel over whole image
            const int x = i + area_interp.ii;
            const int y = j + area_interp.ij;
            const int xy = y * w + x;

            // Initialize flow in patch for filtering
            // We use values of trust points and fixed points
            if (trust_points[xy] == 1 || fixed_points[xy] == 1) {
                u1_filter[xy] = u1[xy];
                u2_filter[xy] = u2[xy];
            } else {
                u1_filter[xy] = 0.0;
                u2_filter[xy] = 0.0;

            }
        }
    }


    /*
     * For each pixel in the patch non trustable, do
     * bilateral filtering
     */
    int iter = ITER_BILATERAL_FILTER;
    for (int it = 0; it < iter; it++) {

        for (int j = 0; j < h_patch; j++) {
            for (int i = 0; i < w_patch; i++) {

                // Coordinates of pixel over whole image
                int x = i + patch.ii;
                int y = j + patch.ij;

                int xy = y * w + x;
                // If pixel has not survived the previous prunning or
                // if it has not been fixed, we make interpolation
                if (trust_points[xy] == 0 && fixed_points[xy] == 0) {

                    // Index of points around ij
                    const PatchIndexes index_interp = BiFilt->indexes_filtering[xy];


                    // Variable that contains the precalculated weights
                    const float *weights = BiFilt->weights_filtering[xy].weight;


                    const int w_neighbor = index_interp.ei - index_interp.ii;
                    const int h_neighbor = index_interp.ej - index_interp.ij;

                    float numerator_u1 = 0.0;
                    float numerator_u2 = 0.0;
                    float denominator = 0.0;


                    for (int idx_j = 0; idx_j < h_neighbor; idx_j++) {
                        for (int idx_i = 0; idx_i < w_neighbor; idx_i++) {


                            const int idx_x = idx_i + index_interp.ii;
                            const int idx_y = idx_j + index_interp.ij;
                            const int idx_xy = idx_y * w + idx_x;
                            const int idx_ij = idx_j * w_neighbor + idx_i;

                            numerator_u1 += u1_filter[idx_xy] * weights[idx_ij];
                            numerator_u2 += u2_filter[idx_xy] * weights[idx_ij];
                            denominator += weights[idx_ij];

                        }
                    }

                    const float new_flow_u1 = numerator_u1 / denominator;
                    const float new_flow_u2 = numerator_u2 / denominator;
                    u1_filter[i] = new_flow_u1;
                    u2_filter[i] = new_flow_u2;
                }
            }
        }
    }

    // Save filtering in flow variable
    for (int j = 0; j < h_patch; j++) {
        for (int i = 0; i < w_patch; i++) {
            int x = i + patch.ii;
            int y = j + patch.ij;
            int xy = y * w + x;

            if (trust_points[xy] == 0 && fixed_points[xy] == 0) {
                u1[xy] = u1_filter[xy];
                u2[xy] = u2_filter[xy];
            }
        }
    }
}

// Insert n_neigh-connected candidates into the priority queue with their energies.
void insert_candidates(
        pq_cand &queue,
        float *ene_val,
        OpticalFlowData *ofD,
        const int i,
        const int j,
        const float ener_N,
        const int w,
        const int h) {

    int n_neigh = 4;
    int neighborhood[8][2] = {
            {0,  1},
            {0,  -1},
            {1,  0},
            {-1, 0},
            {1,  1},
            {1,  -1},
            {-1, 1},
            {-1, -1}};
    // w, h as params in the function call

    //const int w = ofD->params.w;
    //const int h = ofD->params.h;
    const float *sal = ofD->saliency;

    for (int k = 0; k < n_neigh; k++) {
        int px = i + neighborhood[k][0];
        int py = j + neighborhood[k][1];

        if (px >= 0 && px < w && py >= 0 && py < h) {
            float new_ener = ener_N * sal[py * w + px];

            //printf("Ener_N: %f  Sim: %f \n", ener_N, ene_val[py*w + px]);
            if (!ofD->fixed_points[py * w + px] && new_ener < ene_val[py * w + px]) {

                ene_val[py * w + px] = ener_N;
                SparseOF element;
                element.i = px;  // column
                element.j = py;  // row
                element.u = ofD->u1[py * w + px];
                element.v = ofD->u2[py * w + px];
                element.sim_node = new_ener;
                if(ofD->params.val_method >= 8) {
                    element.occluded = ofD->chi[py * w + px];
                }
                queue.push(element);
            }
        }
    }
}


// These types of comments should be removed before publishing!
// TODO: Esto esta fatal. Si finalmenente funciona lo de los pesos arreglarlo para
// que faldoi sea independiente y este dentro de energy_model.cpp
inline void get_relative_index_weight(
        int *iiw,  // initial column
        int *ijw,  // initial row
        const int wr,
        const int i,
        const int j
) {

    (*iiw) = (((i - wr) < 0) ? -(i - wr) : 0);
    (*ijw) = (((j - wr) < 0) ? -(j - wr) : 0);
    assert(*iiw >= 0);
    assert(*ijw >= 0);
}

static void get_index_weight(
        int method,
        SpecificOFStuff *ofS,
        const int wr,
        int i,
        int j) {

    int iiw, ijw;
    if (method == M_TVL1_W || method == M_NLTVCSAD_W || method == M_NLTVL1_W || method == M_TVCSAD_W) {
        get_relative_index_weight(&iiw, &ijw, wr, i, j);
    }
    switch (method) {
        case M_TVL1_W:
            ofS->tvl2w.iiw = iiw;
            ofS->tvl2w.ijw = ijw;
            break;
        case M_NLTVCSAD_W:
            ofS->nltvcsadw.iiw = iiw;
            ofS->nltvcsadw.ijw = ijw;
            break;
        case M_NLTVL1_W:
            ofS->nltvl1w.iiw = iiw;
            ofS->nltvl1w.ijw = ijw;
            break;
        case M_TVCSAD_W:
            ofS->tvcsadw.iiw = iiw;
            ofS->tvcsadw.ijw = ijw;
            break;
        default:
            break;
    }
}


// Copy over ofD->u1 and ofD->u2 the presented values in out.
inline void copy_fixed_coordinates(
        OpticalFlowData *ofD,
        const float *out,
        const PatchIndexes &index,
        const int w,
        const int h
) {
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    //Added changes for subimages
    //int w = ofD->params.w;
    //int h = ofD->params.h;
    int *fixed = ofD->fixed_points;

    for (int l = index.ij; l < index.ej; l++) {
        for (int k = index.ii; k < index.ei; k++) {
            // Copy only fixed values from the patch
            const int i = l * w + k;
            if (fixed[i] == 1) {
                u1[i] = out[i];
                u2[i] = out[w * h + i];
                assert(isfinite(u1[i]));
                assert(isfinite(u2[i]));
            }
        }
    }
}


// Check if there is at least one pixel that hasn't survived to the prunning.
int check_trustable_patch(
        OpticalFlowData *ofD,
        const PatchIndexes &index,
        const int w) {
    // w, h as params in the function call

    //const int w = ofD->params.w;

    int *fixed = ofD->trust_points;

    for (int l = index.ij; l < index.ej; l++)
        for (int k = index.ii; k < index.ei; k++) {
            // Return 0 if it detects that at least one point it is not fixed
            const int i = l * w + k;
            // If the pixel it is not trustable.
            if (fixed[i] == 0) {
                return 0;
            }
        }
    return 1;
}


static void add_neighbors(
        const float *i0,
        const float *i1,
        const float *i_1,
        float *ene_val,
        OpticalFlowData *ofD,
        SpecificOFStuff *ofS,
        pq_cand *queue,
        const int i,
        const int j,
        const int iteration,
        float *out,
        float *out_occ,
        BilateralFilterData *BiFilt,
        const int w,
        const int h
) {
    // Added w, h in as function params

    //const int w = ofD->params.w;
    //const int h = ofD->params.h;
    const int wr = ofD->params.w_radio;
    float ener_N;


    const PatchIndexes index = get_index_patch(wr, w, h, i, j, 1);
    int method = ofD->params.val_method; // used to include no occ.

    // These types of comments should be removed as well
    // TODO: Fix weights's stuff
    get_index_weight(method, ofS, wr, i, j);

    // In first iteration, Poisson interpolation
    if (iteration == 0) {
        // Interpolate by poisson on initialization
        copy_fixed_coordinates(ofD, out, index, w, h);
        // Poisson Interpolation (4wr x 4wr + 1)
        interpolate_poisson(ofD, index);

    } else {
        // Interpolate by bilateral filtering if some points do not survive to prunning
        if (check_trustable_patch(ofD, index, w) == 0) {

            copy_fixed_coordinates(ofD, out, index, w, h);
            bilateral_filter(ofD, BiFilt, index, w, h);

        }
    }

    // Optical flow method on patch (2*wr x 2wr + 1)
    of_estimation(ofS, ofD, &ener_N, i0, i1, i_1, index, w, h);

    // Insert new candidates to the queue
    insert_candidates(*queue, ene_val, ofD, i, j, ener_N, w, h);

    // It is a strange step, if the energy over the patch is lower thant the
    // stored energy, we put the new one, if it's not, we leave the old one.
    if (ene_val[j * w + i] > ener_N) {
        out[j * w + i] = ofD->u1[j * w + i];
        out[w * h + j * w + i] = ofD->u2[j * w + i];
        ene_val[j * w + i] = ener_N;
        // Only if 'occlusions'
        if (method >= 8) {
            out_occ[j * w + i] = ofD->chi[j * w + i];
        }
    }
}

void insert_initial_seeds(
        const float *i0,
        const float *i1,
        const float *i_1,
        float *in_flow,
        pq_cand *queue,
        OpticalFlowData *ofD,
        SpecificOFStuff *ofS,
        float *ene_val,
        float *out_flow,
        float *out_occ,
        BilateralFilterData *BiFilt,
        const int w,
        const int h) {
    // w, h as params in the function call

    //const int w = ofD->params.w;
    //const int h = ofD->params.h;
    const int wr = ofD->params.w_radio;


    // Set to the initial conditions all the stuff
    for (int i = 0; i < w * h; i++) {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out_flow[i] = NAN;
        out_flow[w * h + i] = NAN;
        out_occ[i] = 0;
        ofD->trust_points[i] = 1;
    }

    ofD->params.w_radio = 1;
    // Fix the initial seeds.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {

            // Indicates the initial seed in the similarity map
            if (isfinite(in_flow[j * w + i]) && isfinite(in_flow[w * h + j * w + i])) {

                out_flow[j * w + i] = in_flow[j * w + i];
                out_flow[w * h + j * w + i] = in_flow[w * h + j * w + i];
                ofD->fixed_points[j * w + i] = 1;

                // Add_neigbors 0 means that during the propagation interpolates the patch
                // based on the energy.
                add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, 0, out_flow, out_occ, BiFilt, w, h);

                // These values may have been modified in the previous function
                out_flow[j * w + i] = in_flow[j * w + i];
                out_flow[w * h + j * w + i] = in_flow[w * h + j * w + i];
                ofD->fixed_points[j * w + i] = 1;
                ene_val[j * w + i] = 0.0;
            }
        }
    ofD->params.w_radio = wr;
}


// Insert each pixel into the queue as possible candidate. Its related energy comes
// from the energy store at the moment that the pixel was fixed.
void insert_potential_candidates(
        const float *in,
        OpticalFlowData *ofD,
        pq_cand &queue,
        const float *ene_val,
        const float *out_occ,
        const int w,
        const int h
) {
    // Added changes for subimages

    // Note: in and out are the same pointer
    //const int w = ofD->params.w;
    //const int h = ofD->params.h;
    // Fixed the initial seeds.
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            // Indicates the initial seed in the similarity map
            if (isfinite(in[j * w + i]) && isfinite(in[w * h + j * w + i])) {

                SparseOF element;
                element.i = i;      // column
                element.j = j;      // row
                element.u = in[j * w + i];
                element.v = in[w * h + j * w + i];
                // Obs: Notice that ene_val contains (en)*saliency
                element.sim_node = ene_val[j * w + i];
                if (ofD->params.val_method >= 8) {
                    element.occluded = out_occ[j * w + i];
                }
                assert(isfinite(ene_val[j * w + i]));
                queue.push(element);
            }
        }
    }
}


// Initialize the data to prepare everything for the region growing
void prepare_data_for_growing(
        OpticalFlowData *ofD,
        float *ene_val,
        float *out,
        const int w,
        const int h) {
    // Added w, h as params in the function call

    //int w = ofD->params.w;
    //int h = ofD->params.h;

    // Set to the initial conditions all the stuff
    for (int i = 0; i < w * h; i++) {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out[i] = NAN;
        out[w * h + i] = NAN;
    }
}


void local_growing(
        const float *i0,
        const float *i1,
        const float *i_1,
        pq_cand *queue,
        SpecificOFStuff *ofS,
        OpticalFlowData *ofD,
        int iteration,
        float *ene_val,
        float *out_flow,
        float *out_occ,
        BilateralFilterData *BiFilt,
        bool fwd_or_bwd,
        const int w,
        const int h
) {
    // Added w, h in as parameters

    int fixed = 0;
    vector<int> percent_print = {30, 70, 80, 95, 100};
    //const int w = ofD->params.w;
    //const int h = ofD->params.h;
    const int size = w * h;
    printf("Queue size at start = %d\n", (int) queue->size());
    while (!queue->empty()) {

        SparseOF element = queue->top();
        int i = element.i;
        int j = element.j;
        // While the queue is not empty, take an element to process
        queue->pop();

        if (!ofD->fixed_points[j * w + i]) {
            fixed++;
            assert(isfinite(element.sim_node));
            float u = element.u;
            float v = element.v;
            float energy = element.sim_node;
            float occlusion;

            if (ofD->params.val_method >= 8) {
                occlusion = element.occluded;
            } else {
                occlusion = 0.0;
            }

            if (!isfinite(u)) {
                printf("U1 = %f\n", u);
            }
            if (!isfinite(v)) {
                printf("U2 = %f\n", v);
            }

            ofD->fixed_points[j * w + i] = 1;

            out_flow[j * w + i] = u;
            out_flow[w * h + j * w + i] = v;
            ene_val[j * w + i] = energy;

            out_occ[j * w + i] = occlusion;

            // TODO: copied so those values influence the minimization
            // ofD->u1[j*w + i] = u;
            // ofD->u2[j*w + i] = v;

            add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, iteration, out_flow, out_occ, BiFilt, w, h);

            float percent = 100 * fixed * 1.0 / size * 1.0;

            if (SAVE_RESULTS == 1) {
                for (int k = 0; k < 4; k++) {
                    if (percent > percent_print[k] && percent < percent_print[k + 1]) {
                        string filename_flow = " ";
                        string filename_occ = " ";
                        if (fwd_or_bwd) {
                            filename_flow =
                                    "../Results/Partial_results/partial_results_fwd_" + to_string(percent_print[k]) +
                                    "_iter_" + to_string(iteration) + ".flo";
                            iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);
                            filename_occ =
                                    "../Results/Partial_results/partial_results_fwd_" + to_string(percent_print[k]) +
                                    "_iter_" + to_string(iteration) + "_occ.png";
                            auto *out_occ_int = new int[w * h];

                            for (int i = 0; i < w * h; i++) {

                                out_occ_int[i] = out_occ[i];
                            }

                            iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);

                            percent_print[k] = 200;
                        }
                    }
                }
            }
        }
    }
    if (SAVE_RESULTS == 1) {
        if (fwd_or_bwd) {
            string filename_flow =
                    "../Results/Partial_results/partial_results_fwd_100_iter_" + to_string(iteration) + ".flo";
            iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);
            string filename_occ =
                    "../Results/Partial_results/partial_results_fwd_100_iter_" + to_string(iteration) + "_occ.png";
            auto *out_occ_int = new int[w * h];

            for (int i = 0; i < w * h; i++) {

                out_occ_int[i] = out_occ[i];
            }

            iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);

        }
    }
}


void match_growing_variational(
        float *go,
        float *ba,
        float *i0,
        float *i1,
        float *i_1,
        float *i2,
        float *sal_go,
        float *sal_ba,
        Parameters params,
        float *ene_val,
        float *out_flow,
        float *out_occ
) {
    using namespace chrono;  // debug
    auto clk1 = system_clock::now();
    int w = params.w;
    int h = params.h;

    printf("Initializing stuff\n");
    // Initialize all the stuff for optical flow computation
    // Optical flow t, t+1
    OpticalFlowData ofGo = init_Optical_Flow_Data(sal_go, params, w, h);
    auto *oft0 = new float[w * h * 2];
    auto *ene_Go = new float[w * h];
    auto *occ_Go = new float[w * h];


    // Optical flow t+1, t
    OpticalFlowData ofBa = init_Optical_Flow_Data(sal_ba, params, w, h);
    auto *oft1 = new float[w * h * 2];
    auto *ene_Ba = new float[w * h];
    auto *occ_Ba = new float[w * h];

    // Create queues
    // Partitioning into 6 'subimages', we will need 12 queues, 6 fwd and 6 bwd
    // Nomenclature: queueGo_1, ..., queueGo_6 and queueBa_1, ..., queueBa_6

    pq_cand queue_Go;
    pq_cand queue_Ba;

    // Initialize all the auxiliar data.
    SpecificOFStuff stuffGo;
    SpecificOFStuff stuffBa;
    initialize_auxiliar_stuff(stuffGo, ofGo, w, h);
    initialize_auxiliar_stuff(stuffBa, ofBa, w, h);

    // i0n, i1n, i_1n, i2n are a gray and smooth version of i0, i1, i_1, i2
    float *i0n = nullptr;
    float *i1n = nullptr;
    float *i_1n = nullptr;
    float *i2n = nullptr;
    // Prepare data based on the functional chosen (energy_model.cpp)
    prepare_stuff(&stuffGo, &ofGo, &stuffBa, &ofBa, i0, i1, i_1, i2, params.pd, &i0n, &i1n, &i_1n, &i2n, w, h);

    // Initialize weights for bilateral filtering
    auto BiFilt_Go = init_weights_bilateral(i0n, w, h);
    auto BiFilt_Ba = init_weights_bilateral(i1n, w, h);

    printf("Finished initializing stuff\n");


    auto clk2 = system_clock::now(); // DEBUG
    duration<double> elapsed_secs2 = clk2 - clk1; // DEBUG
    cout << "(match growing) Initializing everything took "
         << elapsed_secs2.count() << endl;

    // Insert initial seeds to queues
    printf("Inserting initial seeds\n");

    // We left occ_Go and occ_Ba as initialized to avoid further changes (just do not use
    // them)
    auto future_nfixed_go = async(launch::async,
                                  [&] {
                                      return insert_initial_seeds(i0n, i1n, i_1n, go, &queue_Go, &ofGo, &stuffGo, ene_Go,
                                                                  oft0, occ_Go, BiFilt_Go, w, h);
                                  });
	
    insert_initial_seeds(i1n, i0n, i2n, ba, &queue_Ba, &ofBa, &stuffBa, ene_Ba,
                         oft1, occ_Ba, BiFilt_Ba, w, h);
    future_nfixed_go.get();


// OpenMP (tests to optimise performance at HPC)
/*
#pragma omp parallel
{
    #pragma omp single
    {
        #pragma omp task
        {
            printf("(initial_seeds-forward) Hello from thread num %d!\n", omp_get_thread_num());
	    nfixed = insert_initial_seeds(i0n, i1n, i_1n, go, &queueGo, &ofGo, &stuffGo, ene_Go, oft0, occ_Go, BiFilt_Go);
        }
	#pragma omp task
	{
	    printf("(initial_seeds-backward) Hello from thread num %d!\n", omp_get_thread_num());
            nfixed = insert_initial_seeds(i1n, i0n, i2n, ba, &queueBa, &ofBa, &stuffBa, ene_Ba, oft1, occ_Ba, BiFilt_Ba);
	} 
    }
}
*/

/*    
    omp_set_dynamic(0);
    omp_set_nested(1);
   
#pragma omp parallel sections
{
    #pragma omp section
    {
        printf("Inserting initial forward seeds...\n");
	insert_initial_seeds(i0n, i1n, i_1n, go, &queueGo, &ofGo, &stuffGo, ene_Go, oft0, occ_Go, BiFilt_Go);
    }
    
    #pragma omp section
    {
        printf("Inserting initial backward seeds...\n");
	insert_initial_seeds(i1n, i0n, i2n, ba, &queueBa, &ofBa, &stuffBa, ene_Ba, oft1, occ_Ba, BiFilt_Ba); 
    }
}   
*/
    printf("Finished inserting initial seeds\n");

    
    /*auto clk3 = system_clock::now(); // DEBUG
    duration<double> elapsed_secs3 = clk3 - clk2; // DEBUG
    cout << "(match growing) inserting initial seeds took "
         << elapsed_secs3.count() << endl;
    */
    const int iter = params.iterations_of;  //LOCAL_ITER;
    // Variables for pruning
    float tol[2] = {FB_TOL, TU_TOL};
    int p[2] = {1, 0};

    // Initialise the partitions
    // Second and third iteration (i == 1, 2) ==> h_parts x v_parts, v_parts x h_parts
//
//    int *sub_w = nullptr, *sub_w_r = nullptr;           // Array with the width of each partition
//    int *sub_h = nullptr, *sub_h_r = nullptr;           //   "    "    "  height "  "      "
//    int *off_x = nullptr, *off_x_r = nullptr;           //   "    "    "  width offset   "   "   "
//    int *off_y = nullptr, *off_y_r = nullptr;           //   "    "    "  height   "     "   "   "
//
//    // Queues (two per partition: 1 fwd + 1 bwd)
//    std::vector<pq_cand> queue_Go_p, queue_Go_p_r;
//    std::vector<pq_cand> queue_Ba_p, queue_Ba_p_r;
//
//    // Specific stuff (functional-specific)
//    std::vector<SpecificOFStuff> stuffGo_sub, stuffGo_sub_r;
//    std::vector<SpecificOFStuff> stuffBa_sub, stuffBa_sub_r;
//
//    // Optical Flow Data (common parameters for OF estimation)
//    std::vector<OpticalFlowData> ofGo_sub, ofGo_sub_r;
//    std::vector<OpticalFlowData> ofBa_sub, ofBa_sub_r;

    // Initialise partitions
    std::vector<PartitionData*> p_data;
    std::vector<PartitionData*> p_data_r;
    // Second iteration (i == 1)
    init_subimage_partitions(i0, i1, i_1, i2, sal_go, sal_ba, w, h, params.h_parts, params.v_parts, &p_data, params);

    // Third iteration (i == 2) ==> v_parts x h_parts
    init_subimage_partitions(i0, i1, i_1, i2, sal_go, sal_ba, w, h, params.v_parts, params.h_parts, &p_data_r, params);

    // TODO: remove all the redundant code to clean everything
    // Frames
    float *i0_sub = nullptr, *i0_sub_r = nullptr;
    float *i1_sub = nullptr, *i1_sub_r = nullptr;
    float *i_1_sub = nullptr, *i_1_sub_r = nullptr;
    float *i2_sub = nullptr, *i2_sub_r = nullptr;

    // OF, energy, occlusions
    float *oft0_sub = nullptr, *oft0_sub_r = nullptr;
    float *oft1_sub = nullptr, *oft1_sub_r = nullptr;
    float *ene_Go_sub = nullptr, *ene_Go_sub_r = nullptr;
    float *ene_Ba_sub = nullptr, *ene_Ba_sub_r = nullptr;
    float *occ_Go_sub = nullptr, *occ_Go_sub_r = nullptr;
    float *occ_Ba_sub = nullptr, *occ_Ba_sub_r = nullptr;


    // Optical Flow (fwd and bwd)
    std::vector<float *> oft0_p, oft0_p_r;
    std::vector<float *> oft1_p, oft1_p_r;

    // Energy (fwd and bwd)
    std::vector<float *> ene_Go_p, ene_Go_p_r;
    std::vector<float *> ene_Ba_p, ene_Ba_p_r;

    // Occlusions (fwd and bwd)
    std::vector<float *> occ_Go_p, occ_Go_p_r;
    std::vector<float *> occ_Ba_p, occ_Ba_p_r;

    // Pointers to frames and their normalized versions (two different partitions)
    std::vector<float *> i0_p, i1_p, i_1_p, i2_p, i0_p_r, i1_p_r, i_1_p_r, i2_p_r;
    std::vector<float *> i0n_p, i1n_p, i_1n_p, i2n_p, i0n_p_r, i1n_p_r, i_1n_p_r, i2n_p_r;

    // Bilateral filters computed for each partition (to have correct indexes and weights)
    std::vector<BilateralFilterData*> BiFilt_Go_p, BiFilt_Go_p_r;
    std::vector<BilateralFilterData*> BiFilt_Ba_p, BiFilt_Ba_p_r;

    // OpticalFlowData
    std::vector<OpticalFlowData> ofGo_p, ofGo_p_r;
    std::vector<OpticalFlowData> ofBa_p, ofBa_p_r;

    // SpecificOFStuff
    std::vector<SpecificOFStuff> stuffGo_p, stuffGo_p_r;
    std::vector<SpecificOFStuff> stuffBa_p, stuffBa_p_r;

    if (params.split_img == 1)
    {
        // Fill partitions from original frames
//        i0_sub = fill_subimage_partition(i0, w, params.pd, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//        i0_sub_r= fill_subimage_partition(i0, w, params.pd, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//
//        i1_sub = fill_subimage_partition(i1, w, params.pd, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//        i1_sub_r = fill_subimage_partition(i1, w, params.pd, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//
//        i_1_sub = fill_subimage_partition(i_1, w, params.pd, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//        i_1_sub_r = fill_subimage_partition(i_1, w, params.pd, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//
//        i2_sub = fill_subimage_partition(i2, w, params.pd, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//        i2_sub_r = fill_subimage_partition(i2, w, params.pd, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);

        // Initialise OpticalFlowData and SpecificOFStuff
        for (int prt = 0; prt < params.h_parts * params.v_parts; prt++)
        {

        }

    }

    for (int i = 0; i < iter; i++) {
       // auto clk4 = system_clock::now();  // DEBUG
        printf("Iteration: %d\n", i);

        // Estimate local minimization (I0-I1)
        // First iteration work on the whole image
        if (params.split_img == 0 || (params.split_img == 1 && i == 0))
        {
            auto growing_fwd = async(launch::async,
                                     [&] {
                                         local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, i, ene_Go, oft0,
                                                       occ_Go, BiFilt_Go, true, w, h);
                                     });

            // auto clk5 = system_clock::now(); // DEBUG
            //duration<double> elapsed_secs4 = clk5- clk4; // DEBUG
            //cout << "(match growing) ASYNCH NOT TRUELocal iteration " << i << " => local_growing (I0-I1) took "
            //   << elapsed_secs4.count() << endl;


            // Estimate local minimization (I1-I0)
            local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba, false, w, h);
            //auto clk5 = system_clock::now(); // DEBUG
            //duration<double> elapsed_secs5 = clk5- clk4; // DEBUG
            //cout << "(match growing) Local iteration " << i << " => local_growing (I1-I0) took "
            //     << elapsed_secs5.count() << endl;

            growing_fwd.get();  // HERE the first growing is retrieved
            //auto clk_extra = system_clock::now();
            //duration<double> elapsed_extra = clk_extra- clk4; // DEBUG
            //cout << "(match growing) ASYNC Local iteration " << i << " => local_growing (I0-I1) took "
            //     << elapsed_secs5.count() << endl;


            /*
                #pragma omp parallel
                {
                    #pragma omp single
                    {
                        #pragma omp task
                    {
                            printf("(local_growing-forward) Hello from thread num %d!\n", omp_get_thread_num());
                    local_growing(i0n, i1n, i_1n, &queueGo, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go,
                                                                    BiFilt_Go, true);
                    }
                    #pragma omp task
                    {
                    printf("(local_growing-backward) Hello from thread num %d!\n", omp_get_thread_num());
                        local_growing(i1n, i0n, i2n, &queueBa, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba,
                                        BiFilt_Ba, false);
                    }
                    }
                }
            */

            /*
            #pragma omp parallel sections
            {
            #pragma omp section
            {
                printf("Starting forward local growing...\n");
                local_growing(i0n, i1n, i_1n, &queueGo, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go,
                                                           BiFilt_Go, true);
            }

            #pragma omp section
            {
                printf("Starting backward local growing...\n");
                local_growing(i1n, i0n, i2n, &queueBa, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba,
                                   BiFilt_Ba, false);
            }
            }
            */
            p_data[0]->ofGo.u1[0] = ofGo.u1[0];
            p_data[0]->stuffGo.tvl2.xi11[0] = stuffGo.tvl2.xi11[0];

            // Pruning method
            pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

            /*auto clk7 = system_clock::now(); // DEBUG
            duration<double> elapsed_secs6 = clk7- clk5; // DEBUG
            cout << "(match growing) Local iteration " << i << " => pruning method took "
                 << elapsed_secs6.count() << endl;

        */
            // Delete not trustable candidates based on the previous pruning
            delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
            delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

            /*
                auto clk8 = system_clock::now(); // DEBUG
                duration<double> elapsed_secs7 = clk8- clk7; // DEBUG
                cout << "(match growing) Local iteration " << i << " => delete non-trustable candidates "
                     << elapsed_secs7.count() << endl;

                */
            // Insert each pixel into the queue as possible candidate
            insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, occ_Go, w, h);
            insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, occ_Ba, w, h);


            /*auto clk9 = system_clock::now(); // DEBUG
                duration<double> elapsed_secs8 = clk9- clk8; // DEBUG
                cout << "(match growing) Local iteration " << i << " => insert potential candidates "
                     << elapsed_secs8.count() << endl;
            */
            prepare_data_for_growing(&ofGo, ene_Go, oft0, w, h);
            prepare_data_for_growing(&ofBa, ene_Ba, oft1, w, h);

            /*auto clk10 = system_clock::now(); // DEBUG
            duration<double> elapsed_secs9 = clk10- clk9; // DEBUG
            cout << "(match growing) Local iteration " << i << " => prepare data for growing "
                 << elapsed_secs9.count() << endl;

            auto clk11 = system_clock::now(); // DEBUG
            duration<double> elapsed_secs10 = clk11- clk4; // DEBUG
            cout << "(match growing) Local iteration " << i << " => all iteration's tasks took "
                 << elapsed_secs10.count() << endl;*/
        }
        else if ((i > 0 && i <= iter - 1) && params.split_img == 1)
        {
            // Common stuff to any iteration from 2nd to last - 1

            if (i % 2 != 0 && i <= iter - 1)           // part. grid: h_parts x v_parts
            {
                // Get HERE oft0, oft1, eneGo, eneBa, ... UPDATED values
                // Fill with last iteration's OF
//                oft0_sub = fill_subimage_partition(oft0, w, 2, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                oft0_sub_r = fill_subimage_partition(oft0, w, 2, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//                oft1_sub = fill_subimage_partition(oft1, w, 2, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                oft1_sub_r = fill_subimage_partition(oft1, w, 2, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//
//                // Fill with last iteration's energy
//                ene_Go_sub = fill_subimage_partition(ene_Go, w, 1, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                ene_Go_sub_r = fill_subimage_partition(ene_Go, w, 1, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//                ene_Ba_sub = fill_subimage_partition(ene_Ba, w, 1, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                ene_Ba_sub_r = fill_subimage_partition(ene_Ba, w, 1, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//
//                // Fill with last iteration's occlusions
//                occ_Go_sub = fill_subimage_partition(occ_Go, w, 1, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                occ_Go_sub_r = fill_subimage_partition(occ_Go, w, 1, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//                occ_Ba_sub = fill_subimage_partition(occ_Ba, w, 1, params.h_parts, params.v_parts, sub_w, sub_h, off_x, off_y);
//                occ_Ba_sub_r = fill_subimage_partition(occ_Ba, w, 1, params.v_parts, params.h_parts, sub_w_r, sub_h_r, off_x_r, off_y_r);
//

                //#pragma omp parallel for
                for (int prt = 0; prt < params.h_parts * params.v_parts; prt++)
                {
                    // Temporal pointers to correct part of partition array (created above)
                    // Frames
                    /*
                    float *i0n_p = get_subimage_partition(i0n_sub, prt, sub_w, sub_h, off_x, off_y, 1);
                    float *i1n_p = get_subimage_partition(i1n_sub, prt, sub_w, sub_h, off_x, off_y, 1);
                    float *i_1n_p = get_subimage_partition(i_1n_sub, prt, sub_w, sub_h, off_x, off_y, 1);
                    float *i2n_p = get_subimage_partition(i2n_sub, prt, sub_w, sub_h, off_x, off_y, 1);
                    */

                    // Fwd and bwd OF
//                    get_subimage_partition(oft0_sub, &oft0_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//                    get_subimage_partition(oft1_sub, &oft1_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//
//                    // Fwd and bwd energy
//                    get_subimage_partition(ene_Go_sub, &ene_Go_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//                    get_subimage_partition(ene_Ba_sub, &ene_Ba_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//
//                    // Fwd and bwd occlusions
//                    get_subimage_partition(occ_Go_sub, &occ_Go_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//                    get_subimage_partition(occ_Ba_sub, &occ_Ba_p[prt], prt, sub_w, sub_h, 1, off_x, off_y);
//
//
//                    // 1. Local growing based on updated oft0 and oft1 of the first iteration
//                    // FWD
//                    local_growing(i0n_p[prt], i1n_p[prt], i_1n_p[prt], &queue_Go_p[prt], &stuffGo_p[prt],
//                                  &ofGo_p[prt], i, ene_Go_p[prt], oft0_p[prt], occ_Go_p[prt], BiFilt_Go_p[prt],
//                                  true, sub_w[prt], sub_h[prt]);
//                    // BWD
//                    local_growing(i1n_p[prt], i0n_p[prt], i2n_p[prt], &queue_Ba_p[prt], &stuffBa_sub[prt],
//                                  &ofBa_sub[prt], i, ene_Ba_p[prt], oft1_p[prt],
//                                  occ_Ba_p[prt], BiFilt_Ba_p[prt], false, sub_w[prt], sub_h[prt]);

                }

                    //TODO: call here function to get back 'whole-image' info'
                    // NOTE: these two functions (one per 'pointers to float' and the other for the struct)
                    // are the 'inverse' functions of those used in the next TODO stament

                    // 2. Pruning
                    pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

                    // 3. Delete non-trustable
                    delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
                    delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

                // TODO: go back to partitions
                // ('get_subimage...' for oft0, oft1, ene_Go and ene_Ba
                //  ofGO, ofBa...through a different 'new' function that maps indices

                for (int prt = 0; prt < params.h_parts * params.v_parts; prt++)
                {
                    // 4. insert potential candidates
//                    insert_potential_candidates(oft0_p[prt], &ofGo_p[prt], queue_Go_p[prt], ene_Go_p[prt],
//                                                occ_Go_p[prt], sub_w[prt], sub_h[prt]);
//                    insert_potential_candidates(oft1_p[prt], &ofBa_p[prt], queue_Ba_p[prt], ene_Ba_p[prt],
//                                                occ_Ba_p[prt], sub_w[prt], sub_h[prt]);
//
//                    // 5. prepare data for growing
//                    prepare_data_for_growing(&ofGo_p[prt], ene_Go_p[prt], oft0_p[prt], sub_w[prt], sub_h[prt]);
//                    prepare_data_for_growing(&ofBa_p[prt], ene_Ba_p[prt], oft1_p[prt], sub_w[prt], sub_h[prt]);
                }

            }
            else if (i % 2 == 0 && i <= iter - 1)        // part. grid: v_parts x h_parts (previous one 'reversed' (transposed))
            {
                //#pragma omp parallel for
                for (int prt = 0; prt < params.v_parts * params.h_parts; prt++) {
                    // Temporal pointers to correct part of partition array (created above)
                    // Frames
                    /*float *i0n_p_r = get_subimage_partition(i0n_sub_r, prt, sub_w_r, sub_h_r, off_x_r, off_y_r, 1);
                    float *i1n_p_r = get_subimage_partition(i1n_sub_r, prt, sub_w_r, sub_h_r, off_x_r, off_y_r, 1);
                    float *i_1n_p_r = get_subimage_partition(i_1n_sub_r, prt, sub_w_r, sub_h_r, off_x_r, off_y_r, 1);
                    float *i2n_p_r = get_subimage_partition(i2n_sub_r, prt, sub_w_r, sub_h_r, off_x_r, off_y_r, 1);
                     */

                    // Fwd and bwd OF
//                    get_subimage_partition(oft0_sub_r, &oft0_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//                    get_subimage_partition(oft1_sub_r, &oft1_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//
//                    // Fwd and bwd energy
//                    get_subimage_partition(ene_Go_sub_r, &ene_Go_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//                    get_subimage_partition(ene_Ba_sub_r, &ene_Ba_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//
//                    // Fwd and bwd occlusions
//                    get_subimage_partition(occ_Go_sub_r, &occ_Go_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//                    get_subimage_partition(occ_Ba_sub_r, &occ_Ba_p_r[prt], prt, sub_w_r, sub_h_r, 1, off_x_r, off_y_r);
//
//
//                    // 1. Local growing based on updated oft0 and oft1 of the first iteration
//                    // FWD
//                    local_growing(i0n_p_r[prt], i1n_p_r[prt], i_1n_p_r[prt], &queue_Go_p_r[prt], &stuffGo_p_r[prt],
//                                  &ofGo_p_r[prt], i, ene_Go_p_r[prt], oft0_p_r[prt], occ_Go_p_r[prt],
//                                  BiFilt_Go_p_r[prt],
//                                  true, sub_w_r[prt], sub_h_r[prt]);
//                    // BWD
//                    local_growing(i1n_p_r[prt], i0n_p_r[prt], i2n_p_r[prt], &queue_Ba_p_r[prt], &stuffBa_p_r[prt],
//                                  &ofBa_p_r[prt], i, ene_Ba_p_r[prt], oft1_p_r[prt], occ_Ba_p_r[prt],
//                                  BiFilt_Ba_p_r[prt],
//                                  false, sub_w_r[prt], sub_h_r[prt]);
                }

                    //TODO: call here function to get back 'whole-image' info'
                    // The next iteration will use the whole image so we do not need to create again the
                    // partitions

                    // 2. Pruning
                    pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

                    // 3. Delete non-trustable
                    delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
                    delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

                    // 4. insert potential candidates
                    insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, occ_Go, w, h);
                    insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, occ_Ba, w, h);

                    // 5. prepare data for growing
                    prepare_data_for_growing(&ofGo, ene_Go, oft0, w, h);
                    prepare_data_for_growing(&ofBa, ene_Ba, oft1, w, h);

            }

        }
    }
    //auto last_growing = system_clock::now();    //DEBUG


    printf("Last growing\n");
    local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, iter, ene_Go, oft0, occ_Go, BiFilt_Go, true, w, h);


    /*auto clk_end = system_clock::now(); // DEBUG
    duration<double> elapsed_secs = clk_end- last_growing; // DEBUG
    cout << "(match growing) Last growing took "
         << elapsed_secs.count() << endl;
    */	
    // Copy the result t, t+1 as output.
    memcpy(out_flow, oft0, sizeof(float) * w * h * 2);
    memcpy(ene_val, ene_Go, sizeof(float) * w * h);
    memcpy(out_occ, occ_Go, sizeof(float) * w * h);


    free_auxiliar_stuff(&stuffGo, &ofGo);
    free_auxiliar_stuff(&stuffBa, &ofBa);


    delete[] i1n;
    delete[] i2n;
    delete[] i0n;
    delete[] i_1n;

    delete[] ofGo.u1;
    delete[] ofBa.u1;

    delete[] ofGo.fixed_points;
    delete[] ofBa.fixed_points;

    delete[] ofGo.trust_points;
    delete[] ofBa.trust_points;

    delete[] oft0;
    delete[] oft1;

    delete[] ene_Go;
    delete[] ene_Ba;

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////MAIN/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Main function that expands sparse flow
int main(int argc, char *argv[]) {


    using namespace chrono;

    system_clock::time_point today = system_clock::now();
    time_t tt;

    tt = system_clock::to_time_t(today);
    cerr << "Starting  date: " << ctime(&tt);


    /*auto clk1 = system_clock::now(); // DEBUG*/
    // Process input optional parameters
    vector<string> args(argv, argv + argc);
    auto windows_ratio = pick_option(args, "wr",
                                     to_string(PAR_DEFAULT_WINSIZE));           // Windows ratio
    auto var_reg = pick_option(args, "m", to_string(M_TVL1));                   // Methods (default: tvl1)
    auto file_params = pick_option(args, "p", "");                              // File of parameters
    auto local_iters = pick_option(args, "loc_it",
                                   to_string(LOCAL_ITER));                      // Local Faldoi number of iterations
    auto max_iters_patch = pick_option(args, "max_pch_it",
                                       to_string(MAX_ITERATIONS_LOCAL));        // Iteration per patch of win_ratio
    auto split_img = pick_option(args, "split_img", to_string(PARTITIONING));   // Whether to split into subimages or not
    auto hor_parts = pick_option(args, "h_parts", to_string(HOR_PARTS));        // Number of horizontal slices (partition)
    auto ver_parts = pick_option(args, "v_parts", to_string(VER_PARTS));        // "       " vertical      "        "


    // [F]: CHANGE: this does not work for methods without occlusions!
    //if (args.size() != 7 && args.size() != 9) {
    if (args.size() < 6 || args.size() > 9) {
        // Without occlusions
        fprintf(stderr, "Without occlusions (nº of params: 5 or 7 + 1 (own function name)):\n");
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff"
                " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                " [-split_img split_image] [-h_parts horiz_parts]"
                " [-v_parts vert_parts] \n", args.size(), args[0].c_str());
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff sal0.tiff sal1.tiff"
                " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                " [-split_img split_image] [-h_parts horiz_parts]"
                " [-v_parts vert_parts] \n", args.size(), args[0].c_str());
        fprintf(stderr, "\n");
        // With occlusions
        fprintf(stderr, "With occlusions (nº of params: 7 or 9 + 1 (own function name)):\n");
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png"
                " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                " [-split_img split_image] [-h_parts horiz_parts]"
                " [-v_parts vert_parts] \n", args.size(), args[0].c_str());
        fprintf(stderr,
                "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png sal0.tiff sal1.tiff"
                " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                " [-split_img split_image] [-h_parts horiz_parts]"
                " [-v_parts vert_parts] \n", args.size(), args[0].c_str());
        return 1;
    }

    // filename that contains all the images to use
    string filename_i_1, filename_i0, filename_i1, filename_i2;

    // Read txt file of images
    const string &filename_images = args[1];
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
                } else {
                    if (num_files == 4) {
                        filename_i2 = line;
                    }
                }
            }
        }
    }
    infile.close();

    if (num_files == 3) {
        fprintf(stderr, "ERROR: 3 images given as input\n");
        // [F]
        // In the example the order of frames is I0 => I1 => I2 ==> I3 that correspond to
        // I-1, I0, I1 and I2 as they are consecutive in time.
        // But in fact, they are read as prompted (I0 => I1 ==> I_1 ==> I2) so the
        // example data file containing the path to the images is not properly formatted
        fprintf(stderr, "Without occlusions:\n");
        fprintf(stderr, "Usage: 2 images in the following order: I0, I1\n");
        fprintf(stderr, "With occlusions:\n");
        fprintf(stderr, "Usage: 4 images in the following order: I0, I1, I-1, I2\n");

        return 1;
    }

    // [F]: CHANGE: add the case of 6 and 8 parameters and check how to set occ param when no occ are wanted==>
    // probably need to change match_gro.... as well.
    // Save other arguments
    const string &filename_go = args[2];
    const string &filename_ba = args[3];
    const string &filename_out = args[4];
    const string &filename_sim = args[5];
    // Case w/o occ or saliency has all params already (args.size()==6)

    if (args.size() == 7 || args.size() == 9) // case with occlusions
    {
        const string &filename_occ = args[6];
    }
    const string &filename_occ = "";
    const char *filename_sal0 = nullptr;
    const char *filename_sal1 = nullptr;

    if (args.size() == 9) {
        filename_sal0 = args[7].c_str();
        filename_sal1 = args[8].c_str();
    } else if (args.size() == 8) {
        filename_sal0 = args[6].c_str();
        filename_sal1 = args[7].c_str();
    }

    // Optional arguments
    int w_radio = stoi(windows_ratio);
    int val_method = stoi(var_reg);
    int loc_it = stoi(local_iters);
    int max_it_pch = stoi(max_iters_patch);
    int sp_img = stoi(split_img);
    int h_prts = stoi(hor_parts);
    int v_prts = stoi(ver_parts);

    // Open input images and .flo
    // pd: number of channels
    int w[8], h[8], pd[6];

    // Convert relative path to absolute path (so it works with any input)
    // Added instruction to avoid failure when working with relative paths with
    // "~/folder/subfolder/..." (to be done in the future, some errors were reported
    // by iio.c (it fails for some reason)
    /*std::string filename_i_1_abs = path_abs2rel(filename_i_1);
    std::string filename_i2_abs = path_abs2rel(filename_i2);
    std::string filename_i0_abs = path_abs2rel(filename_i0);
    std::string filename_i1_abs = path_abs2rel(filename_i1);*/

    // Frame t-1 and t+2
    float *i_1 = nullptr;
    float *i2 = nullptr;
    if (num_files == 4) {
        //std::printf("Local faldoi line 1220\n");
//        std::cout << "I_1 filename (c_str): " << filename_i_1.c_str() << std::endl;
//        std::cout << "I2 filename (c_str): " << filename_i2.c_str() << std::endl;
        i_1 = iio_read_image_float_split(filename_i_1.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i2.c_str(), w + 7, h + 7, pd + 5);
    } else {
//        std::printf("This should NOT be printed\n");
        i_1 = iio_read_image_float_split(filename_i0.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i1.c_str(), w + 7, h + 7, pd + 5);
    }
//    std::printf("Local faldoi line 1227\n");
    // Frames t and t+1
    /*std::cout << "I0filename (c_str): " << filename_i0.c_str() << std::endl;
    std::cout << "I1 filename (c_str): " << filename_i1.c_str() << std::endl;*/
    float *i0 = iio_read_image_float_split(filename_i0.c_str(), w + 0, h + 0, pd + 0);
    float *i1 = iio_read_image_float_split(filename_i1.c_str(), w + 1, h + 1, pd + 1);

    // Sparse Optical flow forward and backward
    float *go = iio_read_image_float_split(filename_go.c_str(), w + 2, h + 2, pd + 2);
    float *ba = iio_read_image_float_split(filename_ba.c_str(), w + 3, h + 3, pd + 3);

    // Ensure dimensions match in images
    if (num_files == 4) {
        if (w[0] != w[1] || h[0] != h[1] || w[0] != w[6] || h[0] != h[6])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

        if (w[0] != w[7] || h[0] != h[7] || w[1] != w[6] || h[1] != h[6])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

        if (w[1] != w[7] || h[1] != h[7] || w[6] != w[7] || h[6] != h[7])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

    } else {
        if (w[0] != w[1] || h[0] != h[1])
            return fprintf(stderr, "ERROR: input images size mismatch\n");
    }

    // Ensure dimensions match in flow
    if (w[2] != w[3] || h[2] != h[3] || pd[2] != 2 || pd[2] != pd[3])
        return fprintf(stderr, "ERROR: input flow field size mismatch\n");
    //std::printf("Local faldoi line 1249\n");

    // Load or compute saliency
    float *sal0 = nullptr;
    float *sal1 = nullptr;
    if (args.size() == 9 || args.size() == 8) {
        sal0 = iio_read_image_float(filename_sal0, w + 4, h + 4);
        sal1 = iio_read_image_float(filename_sal1, w + 5, h + 5);
        fprintf(stderr, "Reading given saliency values\n");
    } else {
        fprintf(stderr, "Saliency values not given\n");
        sal0 = new float[w[0] * h[0]];
        sal1 = new float[w[0] * h[0]];
        for (int i = 0; i < w[0] * h[0]; i++) {
            sal0[i] = 1.0;
            sal1[i] = 1.0;
        }
    }
    //TODO: this for loop is redundant (the same one as in the 'else' above)
    /*for (int i = 0; i < w[0] * h[0]; i++) {
        sal0[i] = 1.0;
        sal1[i] = 1.0;
    }
     */
    for (int i = 0; i < w[0] * h[0]; i++) {
        assert(isfinite(sal0[i]));
        assert(isfinite(sal1[i]));
        if (sal0[i] < 0.0)
            fprintf(stderr, "sal0 index: %d\n", i);
        if (sal1[i] < 0.0)
            fprintf(stderr, "sal1 index: %d\n", i);
    }
//    std::printf("Local faldoi line 1280\n");

    // Initialize output optical flow and energy
    auto *out_flow = new float[w[0] * h[0] * 2];
    auto *out_occ = new float[w[0] * h[0]];
    auto *ene_val = new float[w[0] * h[0]];

    for (int i = 0; i < w[0] * h[0] * 2; i++) {
        out_flow[i] = NAN;
    }

    // Print method used
    if (num_files == 2 && val_method == M_TVL1_OCC) {
        // If only two images given for occ, something not working
        // TODO: when new methods with occlusions are implemented, add them here
        switch (val_method) {
            case M_TVL1_OCC:
                fprintf(stderr, "Since only two images given, method is changed to TV-l2 coupled\n");
                val_method = M_TVL1;
                break;
            default:
                fprintf(stderr, "Method unknown\n");
                break;
        }
    } else {
        // If four images given for without occ, two not needed
        if (num_files == 4 && val_method >= 0 && val_method <= 7) {
            fprintf(stderr, "Only two of the four images given will be used, according to method selected\n");
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
                    fprintf(stderr, "TV-l1 occlusions\n");
                    break;
                default:
                    break;
            }
        }
    }

    // Initialize parameters
    int step_alg = LOCAL_STEP;
    Parameters params = init_params(file_params, step_alg);
    params.w = w[0];
    params.h = h[0];
    params.pd = pd[0];
    params.w_radio = w_radio;
    params.val_method = val_method;
    params.iterations_of = loc_it;
    params.max_iter_patch = max_it_pch;
    params.split_img = sp_img;
    params.h_parts = h_prts;
    params.v_parts = v_prts;
    cerr << params;


/*    auto clk2 = system_clock::now(); // DEBUG
    duration<double> elapsed_secs = clk2 - clk1; // DEBUG
    cout << "(local_faldoi) Reading arguments (+images), preparing them and defining params for growing took "
         << elapsed_secs.count() << endl;*/

    auto clk1 = system_clock::now(); // DEBUG
    // Match growing algorithm
    match_growing_variational(go, ba, i0, i1, i_1, i2, sal0, sal1, params, ene_val, out_flow, out_occ);

    auto clk2 = system_clock::now(); // DEBUG
    duration<double> elapsed_secs2 = clk2 - clk1; // DEBUG
    cout << "(local_faldoi) Match growing variational took "
         << elapsed_secs2.count() << endl;

    // Save results

    iio_save_image_float_split(filename_out.c_str(), out_flow, w[0], h[0], 2);

    iio_save_image_float(filename_sim.c_str(), ene_val, w[0], h[0]);

    // Properly define occlusion mask
    if (args.size() == 7 || args.size() == 9) {
        //iio_save_image_float(filename_occ.c_str(), out_occ, w[0], h[0]);
        auto *out_occ_int = new int[w[0] * h[0]];

        for (int i = 0; i < w[0] * h[0]; i++) {

            out_occ_int[i] = out_occ[i];

        }

        iio_save_image_int(filename_occ.c_str(), out_occ_int, w[0], h[0]);
    }
    // Else, define it as null (to avoid compilation errors)
    int *out_occ_int = nullptr;

/*    auto clk4 = system_clock::now(); // DEBUG
    duration<double> elapsed_secs3 = clk4 - clk3; // DEBUG
    cout << "(local_faldoi) Saving all files to disk took "
         << elapsed_secs3.count() << endl;*/

    // Cleanup and exit
    free(i_1);
    free(i0);
    free(i1);
    free(i2);

    free(go);
    free(ba);

    if (args.size() == 8 || args.size() == 9) {
        free(sal0);
        free(sal1);
    } else {
        delete[] sal0;
        delete[] sal1;
    }

    delete[] out_flow;
    delete[] ene_val;

    if (args.size() == 7 || args.size() == 9) {
        delete[] out_occ;

        delete[] out_occ_int;
    }
    today = system_clock::now();

    tt = system_clock::to_time_t(today);
    cerr << "Finishing date: " << ctime(&tt);
    return 0;
}

#endif
