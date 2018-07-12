// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// TODO: add copyright from this year 2018 and gmail address
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


// Add 'about' and param explanation
static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return INFINITY;
    return x[(i + j * w) * pd + l];
}


// Add 'about' and param explanation
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


// Add 'about' and param explanation
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
// Add param explanation
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


// Add 'about'
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


// Add 'about' and param explanation
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
void  init_subimage_partitions(
        const float *i0,                                // Source image I0 (t)
        const float *i1,                                // Second image I1 (t+1)
        const float *i_1,                               // Previous image I-1 (t-1)
        const float *i2,                                // Third frame I2 (t+2)
        const float *i0n,                               // I0 normalised
        const float *i1n,                               // I1 normalised
        const float *i_1n,                              // I-1 normalised
        const float *i2n,                               // I2 normalised
        BilateralFilterData* BiFilt_Go,                 // Bilateral filter (fwd)
        BilateralFilterData* BiFilt_Ba,                 // Bilateral filter (bwd)
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
    if (rem_width > 0) {
        for (int i = 0; i < v_parts; i++) {
            sub_w[i * h_parts + h_parts - 1] += rem_width;
        }
    }

    // Vertically
    if (rem_height > 0) {
        for (int j = 0; j < h_parts; j++) {
            sub_h[(v_parts - 1) * h_parts + j] += rem_height;
        }
    }

    // Define partitions' offset w.r.t original image
    auto *off_x = new int[num_partitions];
    auto *off_y = new int[num_partitions];
    std::fill_n(off_x, num_partitions, 0);
    std::fill_n(off_y, num_partitions, 0);

    for (int v = 0; v < v_parts; v++)
        for (int h = 0; h < h_parts; h++) {
            if (h > 0) {
                off_x[v * h_parts + h] = h * sub_w[0];
            }
            if (v > 0) {
                off_y[v * h_parts + h] = v * sub_h[0];
            }
        }

    // No return value, just update the structs' fields via pointer
    for (int p = 0; p < num_partitions; p++) {
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

    // Fill images
    const int n_channels = params.pd;

    // Check that dimensions match
    int total_size = 0;
    for (unsigned p = 0; p < num_partitions; p++) {
        total_size += p_data->at(p)->width * p_data->at(p)->height;
    }
    assert(total_size == w_src * h_src);

    for (unsigned p = 0; p < num_partitions; p++) {
        int size = p_data->at(p)->width * p_data->at(p)->height;
        int size2 = size * n_channels;
        auto *i0_p = new float[size2];
        auto *i1_p = new float[size2];
        auto *i_1_p = new float[size2];
        auto *i2_p = new float[size2];
        // Same, but for normalised versions + dummy variables for 'prepare_stuff'
        // NOTE: prepare_stuff updates other fields but the normalised values should be copied from whole-image normalised
        auto *i0n_p = new float[size];
        auto *i1n_p = new float[size];
        auto *i_1n_p = new float[size];
        auto *i2n_p = new float[size];

        auto *i0n_dum = new float[size];
        auto *i1n_dum = new float[size];
        auto *i_1n_dum = new float[size];
        auto *i2n_dum = new float[size];


        for (int k = 0; k < n_channels; k++) {
            for (int j = 0; j < p_data->at(p)->height; j++)
                for (int i = 0; i < p_data->at(p)->width; i++) {
                    //int m = (j * p_data->at(p)->width + i) * n_channels + k;
                    int m = j * p_data->at(p)->width + i + k * p_data->at(p)->width * p_data->at(p)->height;
                    // idx of the subarray is computed as follows (from image array):
                    // idx = (y + j) * p + x + i
                    // where:   x, y (subimage offsets, top-left corner)
                    //          i, j indices (for sub_width[p] cols, sub_height[p] rows)
                    //          p is the width of the image (q is the height)
                    // int idx = (offset_y[p] + j) * w_src + offset_x[p] + i; // equiv for 2D
                    //int idx = ((p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i) * n_channels + k;  // " " 3D
                    int idx = (p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i + k * w_src * h_src;

                    if (k == 0) {
                        i0_p[m] = i0[idx];
                        i1_p[m] = i1[idx];
                        i_1_p[m] = i_1[idx];
                        i2_p[m] = i2[idx];

                        i0n_p[m] = i0n[idx];
                        i1n_p[m] = i1n[idx];

                        // Careful, i_1n and i2n are null for methods w/o occlusions modelling
                        if (params.val_method >= 8) {
                            i_1n_p[m] = i_1n[idx];
                            i2n_p[m] = i2n[idx];
                        }
                    } else {
                        i0_p[m] = i0[idx];
                        i1_p[m] = i1[idx];
                        i_1_p[m] = i_1[idx];
                        i2_p[m] = i2[idx];
                    }
                }
        }
        // Update vector of structs with correct values
        // Frames
        p_data->at(p)->i0 = i0_p;
        p_data->at(p)->i1 = i1_p;
        p_data->at(p)->i_1 = i_1_p;
        p_data->at(p)->i2 = i2_p;

        // Normalised frames
        p_data->at(p)->i0n = i0n_p;
        p_data->at(p)->i1n = i1n_p;
        p_data->at(p)->i_1n = i_1n_p;
        p_data->at(p)->i2n = i2n_p;

        // Optical flow data struct
        p_data->at(p)->ofGo = init_Optical_Flow_Data(sal_go, params, sub_w[p], sub_h[p]);
        p_data->at(p)->ofBa = init_Optical_Flow_Data(sal_ba, params, sub_w[p], sub_h[p]);

        // Initialise Specific OF stuff
        initialize_auxiliar_stuff(p_data->at(p)->stuffGo, p_data->at(p)->ofGo, sub_w[p], sub_h[p]);
        initialize_auxiliar_stuff(p_data->at(p)->stuffBa, p_data->at(p)->ofBa, sub_w[p], sub_h[p]);

        // Prepare auxiliar stuff
        prepare_stuff(&p_data->at(p)->stuffGo, &p_data->at(p)->ofGo, &p_data->at(p)->stuffBa, &p_data->at(p)->ofBa,
                      i0_p, i1_p, i_1_p, i2_p, params.pd, &i0n_dum, &i1n_dum, &i_1n_dum, &i2n_dum, sub_w[p], sub_h[p]);

        // TODO: merge fwd and bwd into one block of code (if possible)
        // At least try to put the 'partition' loop inside the w_src, h_src loop to reduce iterations (may be tricky)
        // BiFilt_Go (fwd)
        for (unsigned p = 0; p < num_partitions; p++) {
            // Create new BilateralFilterData object
            auto *Filter_data = new BilateralFilterData[1];
            Filter_data->weights_filtering = new Weights_Bilateral[p_data->at(p)->width * p_data->at(p)->height];

            Filter_data->indexes_filtering.resize(p_data->at(p)->width * p_data->at(p)->height);
            int min_i = p_data->at(p)->off_x;
            int max_i = min_i + p_data->at(p)->width - 1;
            int min_j = p_data->at(p)->off_y;
            int max_j = min_j + p_data->at(p)->height - 1;

            auto count = -1;
            // Loop through the original (complete image) BiFilt_Go and assign to p the corresponding indices
            for (int j = 0; j < h_src; j++)
                for (int i = 0; i < w_src; i++) {
                    const int ij = j * w_src + i;
                    auto neighbour = BiFilt_Go->indexes_filtering[ij];

                    // Check if this index belongs to the current 'p' partition
                    // If so, proceeed with assigning each value, otherwise check next iteration ('continue')
                    if ((neighbour.i >= min_i && neighbour.i <= max_i) &&
                        (neighbour.j >= min_j && neighbour.j <= max_j)) {
                        count++;
                        Filter_data->weights_filtering[count].weight = BiFilt_Go->weights_filtering[ij].weight;
                        Filter_data->indexes_filtering[count] = neighbour;

                    } else {
                        // Check following iter
                        continue;
                    }
                }

            // Assign to p_data
            p_data->at(p)->BiFilt_Go = Filter_data;
        }

        // BiFilt_Go (bwd)
        for (unsigned p = 0; p < num_partitions; p++) {
            // Create new BilateralFilterData object
            auto *Filter_data = new BilateralFilterData[1];
            Filter_data->weights_filtering = new Weights_Bilateral[p_data->at(p)->width * p_data->at(p)->height];

            Filter_data->indexes_filtering.resize(p_data->at(p)->width * p_data->at(p)->height);
            int min_i = p_data->at(p)->off_x;
            int max_i = min_i + p_data->at(p)->width - 1;
            int min_j = p_data->at(p)->off_y;
            int max_j = min_j + p_data->at(p)->height - 1;

            auto count = -1;
            for (int j = 0; j < h_src; j++)
                for (int i = 0; i < w_src; i++) {
                    const int ij = j * w_src + i;
                    auto neighbour = BiFilt_Ba->indexes_filtering[ij];

                    if ((neighbour.i >= min_i && neighbour.i <= max_i) &&
                        (neighbour.j >= min_j && neighbour.j <= max_j)) {
                        count++;
                        Filter_data->weights_filtering[count].weight = BiFilt_Ba->weights_filtering[ij].weight;
                        Filter_data->indexes_filtering[count] = neighbour;

                    } else {
                        // Check following iter
                        continue;
                    }
                }
            p_data->at(p)->BiFilt_Ba = Filter_data;
        }
    }
}


// Add 'about' and param explanation
void update_of_data(
        OpticalFlowData *& ofGo,
        OpticalFlowData *& ofBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        const int n_ch,
        bool img_to_part
) {
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones (Optical Flow Data)
        if (n_ch == 0) {
            // Update all variables
            // Chi
            p_data->ofGo.chi[idx_par] = ofGo->chi[idx_img];
            p_data->ofBa.chi[idx_par] = ofBa->chi[idx_img];

            // Fixed points and trust points
            p_data->ofGo.fixed_points[idx_par] = ofGo->fixed_points[idx_img];
            p_data->ofBa.fixed_points[idx_par] = ofBa->fixed_points[idx_img];
            p_data->ofGo.trust_points[idx_par] = ofGo->trust_points[idx_img];
            p_data->ofBa.trust_points[idx_par] = ofBa->trust_points[idx_img];

            // Saliency
            p_data->ofGo.saliency[idx_par] = ofGo->saliency[idx_img];
            p_data->ofBa.saliency[idx_par] = ofBa->saliency[idx_img];

            // OF fields
            p_data->ofGo.u1[idx_par] = ofGo->u1[idx_img];
            p_data->ofGo.u2[idx_par] = ofGo->u2[idx_img];
            p_data->ofBa.u1[idx_par] = ofBa->u1[idx_img];
            p_data->ofBa.u2[idx_par] = ofBa->u2[idx_img];

            p_data->ofGo.u1_ba[idx_par] = ofGo->u1_ba[idx_img];
            p_data->ofGo.u2_ba[idx_par] = ofGo->u2_ba[idx_img];
            p_data->ofBa.u1_ba[idx_par] = ofBa->u1_ba[idx_img];
            p_data->ofBa.u2_ba[idx_par] = ofBa->u2_ba[idx_img];

            // Filters
            p_data->ofGo.u1_filter[idx_par] = ofGo->u1_filter[idx_img];
            p_data->ofGo.u2_filter[idx_par] = ofGo->u2_filter[idx_img];
            p_data->ofBa.u1_filter[idx_par] = ofBa->u1_filter[idx_img];
            p_data->ofBa.u2_filter[idx_par] = ofBa->u2_filter[idx_img];

        } else {
            // Only update variables with 2 channels
            // OF fields
            p_data->ofGo.u1[idx_par] = ofGo->u1[idx_img];
            //p_data->ofGo.u2[idx_par] = ofGo->u2[idx_img];  // access out of range (only w*h elems)
            p_data->ofBa.u1[idx_par] = ofBa->u1[idx_img];
            //p_data->ofBa.u2[idx_par] = ofBa->u2[idx_img];  // access out of range (only w*h elems)

            p_data->ofGo.u1_ba[idx_par] = ofGo->u1_ba[idx_img];
            //p_data->ofGo.u2_ba[idx_par] = ofGo->u2_ba[idx_img];  // access out of range (only w*h elems)
            p_data->ofBa.u1_ba[idx_par] = ofBa->u1_ba[idx_img];
            //p_data->ofBa.u2_ba[idx_par] = ofBa->u2_ba[idx_img];  // access out of range (only w*h elems)

            // Filters
            p_data->ofGo.u1_filter[idx_par] = ofGo->u1_filter[idx_img];
            //p_data->ofGo.u2_filter[idx_par] = ofGo->u2_filter[idx_img];  // access out of range (only w*h elems)
            p_data->ofBa.u1_filter[idx_par] = ofBa->u1_filter[idx_img];
            //p_data->ofBa.u2_filter[idx_par] = ofBa->u2_filter[idx_img];  // access out of range (only w*h elems)
        }
    }
    else
    {
        // Copy partition-wise variables to corresponding image-wise variables
        if (n_ch == 0) {
            // Update all variables
            // Chi
            ofGo->chi[idx_img] = p_data->ofGo.chi[idx_par];
            ofBa->chi[idx_img] = p_data->ofBa.chi[idx_par];

            // Fixed points and trust points
            ofGo->fixed_points[idx_img] = p_data->ofGo.fixed_points[idx_par];
            ofBa->fixed_points[idx_img] = p_data->ofBa.fixed_points[idx_par];
            ofGo->trust_points[idx_img] = p_data->ofGo.trust_points[idx_par];
            ofBa->trust_points[idx_img] = p_data->ofBa.trust_points[idx_par];

            // Saliency
            ofGo->saliency[idx_img] = p_data->ofGo.saliency[idx_par];
            ofBa->saliency[idx_img] = p_data->ofBa.saliency[idx_par];

            // OF fields
            ofGo->u1[idx_img] = p_data->ofGo.u1[idx_par];
            ofGo->u2[idx_img] = p_data->ofGo.u2[idx_par];
            ofBa->u1[idx_img] = p_data->ofBa.u1[idx_par];
            ofBa->u2[idx_img] = p_data->ofBa.u2[idx_par];

            ofGo->u1_ba[idx_img] = p_data->ofGo.u1_ba[idx_par];
            ofGo->u2_ba[idx_img] = p_data->ofGo.u2_ba[idx_par];
            ofBa->u1_ba[idx_img] = p_data->ofBa.u1_ba[idx_par];
            ofBa->u2_ba[idx_img] = p_data->ofBa.u2_ba[idx_par];

            // Filters
            ofGo->u1_filter[idx_img] = p_data->ofGo.u1_filter[idx_par];
            ofGo->u2_filter[idx_img] = p_data->ofGo.u2_filter[idx_par];
            ofBa->u1_filter[idx_img] = p_data->ofBa.u1_filter[idx_par];
            ofBa->u2_filter[idx_img] = p_data->ofBa.u2_filter[idx_par];

        } else {
            // Only update variables with 2 channels
            // OF fields
            ofGo->u1[idx_img] = p_data->ofGo.u1[idx_par];
            //ofGo->u2[idx_img] = p_data->ofGo.u2[idx_par];  // access out of range (only w*h elems)
            ofBa->u1[idx_img] = p_data->ofBa.u1[idx_par];
            //ofBa->u2[idx_img] = p_data->ofBa.u2[idx_par];  // "   "   "   "

            ofGo->u1_ba[idx_img] = p_data->ofGo.u1_ba[idx_par];
            //ofGo->u2_ba[idx_img] = p_data->ofGo.u2_ba[idx_par];  // "   "   "   "
            ofBa->u1_ba[idx_img] = p_data->ofBa.u1_ba[idx_par];
            //ofBa->u2_ba[idx_img] = p_data->ofBa.u2_ba[idx_par];  // "   "   "   "

            // Filters
            ofGo->u1_filter[idx_img] = p_data->ofGo.u1_filter[idx_par];
            //ofGo->u2_filter[idx_img] = p_data->ofGo.u2_filter[idx_par];  // "   "   "   "
            ofBa->u1_filter[idx_img] = p_data->ofBa.u1_filter[idx_par];
            //ofBa->u2_filter[idx_img] = p_data->ofBa.u2_filter[idx_par];  // "   "   "   "
        }
    }
}


// Add 'about' and param explanation
void update_tvl2_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Xi
        p_data->stuffGo.tvl2.xi11[idx_par] = stuffGo->tvl2.xi11[idx_img];
        p_data->stuffGo.tvl2.xi12[idx_par] = stuffGo->tvl2.xi12[idx_img];
        p_data->stuffGo.tvl2.xi21[idx_par] = stuffGo->tvl2.xi21[idx_img];
        p_data->stuffGo.tvl2.xi22[idx_par] = stuffGo->tvl2.xi22[idx_img];

        p_data->stuffBa.tvl2.xi11[idx_par] = stuffBa->tvl2.xi11[idx_img];
        p_data->stuffBa.tvl2.xi12[idx_par] = stuffBa->tvl2.xi12[idx_img];
        p_data->stuffBa.tvl2.xi21[idx_par] = stuffBa->tvl2.xi21[idx_img];
        p_data->stuffBa.tvl2.xi22[idx_par] = stuffBa->tvl2.xi22[idx_img];

        // u1, u2
        p_data->stuffGo.tvl2.u1x[idx_par] = stuffGo->tvl2.u1x[idx_img];
        p_data->stuffGo.tvl2.u1y[idx_par] = stuffGo->tvl2.u1y[idx_img];
        p_data->stuffGo.tvl2.u2x[idx_par] = stuffGo->tvl2.u2x[idx_img];
        p_data->stuffGo.tvl2.u2y[idx_par] = stuffGo->tvl2.u2y[idx_img];

        p_data->stuffBa.tvl2.u1x[idx_par] = stuffBa->tvl2.u1x[idx_img];
        p_data->stuffBa.tvl2.u1y[idx_par] = stuffBa->tvl2.u1y[idx_img];
        p_data->stuffBa.tvl2.u2x[idx_par] = stuffBa->tvl2.u2x[idx_img];
        p_data->stuffBa.tvl2.u2y[idx_par] = stuffBa->tvl2.u2y[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.tvl2.v1[idx_par] = stuffGo->tvl2.v1[idx_img];
        p_data->stuffGo.tvl2.v2[idx_par] = stuffGo->tvl2.v2[idx_img];

        p_data->stuffBa.tvl2.v1[idx_par] = stuffBa->tvl2.v1[idx_img];
        p_data->stuffBa.tvl2.v2[idx_par] = stuffBa->tvl2.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.tvl2.rho_c[idx_par] = stuffGo->tvl2.rho_c[idx_img];
        p_data->stuffBa.tvl2.rho_c[idx_par] = stuffBa->tvl2.rho_c[idx_img];

        p_data->stuffGo.tvl2.grad[idx_par] = stuffGo->tvl2.grad[idx_img];
        p_data->stuffBa.tvl2.grad[idx_par] = stuffBa->tvl2.grad[idx_img];

        p_data->stuffGo.tvl2.u1_[idx_par] = stuffGo->tvl2.u1_[idx_img];
        p_data->stuffGo.tvl2.u2_[idx_par] = stuffGo->tvl2.u2_[idx_img];

        p_data->stuffBa.tvl2.u1_[idx_par] = stuffBa->tvl2.u1_[idx_img];
        p_data->stuffBa.tvl2.u2_[idx_par] = stuffBa->tvl2.u2_[idx_img];

        p_data->stuffGo.tvl2.u1Aux[idx_par] = stuffGo->tvl2.u1Aux[idx_img];
        p_data->stuffGo.tvl2.u2Aux[idx_par] = stuffGo->tvl2.u2Aux[idx_img];

        p_data->stuffBa.tvl2.u1Aux[idx_par] = stuffBa->tvl2.u1Aux[idx_img];
        p_data->stuffBa.tvl2.u2Aux[idx_par] = stuffBa->tvl2.u2Aux[idx_img];

        p_data->stuffGo.tvl2.I1x[idx_par] = stuffGo->tvl2.I1x[idx_img];
        p_data->stuffGo.tvl2.I1y[idx_par] = stuffGo->tvl2.I1y[idx_img];

        p_data->stuffBa.tvl2.I1x[idx_par] = stuffBa->tvl2.I1x[idx_img];
        p_data->stuffBa.tvl2.I1y[idx_par] = stuffBa->tvl2.I1y[idx_img];

        p_data->stuffGo.tvl2.I1wx[idx_par] = stuffGo->tvl2.I1wx[idx_img];
        p_data->stuffGo.tvl2.I1wy[idx_par] = stuffGo->tvl2.I1wy[idx_img];

        p_data->stuffBa.tvl2.I1wx[idx_par] = stuffBa->tvl2.I1wx[idx_img];
        p_data->stuffBa.tvl2.I1wy[idx_par] = stuffBa->tvl2.I1wy[idx_img];

        p_data->stuffGo.tvl2.div_xi1[idx_par] = stuffGo->tvl2.div_xi1[idx_img];
        p_data->stuffGo.tvl2.div_xi2[idx_par] = stuffGo->tvl2.div_xi2[idx_img];

        p_data->stuffBa.tvl2.div_xi1[idx_par] = stuffBa->tvl2.div_xi1[idx_img];
        p_data->stuffBa.tvl2.div_xi2[idx_par] = stuffBa->tvl2.div_xi2[idx_img];

        p_data->stuffGo.tvl2.u_N[idx_par] = stuffGo->tvl2.u_N[idx_img];
        p_data->stuffBa.tvl2.u_N[idx_par] = stuffBa->tvl2.u_N[idx_img];
    }
    else
    {
        // Copy partition-wise variables to corresponding image-wise variables
        // Xi
        stuffGo->tvl2.xi11[idx_img] = p_data->stuffGo.tvl2.xi11[idx_par];
        stuffGo->tvl2.xi12[idx_img] = p_data->stuffGo.tvl2.xi12[idx_par];
        stuffGo->tvl2.xi21[idx_img] = p_data->stuffGo.tvl2.xi21[idx_par];
        stuffGo->tvl2.xi22[idx_img] = p_data->stuffGo.tvl2.xi22[idx_par];

        stuffBa->tvl2.xi11[idx_img] = p_data->stuffBa.tvl2.xi11[idx_par];
        stuffBa->tvl2.xi12[idx_img] = p_data->stuffBa.tvl2.xi12[idx_par];
        stuffBa->tvl2.xi21[idx_img] = p_data->stuffBa.tvl2.xi21[idx_par];
        stuffBa->tvl2.xi22[idx_img] = p_data->stuffBa.tvl2.xi22[idx_par];

        // u1, u2
        stuffGo->tvl2.u1x[idx_img] = p_data->stuffGo.tvl2.u1x[idx_par];
        stuffGo->tvl2.u1y[idx_img] = p_data->stuffGo.tvl2.u1y[idx_par];
        stuffGo->tvl2.u2x[idx_img] = p_data->stuffGo.tvl2.u2x[idx_par];
        stuffGo->tvl2.u2y[idx_img] = p_data->stuffGo.tvl2.u2y[idx_par];

        stuffBa->tvl2.u1x[idx_img] = p_data->stuffBa.tvl2.u1x[idx_par];
        stuffBa->tvl2.u1y[idx_img] = p_data->stuffBa.tvl2.u1y[idx_par];
        stuffBa->tvl2.u2x[idx_img] = p_data->stuffBa.tvl2.u2x[idx_par];
        stuffBa->tvl2.u2y[idx_img] = p_data->stuffBa.tvl2.u2y[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->tvl2.v1[idx_img] = p_data->stuffGo.tvl2.v1[idx_par];
        stuffGo->tvl2.v2[idx_img] = p_data->stuffGo.tvl2.v2[idx_par];

        stuffBa->tvl2.v1[idx_img] = p_data->stuffBa.tvl2.v1[idx_par];
        stuffBa->tvl2.v2[idx_img] = p_data->stuffBa.tvl2.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->tvl2.rho_c[idx_img] = p_data->stuffGo.tvl2.rho_c[idx_par];
        stuffBa->tvl2.rho_c[idx_img] = p_data->stuffBa.tvl2.rho_c[idx_par];

        stuffGo->tvl2.grad[idx_img] = p_data->stuffGo.tvl2.grad[idx_par];
        stuffBa->tvl2.grad[idx_img] = p_data->stuffBa.tvl2.grad[idx_par];

        stuffGo->tvl2.u1_[idx_img] = p_data->stuffGo.tvl2.u1_[idx_par];
        stuffGo->tvl2.u2_[idx_img] = p_data->stuffGo.tvl2.u2_[idx_par];

        stuffBa->tvl2.u1_[idx_img] = p_data->stuffBa.tvl2.u1_[idx_par];
        stuffBa->tvl2.u2_[idx_img] = p_data->stuffBa.tvl2.u2_[idx_par];

        stuffGo->tvl2.u1Aux[idx_img] = p_data->stuffGo.tvl2.u1Aux[idx_par];
        stuffGo->tvl2.u2Aux[idx_img] = p_data->stuffGo.tvl2.u2Aux[idx_par];

        stuffBa->tvl2.u1Aux[idx_img] = p_data->stuffBa.tvl2.u1Aux[idx_par];
        stuffBa->tvl2.u2Aux[idx_img] = p_data->stuffBa.tvl2.u2Aux[idx_par];

        stuffGo->tvl2.I1x[idx_img] = p_data->stuffGo.tvl2.I1x[idx_par];
        stuffGo->tvl2.I1y[idx_img] = p_data->stuffGo.tvl2.I1y[idx_par];

        stuffBa->tvl2.I1x[idx_img] = p_data->stuffBa.tvl2.I1x[idx_par];
        stuffBa->tvl2.I1y[idx_img] = p_data->stuffBa.tvl2.I1y[idx_par];

        stuffGo->tvl2.I1wx[idx_img] = p_data->stuffGo.tvl2.I1wx[idx_par];
        stuffGo->tvl2.I1wy[idx_img] = p_data->stuffGo.tvl2.I1wy[idx_par];

        stuffBa->tvl2.I1wx[idx_img] = p_data->stuffBa.tvl2.I1wx[idx_par];
        stuffBa->tvl2.I1wy[idx_img] = p_data->stuffBa.tvl2.I1wy[idx_par];

        stuffGo->tvl2.div_xi1[idx_img] = p_data->stuffGo.tvl2.div_xi1[idx_par];
        stuffGo->tvl2.div_xi2[idx_img] = p_data->stuffGo.tvl2.div_xi2[idx_par];

        stuffBa->tvl2.div_xi1[idx_img] = p_data->stuffBa.tvl2.div_xi1[idx_par];
        stuffBa->tvl2.div_xi2[idx_img] = p_data->stuffBa.tvl2.div_xi2[idx_par];

        stuffGo->tvl2.u_N[idx_img] = p_data->stuffGo.tvl2.u_N[idx_par];
        stuffBa->tvl2.u_N[idx_img] = p_data->stuffBa.tvl2.u_N[idx_par];
    }
}


// TODO: implement the rest once it works with the TVL2 functional
// Add 'about' and param explanation and implementation
void update_nltvl1_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Dual variables
        p_data->stuffGo.nltvl1.p[idx_par] = stuffGo->nltvl1.p[idx_img];
        p_data->stuffGo.nltvl1.q[idx_par] = stuffGo->nltvl1.q[idx_img];

        p_data->stuffBa.nltvl1.p[idx_par] = stuffBa->nltvl1.p[idx_img];
        p_data->stuffBa.nltvl1.q[idx_par] = stuffBa->nltvl1.q[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.nltvl1.v1[idx_par] = stuffGo->nltvl1.v1[idx_img];
        p_data->stuffGo.nltvl1.v2[idx_par] = stuffGo->nltvl1.v2[idx_img];

        p_data->stuffBa.nltvl1.v1[idx_par] = stuffBa->nltvl1.v1[idx_img];
        p_data->stuffBa.nltvl1.v2[idx_par] = stuffBa->nltvl1.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.nltvl1.rho_c[idx_par] = stuffGo->nltvl1.rho_c[idx_img];
        p_data->stuffBa.nltvl1.rho_c[idx_par] = stuffBa->nltvl1.rho_c[idx_img];

        p_data->stuffGo.nltvl1.grad[idx_par] = stuffGo->nltvl1.grad[idx_img];
        p_data->stuffBa.nltvl1.grad[idx_par] = stuffBa->nltvl1.grad[idx_img];

        p_data->stuffGo.nltvl1.u1_[idx_par] = stuffGo->nltvl1.u1_[idx_img];
        p_data->stuffGo.nltvl1.u2_[idx_par] = stuffGo->nltvl1.u2_[idx_img];

        p_data->stuffBa.nltvl1.u1_[idx_par] = stuffBa->nltvl1.u1_[idx_img];
        p_data->stuffBa.nltvl1.u2_[idx_par] = stuffBa->nltvl1.u2_[idx_img];

        p_data->stuffGo.nltvl1.u1_tmp[idx_par] = stuffGo->nltvl1.u1_tmp[idx_img];
        p_data->stuffGo.nltvl1.u2_tmp[idx_par] = stuffGo->nltvl1.u2_tmp[idx_img];

        p_data->stuffBa.nltvl1.u1_tmp[idx_par] = stuffBa->nltvl1.u1_tmp[idx_img];
        p_data->stuffBa.nltvl1.u2_tmp[idx_par] = stuffBa->nltvl1.u2_tmp[idx_img];

        p_data->stuffGo.nltvl1.I1x[idx_par] = stuffGo->nltvl1.I1x[idx_img];
        p_data->stuffGo.nltvl1.I1y[idx_par] = stuffGo->nltvl1.I1y[idx_img];

        p_data->stuffBa.nltvl1.I1x[idx_par] = stuffBa->nltvl1.I1x[idx_img];
        p_data->stuffBa.nltvl1.I1y[idx_par] = stuffBa->nltvl1.I1y[idx_img];

        p_data->stuffGo.nltvl1.I1wx[idx_par] = stuffGo->nltvl1.I1wx[idx_img];
        p_data->stuffGo.nltvl1.I1wy[idx_par] = stuffGo->nltvl1.I1wy[idx_img];

        p_data->stuffBa.nltvl1.I1wx[idx_par] = stuffBa->nltvl1.I1wx[idx_img];
        p_data->stuffBa.nltvl1.I1wy[idx_par] = stuffBa->nltvl1.I1wy[idx_img];

        p_data->stuffGo.tvl2.div_xi1[idx_par] = stuffGo->tvl2.div_xi1[idx_img];
        p_data->stuffGo.tvl2.div_xi2[idx_par] = stuffGo->tvl2.div_xi2[idx_img];

        p_data->stuffBa.tvl2.div_xi1[idx_par] = stuffBa->tvl2.div_xi1[idx_img];
        p_data->stuffBa.tvl2.div_xi2[idx_par] = stuffBa->tvl2.div_xi2[idx_img];

        p_data->stuffGo.tvl2.u_N[idx_par] = stuffGo->tvl2.u_N[idx_img];
        p_data->stuffBa.tvl2.u_N[idx_par] = stuffBa->tvl2.u_N[idx_img];
    }
    else
    {
        // Copy partition-wise variables to corresponding image-wise variables
        // Dual variables
        stuffGo->nltvl1.p[idx_img] = p_data->stuffGo.nltvl1.p[idx_par];
        stuffGo->nltvl1.q[idx_img] = p_data->stuffGo.nltvl1.q[idx_par];

        stuffBa->nltvl1.p[idx_img] = p_data->stuffBa.nltvl1.p[idx_par];
        stuffBa->nltvl1.q[idx_img] = p_data->stuffBa.nltvl1.q[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->nltvl1.v1[idx_img] = p_data->stuffGo.nltvl1.v1[idx_par];
        stuffGo->nltvl1.v2[idx_img] = p_data->stuffGo.nltvl1.v2[idx_par];

        stuffBa->nltvl1.v1[idx_img] = p_data->stuffBa.nltvl1.v1[idx_par];
        stuffBa->nltvl1.v2[idx_img] = p_data->stuffBa.nltvl1.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->nltvl1.rho_c[idx_img] = p_data->stuffGo.nltvl1.rho_c[idx_par];
        stuffBa->nltvl1.rho_c[idx_img] = p_data->stuffBa.nltvl1.rho_c[idx_par];

        stuffGo->nltvl1.grad[idx_img] = p_data->stuffGo.nltvl1.grad[idx_par];
        stuffBa->nltvl1.grad[idx_img] = p_data->stuffBa.nltvl1.grad[idx_par];

        stuffGo->nltvl1.u1_[idx_img] = p_data->stuffGo.nltvl1.u1_[idx_par];
        stuffGo->nltvl1.u2_[idx_img] = p_data->stuffGo.nltvl1.u2_[idx_par];

        stuffBa->nltvl1.u1_[idx_img] = p_data->stuffBa.nltvl1.u1_[idx_par];
        stuffBa->nltvl1.u2_[idx_img] = p_data->stuffBa.nltvl1.u2_[idx_par];

        stuffGo->nltvl1.u1_tmp[idx_img] = p_data->stuffGo.nltvl1.u1_tmp[idx_par];
        stuffGo->nltvl1.u2_tmp[idx_img] = p_data->stuffGo.nltvl1.u2_tmp[idx_par];

        stuffBa->nltvl1.u1_tmp[idx_img] = p_data->stuffBa.nltvl1.u1_tmp[idx_par];
        stuffBa->nltvl1.u2_tmp[idx_img] = p_data->stuffBa.nltvl1.u2_tmp[idx_par];

        stuffGo->nltvl1.I1x[idx_img] = p_data->stuffGo.nltvl1.I1x[idx_par];
        stuffGo->nltvl1.I1y[idx_img] = p_data->stuffGo.nltvl1.I1y[idx_par];

        stuffBa->nltvl1.I1x[idx_img] = p_data->stuffBa.nltvl1.I1x[idx_par];
        stuffBa->nltvl1.I1y[idx_img] = p_data->stuffBa.nltvl1.I1y[idx_par];

        stuffGo->nltvl1.I1wx[idx_img] = p_data->stuffGo.nltvl1.I1wx[idx_par];
        stuffGo->nltvl1.I1wy[idx_img] = p_data->stuffGo.nltvl1.I1wy[idx_par];

        stuffBa->nltvl1.I1wx[idx_img] = p_data->stuffBa.nltvl1.I1wx[idx_par];
        stuffBa->nltvl1.I1wy[idx_img] = p_data->stuffBa.nltvl1.I1wy[idx_par];

        stuffGo->nltvl1.div_p[idx_img] = p_data->stuffGo.nltvl1.div_p[idx_par];
        stuffGo->nltvl1.div_q[idx_img] = p_data->stuffGo.nltvl1.div_q[idx_par];

        stuffBa->nltvl1.div_p[idx_img] = p_data->stuffBa.nltvl1.div_p[idx_par];
        stuffBa->nltvl1.div_q[idx_img] = p_data->stuffBa.nltvl1.div_q[idx_par];
    }

}


// Add 'about' and param explanation and implementation
void update_tvl2w_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation and implementation
void update_tvl2occ_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {


}


// Add 'about' and param explanation and implementation
void update_nltv1w_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation and implementation
void update_tvcsad_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation and implementation
void update_tvcsadw_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation and implementation
void update_nltvcsad_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation and implementation
void update_nltvcsadw_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {

}


// Add 'about' and param explanation
void update_partitions_structures(
        OpticalFlowData *& ofGo,
        OpticalFlowData *& ofBa,
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        const int n_ch,
        const int idx_img,
        const int idx_par,
        bool img_to_part
) {
    // Update both structures when we are filling the 'first' channel
    // For the second channel, only certain members of OpticalFlowData have to be filled
    if (n_ch == 0) {
        // Update OpticalFlowData which is common to any functional
        update_of_data(ofGo, ofBa, p_data, idx_img, idx_par, n_ch, img_to_part);
        switch (ofGo->params.val_method) {
            case M_NLTVL1:          // NLTV-L1
                update_nltvl1_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVCSAD:          // TV-CSAD
                update_tvcsad_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_NLTVCSAD:        // NLTV-CSAD
                update_nltvcsad_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVL1_W:          // TV-l2 coupled with weights
                update_tvl2w_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_NLTVCSAD_W:      // NLTV-CSAD with weights
                update_nltvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_NLTVL1_W:        // NLTV-L1 with weights
                update_nltv1w_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVCSAD_W:        // TV-CSAD with weights
                update_tvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVL1_OCC:        // TV-l2 with occlusion
                update_tvl2occ_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            default:                // TV-l2 coupled
                update_tvl2_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
        }
    }
    else {
        update_of_data(ofGo, ofBa, p_data, idx_img, idx_par, n_ch, img_to_part);
    }
}


// Update queues (FWD + BWD) with potential new candidates (image-wise to partition-specific)
// Add param explanation
void update_candidate_queues(
        pq_cand &queue_Go,
        pq_cand &queue_Ba,
        const int n_partitions,
        std::vector<PartitionData*> *p_data
) {
    // Local copies of the original queues so we can 'iterate' trough them without creating other std::queue's
    pq_cand tmp_queueGo = queue_Go;
    pq_cand tmp_queueBa = queue_Ba;

    // Forward queue 'queue_Go'
    while (!tmp_queueGo.empty()) {
        SparseOF element = tmp_queueGo.top();
        int i = element.i;
        int j = element.j;

        for (unsigned p = 0; p < n_partitions; p++) {
            // Temporal variables to simplify 'if' statement
            int min_i = p_data->at(p)->off_x;
            int max_i = min_i + p_data->at(p)->width - 1;
            int min_j = p_data->at(p)->off_y;
            int max_j = min_j + p_data->at(p)->height - 1;

            if ((i >= min_i && i <= max_i) && (j >= min_j && j <= max_j)) {
                // Belongs to partition 'p', update the idx so it is partition-relative (not image-specific)
                element.i = i >= p_data->at(p)->width ? i - p_data->at(p)->off_x : i;
                element.j = j >= p_data->at(p)->height ? j - p_data->at(p)->off_y : j;

                // Add to correct queue
                p_data->at(p)->queue_Go.push(element);
                p_data->at(p)->queue_Go_size++;
            }
            // Otherwise, do nothing, check following partition
        }
        // While the queue is not empty, take an element to process
        tmp_queueGo.pop();
    }

    // Backward queue 'queue_Ba' (identical)
    while (!tmp_queueBa.empty()) {
        SparseOF element = tmp_queueBa.top();
        int i = element.i;
        int j = element.j;

        for (unsigned p = 0; p < n_partitions; p++) {
            // Temporal variables to simplify 'if' statement
            int min_i = p_data->at(p)->off_x;
            int max_i = min_i + p_data->at(p)->width - 1;
            int min_j = p_data->at(p)->off_y;
            int max_j = min_j + p_data->at(p)->height - 1;

            if ((i >= min_i && i <= max_i) && (j >= min_j && j <= max_j)) {
                // Belongs to partition 'p', update the idx so it is partition-relative (not image-specific)
                element.i = i >= p_data->at(p)->width ? i - p_data->at(p)->off_x : i;
                element.j = j >= p_data->at(p)->height ? j - p_data->at(p)->off_y : j;
                p_data->at(p)->queue_Ba.push(element);
                p_data->at(p)->queue_Ba_size++;
            }
            // Otherwise, do nothing, check following partition
        }
        // While the queue is not empty, take an element to process
        tmp_queueBa.pop();
    }
}


// Copy values from image-wise variables to each corresponding partition
void image_to_partitions(
        float *oft0,
        float *oft1,
        float *ene_Go,
        float *ene_Ba,
        float *occ_Go,
        float *occ_Ba,
        OpticalFlowData *ofGo,
        OpticalFlowData *ofBa,
        SpecificOFStuff *stuffGo,
        SpecificOFStuff *stuffBa,
        pq_cand &queue_Go,
        pq_cand &queue_Ba,
        const int n_partitions,
        const int w_src,
        const int h_src,
        std::vector<PartitionData*> *p_data,
        bool img_to_part                     // the 'direction' of the update (copy) of params
) {
    // Update variables
    int n_channels = 2;
    for (unsigned p = 0; p < n_partitions; p++) {
        int size = p_data->at(p)->width * p_data->at(p)->height;
        int size2 = size * n_channels;

        // Define temporal variables and initialise them with the same default values as 'prepare_data_for_growing'
        auto *ene_Go_p = new float[size];
        auto *ene_Ba_p = new float[size];
        auto *occ_Go_p = new float[size];
        auto *occ_Ba_p = new float[size];
        auto *oft0_p = new float[size2];
        auto *oft1_p = new float[size2];
        std::fill_n(ene_Go_p, size, INFINITY);
        std::fill_n(ene_Ba_p, size, INFINITY);
        std::fill_n(oft0_p, size2, NAN);
        std::fill_n(oft1_p, size2, NAN);
        // Leave occ as is if they are being computed
        if (ofGo->params.val_method < 8) {
            std::fill_n(occ_Go_p, size, 0);
            std::fill_n(occ_Ba_p, size, 0);
        }

        for (int k = 0; k < n_channels; k++)
            for (int j = 0; j < p_data->at(p)->height; j++)
                for (int i = 0; i < p_data->at(p)->width; i++) {
                    // 'Mapping' indices
                    // NOTE: we have to be careful, we stored multiple channel variables in the usual manner:
                    //  each channel is side by side, with alternate values
                    // In this code, multiple channels are stacked one after the other without alternate values
                    // So, instead of: V1_C1, V1_C2, V2_C1, V2_C1, ..., VN_C1, VN_C2; we have:
                    //                 V1_C1, V2_C1, ..., VN_C1, V1_C2, V2_C2, ..., VN_C2
                    //int m = (j * p_data->at(p)->width + i) * n_channels + k;
                    int m = j * p_data->at(p)->width + i + k * p_data->at(p)->width * p_data->at(p)->height;
                    //int idx = ((p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i) * n_channels + k; // + offset_p;
                    int idx = (p_data->at(p)->off_y + j) * w_src + p_data->at(p)->off_x + i + k * w_src * h_src;

                    if (img_to_part) {
                        if (k < n_channels - 1) {
                            // Update everything here
                            ene_Go_p[m] = ene_Go[idx];
                            ene_Ba_p[m] = ene_Ba[idx];
                            occ_Go_p[m] = occ_Go[idx];
                            occ_Ba_p[m] = occ_Ba[idx];
                            oft0_p[m] = oft0[idx];
                            oft1_p[m] = oft1[idx];

                        } else {
                            // Only update variables with more than a channel
                            oft0_p[m] = oft0[idx];
                            oft1_p[m] = oft1[idx];
                        }
                    } else {
                        if (k < n_channels - 1) {
                            ene_Go[idx] = p_data->at(p)->ene_Go[m];
                            ene_Ba[idx] = p_data->at(p)->ene_Ba[m];
                            occ_Go[idx] = p_data->at(p)->occ_Go[m];
                            occ_Ba[idx] = p_data->at(p)->occ_Ba[m];
                            oft0[idx] = p_data->at(p)->oft0[m];
                            oft1[idx] = p_data->at(p)->oft1[m];

                        } else {

                            oft0[idx] = p_data->at(p)->oft0[m];
                            oft1[idx] = p_data->at(p)->oft1[m];
                        }
                    }
                    update_partitions_structures(ofGo, ofBa, stuffGo, stuffBa, p_data->at(p), k, idx, m, img_to_part);
                }

        if (img_to_part)
        {
            p_data->at(p)->ene_Go = ene_Go_p;
            p_data->at(p)->ene_Ba = ene_Ba_p;
            p_data->at(p)->occ_Go = occ_Go_p;
            p_data->at(p)->occ_Ba = occ_Ba_p;
            p_data->at(p)->oft0 = oft0_p;
            p_data->at(p)->oft1 = oft1_p;
            p_data->at(p)->ofGo.params = ofGo->params;
            p_data->at(p)->ofBa.params = ofBa->params;
        }
        else {
            ofGo->params = p_data->at(p)->ofGo.params;
            ofBa->params = p_data->at(p)->ofBa.params;
        }
    }

    if (img_to_part) {
        // Update queues only if we are going back to partitions (after 'insert_potential_candidates')
        update_candidate_queues(queue_Go, queue_Ba, n_partitions, p_data);
    } else {
        // Otherwise, empty image-wise queues so only the new candidates are inserted
        pq_cand empty_queue_Go, empty_queue_Ba;
        std::swap(queue_Go, empty_queue_Go);
        std::swap(queue_Ba, empty_queue_Ba);
    }
}


bool anyEmptyQueues(
        std::vector<PartitionData*> *p_data,
        const int n_partitions
) {
    int n_empty_queues = 0;
    for (unsigned n = 0; n < n_partitions; n++) {
        if (p_data->at(n)->queue_Go_size == 0) {
            n_empty_queues++;
        }

        if (p_data->at(n)->queue_Ba_size == 0) {
            n_empty_queues++;
        }

        if (n_empty_queues > 0) {
            // Enough to revert back to images
            return true;
        }
        // Otherwise we keep checking other partitions
    }

    return 0 != n_empty_queues;  // Returns 'false' if 'n_emtpy_queues' = 0; 'true' otherwise
}


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
            } else {  // not fixed
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


// Add 'about' and param explanation
void bilateral_filter(
        OpticalFlowData *ofD,
        BilateralFilterData *BiFilt,
        const PatchIndexes &patch,
        const int w,
        const int h
) {

    int *trust_points = ofD->trust_points;
    int *fixed_points = ofD->fixed_points;

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
// Add param explanation
void insert_candidates(
        pq_cand &queue,
        float *ene_val,
        OpticalFlowData *ofD,
        const int i,
        const int j,
        const float ener_N,
        const int w,
        const int h
) {

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

    const float *sal = ofD->saliency;

    for (int k = 0; k < n_neigh; k++) {
        int px = i + neighborhood[k][0];
        int py = j + neighborhood[k][1];

        if (px >= 0 && px < w && py >= 0 && py < h) {
            float new_ener = ener_N * sal[py * w + px];

            //printf("Ener_N: %f  Sim: %f \n", ener_N, ene_val[py*w + px]);
            if (!ofD->fixed_points[py * w + px] && new_ener < ene_val[py * w + px]) {

                ene_val[py * w + px] = ener_N;
                SparseOF element{};
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
// Add 'about' and param explanation
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
        int j
) {

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
// Add param explanation
inline void copy_fixed_coordinates(
        OpticalFlowData *ofD,
        const float *out,
        const PatchIndexes &index,
        const int w,
        const int h
) {
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
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
// Add param explanation
int check_trustable_patch(
        OpticalFlowData *ofD,
        const PatchIndexes &index,
        const int w
) {
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

// Add 'about' and param explanation
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
        const int h
) {
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
// from the energy stored at the moment that the pixel was fixed.
void insert_potential_candidates(
        const float *in,
        OpticalFlowData *ofD,
        pq_cand &queue,
        const float *ene_val,
        const float *out_occ,
        const int w,
        const int h
) {
    // Fixed the initial seeds.
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            // Indicates the initial seed in the similarity map
            if (isfinite(in[j * w + i]) && isfinite(in[w * h + j * w + i])) {

                SparseOF element{};
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
        const int h
) {
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

            // Update output variables with the latest flow values (u, v)
            out_flow[j * w + i] = u;
            out_flow[w * h + j * w + i] = v;
            ene_val[j * w + i] = energy;

            out_occ[j * w + i] = occlusion;

            // Add new neigbours to the queue and compute the local minimization (patch-wise)
            // The output values may be updated if the new energy is better (lower) than the stored one
            add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, iteration, out_flow, out_occ, BiFilt, w, h);

            // From here to the end of the function:
            // Code used to print partial growing results for debugging or further exploration
            // Just set SAVE_RESULTS in parameters.h to 1 after creating the folder Partial_results under ../Results
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

                            for (int l = 0; l < w * h; l++) {

                                out_occ_int[l] = out_occ[l];
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
    using namespace chrono;  // PROFILING
    auto clk_matchGrow_init = system_clock::now();
    int w = params.w;
    int h = params.h;

    printf("Initializing stuff\n");
    // Initialize all the variables specific to the actual optical flow computation
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
    pq_cand queue_Go;
    pq_cand queue_Ba;

    // Initialize all the auxiliar data.
    SpecificOFStuff stuffGo{};
    SpecificOFStuff stuffBa{};
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


    auto clk_init_stuff = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_init_stuff = clk_init_stuff - clk_matchGrow_init; // PROFILING
    cout << "(match growing) Initializing everything took "
         << elapsed_secs_init_stuff.count() << endl;

    // Insert initial seeds to queues
    printf("Inserting initial seeds\n");

    auto clk_seeds = system_clock::now(); // PROFILING
#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            insert_initial_seeds(i0n, i1n, i_1n, go, &queue_Go, &ofGo, &stuffGo, ene_Go, oft0, occ_Go, BiFilt_Go,
                                 w, h);
        } else {
            insert_initial_seeds(i1n, i0n, i2n, ba, &queue_Ba, &ofBa, &stuffBa, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                 w, h);
        }
    }
#endif
    printf("Finished inserting initial seeds\n");

    auto clk_seeds_end = system_clock::now(); // PROFILING
    duration<double> elapsed_secs_seeds = clk_seeds_end - clk_seeds; // PROFILING
    cout << "(match growing) inserting initial seeds took "
         << elapsed_secs_seeds.count() << endl;

    const int iter = params.iterations_of;  // LOCAL_ITER;
    // Variables for pruning
    float tol[2] = {FB_TOL, TU_TOL};
    int p[2] = {1, 0};

    // Create partitions data structures
    std::vector<PartitionData*> p_data;
    std::vector<PartitionData*> p_data_r;

    if (params.split_img == 1) {
        auto clk_init_part = system_clock::now(); // PROFILING
        // Initialise partitions
        // Note: to avoid reinforcing the discontinuities that may be caused by the partitions, we
        // flip the grid/partition to avoid 'cutting' the image twice on successive iterations the same spot for two
        // (i.e.: 3x2 => 2x3)
        // Odd iterations (i == 1, 3, 5, ...)
        init_subimage_partitions(i0, i1, i_1, i2, i0n, i1n, i_1n, i2n, BiFilt_Go, BiFilt_Ba, sal_go, sal_ba, w, h,
                                 params.h_parts, params.v_parts, &p_data, params);

        // Even iterations (i == 2, 4, 6, ...) ==> v_parts x h_parts
        init_subimage_partitions(i0, i1, i_1, i2, i0n, i1n, i_1n, i2n, BiFilt_Go, BiFilt_Ba, sal_go, sal_ba, w, h,
                                 params.v_parts, params.h_parts, &p_data_r, params);

        auto clk_init_part_end = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_init_part = clk_init_part_end - clk_init_part;  // PROFILING
        cout << "(match growing) initialising partitions took "
             << elapsed_secs_init_part.count() << endl;
    }

    // Main local FALDOI loop (i.e.: 'iterated faldoi')
    for (int i = 0; i < iter; i++) {
        auto clk_init_iter = system_clock::now();  // PROFILING
        printf("Iteration: %d\n", i);

        // Estimate local minimization (forward: I0 ==> I1 and backward: I0 <== I1)
        // First iteration works on the whole image
        if (params.split_img == 0 || (params.split_img == 1 && i == 0)) {
#ifdef _OPENMP
#pragma omp parallel for
            for (int k = 0; k < 2; k++) {
                if (k == 0) {
                    // FWD
                    auto clk_start_fwd = system_clock::now();
                    local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go, BiFilt_Go,
                                  true, w, h);
                    auto clk_fwd_grow = system_clock::now();  // PROFILING
                    duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                    cout << "(match growing) FWD local growing (it=" << i << ") took "
                         << elapsed_secs_fwd_grow.count() << endl;
                } else {
                    // BWD
                    auto clk_start_bwd = system_clock::now();
                    local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                  false, w, h);
                    auto clk_bwd_grow = system_clock::now();  // PROFILING
                    duration<double> elapsed_secs_bwd_grow = clk_bwd_grow - clk_start_bwd;  // PROFILING
                    cout << "(match growing) BWD local growing (it=" << i << ") took "
                         << elapsed_secs_bwd_grow.count() << endl;
                }
            }
#endif

            auto clk_loc_grow_end = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_loc_grow = clk_loc_grow_end - clk_init_iter;  // PROFILING
            cout << "(match growing) FWD + BWD local growing (it=" << i << ") took "
                 << elapsed_secs_loc_grow.count() << endl;

            // Pruning method
            pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

            auto clk_pruning = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_prune = clk_pruning - clk_loc_grow_end; // PROFILING
            cout << "(match growing) Local iteration " << i << " => pruning method took "
                 << elapsed_secs_prune.count() << endl;

            // Delete not trustable candidates based on the previous pruning
            delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
            delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

            auto clk_delete_non_trust = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_delete = clk_delete_non_trust - clk_pruning; // PROFILING
            cout << "(match growing) Local iteration " << i << " => delete non-trustable candidates "
                 << elapsed_secs_delete.count() << endl;

            // Insert each pixel into the queue as possible candidate
            insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, occ_Go, w, h);
            insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, occ_Ba, w, h);

            auto clk_insert_cand = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_insert_cand = clk_insert_cand - clk_delete_non_trust; // PROFILING
            cout << "(match growing) Local iteration " << i << " => insert potential candidates "
                 << elapsed_secs_insert_cand.count() << endl;

            prepare_data_for_growing(&ofGo, ene_Go, oft0, w, h);
            prepare_data_for_growing(&ofBa, ene_Ba, oft1, w, h);

            auto clk_prepare_grow = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_prepare_grow = clk_prepare_grow - clk_insert_cand; // PROFILING
            cout << "(match growing) Local iteration " << i << " => prepare data for growing "
                 << elapsed_secs_prepare_grow.count() << endl;

            auto clk_all_tasks = system_clock::now(); // PROFILING
            duration<double> elapsed_secs_all_tasks = clk_all_tasks - clk_init_iter; // PROFILING
            cout << "(match growing) Local iteration " << i << " => all iteration's tasks took "
                 << elapsed_secs_all_tasks.count() << endl;
        }
        else if ((i > 0 && i <= iter - 1) && params.split_img == 1) {
            // Common stuff to any iteration from 2nd to last
            const int n_partitions = params.h_parts * params.v_parts;

            if (i % 2 != 0 && i <= iter - 1)  {  // part. grid: h_parts (cols) x v_parts (rows)
                auto clk_odd_start = system_clock::now();  // PROFILING
                // Update partition-specific variables with the image-specific values from the first iteration
                image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                                    queue_Go, queue_Ba, n_partitions, w, h, &p_data, true);

                auto clk_update_part = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_update_part = clk_update_part - clk_odd_start; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Update partitions (image => part) took "
                     << elapsed_secs_update_part.count() << endl;

                // If at least one seed falls into each partition queue, we can parallelise the code
                if (!anyEmptyQueues(&p_data, n_partitions)) {
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned n = 0; n < n_partitions * 2; n++) {
                        // 1. Local growing based on updated OF variables ofGo and ofBa (fwd/bwd)
                        if (n % 2 == 0) {
                            // FWD
                            unsigned m = n / 2;
                            std::cout << "(FWD) Local growing partition (h x v) => " << m << std::endl;
                            auto clk_start_fwd = system_clock::now();
                            local_growing(p_data.at(m)->i0n, p_data.at(m)->i1n, p_data.at(m)->i_1n,
                                          &(p_data.at(m)->queue_Go),
                                          &(p_data.at(m)->stuffGo), &(p_data.at(m)->ofGo), i, p_data.at(m)->ene_Go,
                                          p_data.at(m)->oft0, p_data.at(m)->occ_Go, p_data.at(m)->BiFilt_Go, true,
                                          p_data.at(m)->width, p_data.at(m)->height);

                            auto clk_fwd_grow = system_clock::now(); // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd; // PROFILING
                            cout << "(match growing) Local iteration " << i << ", partition " << m
                                 << " => FWD growing took "
                                 << elapsed_secs_fwd_grow.count() << endl;
                        } else {
                            // BWD
                            unsigned m = n / 2;  // taking advantage of integer division (otherwise we'd have wrote ==> m = (n-1)/2)
                            std::cout << "(BWD) Local growing partition (h x v) => " << m << std::endl;
                            auto clk_start_bwd = system_clock::now();
                            local_growing(p_data.at(m)->i1n, p_data.at(m)->i0n, p_data.at(m)->i2n,
                                          &(p_data.at(m)->queue_Ba),
                                          &(p_data.at(m)->stuffBa), &(p_data.at(m)->ofBa), i, p_data.at(m)->ene_Ba,
                                          p_data.at(m)->oft1, p_data.at(m)->occ_Ba, p_data.at(m)->BiFilt_Ba, false,
                                          p_data.at(m)->width, p_data.at(m)->height);

                            auto clk_bwd_grow = system_clock::now(); // PROFILING
                            duration<double> elapsed_secs_bwd_grow = clk_bwd_grow - clk_start_bwd; // PROFILING
                            cout << "(match growing) Local iteration " << i << ", partition " << m
                                 << " => BWD growing took "
                                 << elapsed_secs_bwd_grow.count() << endl;
                        }
                    }
#endif
                } else {
                    // Otherwise, we revert to using the whole image for the local growing
                    // Note that this is the quickest and most naive solution since we could only grow those that
                    // are not empty but that would imply signifcant changes to the code as is
                    cout << "Reverted back to whole-image based processing due to lack of seeds (1 or more empty queues)" << endl;
#ifdef _OPENMP
#pragma omp parallel for
                    for (int k = 0; k < 2; k++) {
                        if (k == 0) {
                            // FWD
                            auto clk_start_fwd = system_clock::now();
                            local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go, BiFilt_Go,
                                          true, w, h);
                            auto clk_fwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                            cout << "(match growing) FWD local growing (it=" << i << ") took "
                                 << elapsed_secs_fwd_grow.count() << endl;
                        } else {
                            // BWD
                            auto clk_start_bwd = system_clock::now();
                            local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                          false, w, h);
                            auto clk_bwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_bwd_grow = clk_bwd_grow - clk_start_bwd;  // PROFILING
                            cout << "(match growing) BWD local growing (it=" << i << ") took "
                                 << elapsed_secs_bwd_grow.count() << endl;

                        }
                    }
#endif
                }

                auto clk_grow = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_grow = clk_grow - clk_update_part; // PROFILING
                cout << "(match growing) Local iteration " << i <<" => All FWD + BWD growings took "
                     << elapsed_secs_grow.count() << endl;

                if (!anyEmptyQueues(&p_data, n_partitions)) {
                    // Copy partition growing information back to image-wise variables for pruning
                    image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                                        queue_Go, queue_Ba, n_partitions, w, h, &p_data, false);
                }
                auto clk_update_image = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_update_image = clk_update_image - clk_grow; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Update partitions (part => image) took "
                     << elapsed_secs_update_image.count() << endl;

                // 2. Pruning
                pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

                auto clk_pruning = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_pruning = clk_pruning - clk_update_image; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Pruning method took "
                     << elapsed_secs_pruning.count() << endl;

                // 3. Delete non-trustable
                delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
                delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

                auto clk_delete = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_delete = clk_delete - clk_pruning; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Deleting non-trustable candidates took "
                     << elapsed_secs_delete.count() << endl;

                // 4. insert potential candidates
                insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, occ_Go, w, h);
                insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, occ_Ba, w, h);

                auto clk_insert_cand = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_insert_cand = clk_insert_cand - clk_delete; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Inserting potential candidates took "
                     << elapsed_secs_insert_cand.count() << endl;

                // 5. prepare data for growing
                prepare_data_for_growing(&ofGo, ene_Go, oft0, w, h);
                prepare_data_for_growing(&ofBa, ene_Ba, oft1, w, h);

                auto clk_prepare_grow = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_prepare_grow = clk_prepare_grow - clk_insert_cand; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Preparing data for grwing took "
                     << elapsed_secs_prepare_grow.count() << endl;

                auto clk_all_tasks = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_all_tasks = clk_all_tasks - clk_odd_start; // PROFILING
                cout << "(match growing) Local iteration " << i << " => All iteration's tasks took "
                     << elapsed_secs_all_tasks.count() << endl;
            }
            else if (i % 2 == 0 && i <= iter - 1)        // part. grid: v_parts (cols) x h_parts (rows)
            {
                auto clk_even_start = system_clock::now(); // PROFILING

                image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                                    queue_Go, queue_Ba, n_partitions, w, h, &p_data_r, true);

                auto clk_update_part = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_update_part = clk_update_part - clk_even_start; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Update partitions (image => part) took "
                     << elapsed_secs_update_part.count() << endl;

                // If at least one seed falls into each partition queue, we can parallelise the code
                if (!anyEmptyQueues(&p_data_r, n_partitions)) {
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned n = 0; n < n_partitions * 2; n++) {
                        if (n % 2 == 0) {
                            // 1. Local growing based on updated OF variables, ofGo, ofBa (previous iteration)
                            // FWD
                            unsigned m = n / 2;
                            std::cout << "(FWD) Local growing partition (v x h) => " << m << std::endl;
                            auto clk_start_fwd = system_clock::now();
                            local_growing(p_data_r.at(m)->i0n, p_data_r.at(m)->i1n, p_data_r.at(m)->i_1n,
                                          &(p_data_r.at(m)->queue_Go), &(p_data_r.at(m)->stuffGo),
                                          &(p_data_r.at(m)->ofGo), i, p_data_r.at(m)->ene_Go, p_data_r.at(m)->oft0,
                                          p_data_r.at(m)->occ_Go, p_data_r.at(m)->BiFilt_Go, true,
                                          p_data_r.at(m)->width, p_data_r.at(m)->height);

                            auto clk_fwd_grow = system_clock::now(); // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd; // PROFILING
                            cout << "(match growing) Local iteration " << i << ", partition " << m
                                 << " => FWD growing took "
                                 << elapsed_secs_fwd_grow.count() << endl;

                        } else {
                            // BWD
                            unsigned m = n / 2;  // Valid because the integer division is applied (otherwise m = (n-1)/2 would be used)
                            std::cout << "(BWD) Local growing partition (v x h) => " << m << std::endl;
                            auto clk_start_bwd = system_clock::now();
                            local_growing(p_data_r.at(m)->i1n, p_data_r.at(m)->i0n, p_data_r.at(m)->i2n,
                                          &(p_data_r.at(m)->queue_Ba), &(p_data_r.at(m)->stuffBa),
                                          &(p_data_r.at(m)->ofBa), i, p_data_r.at(m)->ene_Ba, p_data_r.at(m)->oft1,
                                          p_data_r.at(m)->occ_Ba, p_data_r.at(m)->BiFilt_Ba, false,
                                          p_data_r.at(m)->width, p_data_r.at(m)->height);

                            auto clk_bwd_grow = system_clock::now(); // PROFILING
                            duration<double> elapsed_secs_bwd_grow = clk_bwd_grow - clk_start_bwd; // PROFILING
                            cout << "(match growing) Local iteration " << i << ", partition " << m
                                 << " => BWD growing took "
                                 << elapsed_secs_bwd_grow.count() << endl;
                        }
                    }
#endif
                } else {

                    // Otherwise, we revert to using the whole image for the local growing
                    cout << "Reverted back to whole-image based processing due to lack of seeds (1 or more empty queues)" << endl;
#ifdef _OPENMP
#pragma omp parallel for
                    for (int k = 0; k < 2; k++) {
                        if (k == 0) {
                            // FWD
                            auto clk_start_fwd = system_clock::now();
                            local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go,
                                          BiFilt_Go,
                                          true, w, h);
                            auto clk_fwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                            cout << "(match growing) FWD local growing (it=" << i << ") took "
                                 << elapsed_secs_fwd_grow.count() << endl;
                        } else {
                            // BWD
                            auto clk_start_bwd = system_clock::now();
                            local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                          false, w, h);
                            auto clk_bwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_bwd_grow = clk_bwd_grow - clk_start_bwd;  // PROFILING
                            cout << "(match growing) BWD local growing (it=" << i << ") took "
                                 << elapsed_secs_bwd_grow.count() << endl;
                        }
                    }
#endif
                }

                auto clk_grow = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_grow = clk_grow - clk_update_part; // PROFILING
                cout << "(match growing) Local iteration " << i <<" => All FWD + BWD growings took "
                     << elapsed_secs_grow.count() << endl;

                if (!anyEmptyQueues(&p_data_r, n_partitions)) {
                    // Copy partition growing information back to image-wise variables for pruning
                    image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                                        queue_Go, queue_Ba, n_partitions, w, h, &p_data_r, false);
                }

                auto clk_update_image = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_update_image = clk_update_image - clk_grow; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Update partitions (part => image) took "
                     << elapsed_secs_update_image.count() << endl;

                // 2. Pruning
                pruning_method(i0n, i1n, w, h, tol, p, ofGo.trust_points, oft0, ofBa.trust_points, oft1);

                auto clk_pruning = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_pruning = clk_pruning - clk_update_image; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Pruning method took "
                     << elapsed_secs_pruning.count() << endl;

                // 3. Delete non-trustable
                delete_not_trustable_candidates(&ofGo, oft0, ene_Go, w, h);
                delete_not_trustable_candidates(&ofBa, oft1, ene_Ba, w, h);

                auto clk_delete = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_delete = clk_delete - clk_pruning; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Deleting non-trustable candidates took "
                     << elapsed_secs_delete.count() << endl;

                // 4. insert potential candidates
                insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, occ_Go, w, h);
                insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, occ_Ba, w, h);

                auto clk_insert_cand = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_insert_cand = clk_insert_cand - clk_delete; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Inserting potential candidates took "
                     << elapsed_secs_insert_cand.count() << endl;

                // 5. prepare data for growing
                prepare_data_for_growing(&ofGo, ene_Go, oft0, w, h);
                prepare_data_for_growing(&ofBa, ene_Ba, oft1, w, h);

                auto clk_prepare_grow = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_prepare_grow = clk_prepare_grow - clk_insert_cand; // PROFILING
                cout << "(match growing) Local iteration " << i << " => Preparing data for grwing took "
                     << elapsed_secs_prepare_grow.count() << endl;

                auto clk_all_tasks = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_all_tasks = clk_all_tasks - clk_even_start; // PROFILING
                cout << "(match growing) Local iteration " << i << " => All iteration's tasks took "
                     << elapsed_secs_all_tasks.count() << endl;
            }
        }
    }

    printf("Last growing (FWD only)\n");

    if (params.split_img == 1) {
        const int n_partitions = params.h_parts * params.v_parts;

        auto clk_update_last_start = system_clock::now(); // PROFILING
        // NOTE: need to update partition (copy image-wise values to partition-specific variables)
        image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                            queue_Go, queue_Ba, n_partitions, w, h, &p_data, true);

        auto clk_update_last = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_update_part = clk_update_last - clk_update_last_start; // PROFILING
        cout << "(match growing) Last local iteration " << iter << " => Update partitions (image => part) took "
             << elapsed_secs_update_part.count() << endl;

        auto last_growing = system_clock::now();

        if (!anyEmptyQueues(&p_data, n_partitions)) {
#ifdef _OPENMP
#pragma omp parallel for
            for (unsigned m = 0; m < n_partitions; m++) {
                std::cout << "Last local growing partition (h x v) => " << m << std::endl;
                auto clk_start_fwd = system_clock::now();
                // FWD only (n_partitions function calls in parallel)
                local_growing(p_data.at(m)->i0n, p_data.at(m)->i1n, p_data.at(m)->i_1n, &(p_data.at(m)->queue_Go),
                              &(p_data.at(m)->stuffGo), &(p_data.at(m)->ofGo), iter, p_data.at(m)->ene_Go,
                              p_data.at(m)->oft0, p_data.at(m)->occ_Go, p_data.at(m)->BiFilt_Go, true,
                              p_data.at(m)->width, p_data.at(m)->height);

                auto clk_fwd_grow = system_clock::now(); // PROFILING
                duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd; // PROFILING
                cout << "(match growing) Last local iteration " << iter << ", partition " << m
                     << " => FWD growing took " << elapsed_secs_fwd_grow.count() << endl;
            }
#endif
        } else {
            // Revert to whole-image-based growing (slower but more robust with few seeds)
            // FWD only
            cout << "Reverted back to whole-image based processing due to lack of seeds (1 or more empty queues)" << endl;
            local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, iter, ene_Go, oft0, occ_Go, BiFilt_Go, true, w, h);
        }
        auto clk_end = system_clock::now(); // PROFILING
        duration<double> elapsed_secs = clk_end - last_growing; // PROFILING
        cout << "(match growing) Last growing (FWD only) took "
             << elapsed_secs.count() << endl;
        if (!anyEmptyQueues(&p_data, n_partitions)) {
            // Copy partition growing information back to image-wise variables for pruning
            image_to_partitions(oft0, oft1, ene_Go, ene_Ba, occ_Go, occ_Ba, &ofGo, &ofBa, &stuffGo, &stuffBa,
                                queue_Go, queue_Ba, n_partitions, w, h, &p_data, false);
        }

        auto clk_update_image = system_clock::now(); // PROFILING
        duration<double> elapsed_secs_update_image = clk_update_image - clk_end; // PROFILING
        cout << "(match growing) Last local iteration " << iter << " => Update partitions (part => image) took "
             << elapsed_secs_update_image.count() << endl;

    } else {
        auto last_growing = system_clock::now();  // PROFILING
        local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, iter, ene_Go, oft0, occ_Go, BiFilt_Go, true, w, h);

        auto clk_end = system_clock::now(); // PROFILING
        duration<double> elapsed_secs = clk_end - last_growing; // PROFILING
        cout << "(match growing) Last growing (FWD only) took "
             << elapsed_secs.count() << endl;
    }

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
        fprintf(stderr, "Without occlusions:\n");
        fprintf(stderr, "Usage: 2 images in the following order: I0, I1\n");
        fprintf(stderr, "With occlusions:\n");
        fprintf(stderr, "Usage: 4 images in the following order: I0, I1, I-1, I2\n");

        return 1;
    }

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
    // Now absolute paths and relative ones (w.r.t to current dir, i.e.: '../example_data/sample/image_001.png')

    // Frame t-1 and t+2
    float *i_1 = nullptr;
    float *i2 = nullptr;
    if (num_files == 4) {
        i_1 = iio_read_image_float_split(filename_i_1.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i2.c_str(), w + 7, h + 7, pd + 5);
    } else {
        i_1 = iio_read_image_float_split(filename_i0.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i1.c_str(), w + 7, h + 7, pd + 5);
    }
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

    for (int i = 0; i < w[0] * h[0]; i++) {
        sal0[i] = 1.0;
        sal1[i] = 1.0;
    }

    for (int i = 0; i < w[0] * h[0]; i++) {
        assert(isfinite(sal0[i]));
        assert(isfinite(sal1[i]));
        if (sal0[i] < 0.0)
            fprintf(stderr, "sal0 index: %d\n", i);
        if (sal1[i] < 0.0)
            fprintf(stderr, "sal1 index: %d\n", i);
    }

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

    auto clk1 = system_clock::now(); // PROFILING
    // Match growing algorithm
    match_growing_variational(go, ba, i0, i1, i_1, i2, sal0, sal1, params, ene_val, out_flow, out_occ);

    auto clk2 = system_clock::now(); // PROFILING
    duration<double> elapsed_secs2 = clk2 - clk1; // PROFILING
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
