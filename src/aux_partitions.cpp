// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2018, Ferran Pérez <fperez.gamonal@gmail.com>
// All rights reserved.

#ifndef AUX_PARTITIONS
#define AUX_PARTITIONS

#include <cassert>
#include <cmath>

#include "aux_partitions.h"
#include "energy_structures.h"
#include "energy_model.h"

////////////////////////////////////////////////////////////////////////////////
//////////////////LOCAL PARTITION INTO SUBIMAGES////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief               Initialises the subpartitions structures' fields with their default values.
 *
 * @param i0            source frame at time 't'
 * @param i1            second frame at time 't+1'
 * @param i_1           previous frame at time 't-1' (used for occlusions only)
 * @param i2            second frame at time 't+2' (used for occlusions only)
 * @param i0n           normalised (gray and smooth) source frame
 * @param i1n           normalised (gray and smooth) second frame
 * @param i_1n          normalised (gray and smooth) previous frame
 * @param i2n           normalised (gray and smooth) second frame
 * @param BiFilt_Go     bilateral filter values for forward flow
 * @param BiFilt_Ba     bilateral filter values for backward flow
 * @param sal_go        saliency values for forward flow
 * @param sal_ba        saliency values for backward flow
 * @param w_src         width of the input frames
 * @param h_src         height of the input frames
 * @param h_parts       number of horizontal parts of the partition grid (e.g.: a 3x2 grid has h_parts=3)
 * @param v_parts       number of vertical parts of the partition grid (e.g.: a 3x2 grid has h_parts=2)
 * @param p_data        vector of data structs (one per partition) that define all needed parameters (see PartitionData)
 * @param params        struct of basic params initialised earlier in this source file's main function
 *
 * @sa                  image_to_partitions
 */
void  init_subimage_partitions(const float *i0, const float *i1, const float *i_1, const float *i2, const float *i0n,
                               const float *i1n, const float *i_1n, const float *i2n, BilateralFilterData* BiFilt_Go,
                               BilateralFilterData* BiFilt_Ba, float *sal_go, float *sal_ba, const int w_src,
                               const int h_src, const int h_parts, const int v_parts,
                               std::vector<PartitionData*> *p_data, Parameters params)
{
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
        auto *t_pdata = new PartitionData;  // Create new struct of partition data and initialise fields
        t_pdata->idx = p;
        t_pdata->width = sub_w[p];
        t_pdata->height = sub_h[p];
        t_pdata->off_x = off_x[p];
        t_pdata->off_y = off_y[p];
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
                    int m = j * p_data->at(p)->width + i + k * p_data->at(p)->width * p_data->at(p)->height;
                    // idx of the subarray is computed as follows (from image array):
                    // where:   off_x, off_y (subimage offsets, ref.: top-left corner)
                    //          i, j indices (sub_width[p] cols, sub_height[p] rows)
                    //          k index for the channel (always 0 for '1d' variables, '0' or '1' for 2d (i.e.: the flow)
                    //          p is the partition index
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
    }

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
                    continue;  // Check following iteration
                }
            }
        // Assign to p_data
        p_data->at(p)->BiFilt_Go = Filter_data;
    }

    // BiFilt_Ba (bwd)
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
                    continue;  // Check following iteration
                }
            }
        p_data->at(p)->BiFilt_Ba = Filter_data;
    }
}


/**
 * @brief                   update OpticalFlowData struct when using partitions (common to any functional)
 * @details                 copies the values from image-wise variables to partition-wise (or viceversa)
 *
 * @param ofGo              forward OpticalFlowData struct
 * @param ofBa              backward OpticalFlowData struct
 * @param p_data            vector of data structs (one/partition) that define all needed parameters (see PartitionData)
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param n_ch              index of the channel (0 for all variables but the ones containing the flow field or similar)
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 */
void update_of_data(OpticalFlowData *& ofGo, OpticalFlowData *& ofBa, PartitionData *& p_data, const int idx_img,
                    const int idx_par, const int n_ch, bool img_to_part)
{
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
            p_data->ofBa.u1[idx_par] = ofBa->u1[idx_img];

            p_data->ofGo.u1_ba[idx_par] = ofGo->u1_ba[idx_img];
            p_data->ofBa.u1_ba[idx_par] = ofBa->u1_ba[idx_img];

            // Filters
            p_data->ofGo.u1_filter[idx_par] = ofGo->u1_filter[idx_img];
            p_data->ofBa.u1_filter[idx_par] = ofBa->u1_filter[idx_img];
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
            ofBa->u1[idx_img] = p_data->ofBa.u1[idx_par];

            ofGo->u1_ba[idx_img] = p_data->ofGo.u1_ba[idx_par];
            ofBa->u1_ba[idx_img] = p_data->ofBa.u1_ba[idx_par];

            // Filters
            ofGo->u1_filter[idx_img] = p_data->ofGo.u1_filter[idx_par];
            ofBa->u1_filter[idx_img] = p_data->ofBa.u1_filter[idx_par];
        }
    }
}

/// Set of functions to update each SpecificOFStuff struct (functional-dependant) ///

/**
 * @brief                   updates SpecificOFStuff struct for the TVL1 functional
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_tvl2_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                         const int idx_img, const int idx_par, bool img_to_part)
{
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

    } else {
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


/**
 * @brief                   updates SpecificOFStuff struct for the TVL1 functional with weights
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_tvl1w_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                          const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Xi
        p_data->stuffGo.tvl2w.xi11[idx_par] = stuffGo->tvl2w.xi11[idx_img];
        p_data->stuffGo.tvl2w.xi12[idx_par] = stuffGo->tvl2w.xi12[idx_img];
        p_data->stuffGo.tvl2w.xi21[idx_par] = stuffGo->tvl2w.xi21[idx_img];
        p_data->stuffGo.tvl2w.xi22[idx_par] = stuffGo->tvl2w.xi22[idx_img];

        p_data->stuffBa.tvl2w.xi11[idx_par] = stuffBa->tvl2w.xi11[idx_img];
        p_data->stuffBa.tvl2w.xi12[idx_par] = stuffBa->tvl2w.xi12[idx_img];
        p_data->stuffBa.tvl2w.xi21[idx_par] = stuffBa->tvl2w.xi21[idx_img];
        p_data->stuffBa.tvl2w.xi22[idx_par] = stuffBa->tvl2w.xi22[idx_img];

        // u1, u2
        p_data->stuffGo.tvl2w.u1x[idx_par] = stuffGo->tvl2w.u1x[idx_img];
        p_data->stuffGo.tvl2w.u1y[idx_par] = stuffGo->tvl2w.u1y[idx_img];
        p_data->stuffGo.tvl2w.u2x[idx_par] = stuffGo->tvl2w.u2x[idx_img];
        p_data->stuffGo.tvl2w.u2y[idx_par] = stuffGo->tvl2w.u2y[idx_img];

        p_data->stuffBa.tvl2w.u1x[idx_par] = stuffBa->tvl2w.u1x[idx_img];
        p_data->stuffBa.tvl2w.u1y[idx_par] = stuffBa->tvl2w.u1y[idx_img];
        p_data->stuffBa.tvl2w.u2x[idx_par] = stuffBa->tvl2w.u2x[idx_img];
        p_data->stuffBa.tvl2w.u2y[idx_par] = stuffBa->tvl2w.u2y[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.tvl2w.v1[idx_par] = stuffGo->tvl2w.v1[idx_img];
        p_data->stuffGo.tvl2w.v2[idx_par] = stuffGo->tvl2w.v2[idx_img];

        p_data->stuffBa.tvl2w.v1[idx_par] = stuffBa->tvl2w.v1[idx_img];
        p_data->stuffBa.tvl2w.v2[idx_par] = stuffBa->tvl2w.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.tvl2w.rho_c[idx_par] = stuffGo->tvl2w.rho_c[idx_img];
        p_data->stuffBa.tvl2w.rho_c[idx_par] = stuffBa->tvl2w.rho_c[idx_img];

        p_data->stuffGo.tvl2w.grad[idx_par] = stuffGo->tvl2w.grad[idx_img];
        p_data->stuffBa.tvl2w.grad[idx_par] = stuffBa->tvl2w.grad[idx_img];

        p_data->stuffGo.tvl2w.u1_[idx_par] = stuffGo->tvl2w.u1_[idx_img];
        p_data->stuffGo.tvl2w.u2_[idx_par] = stuffGo->tvl2w.u2_[idx_img];

        p_data->stuffBa.tvl2w.u1_[idx_par] = stuffBa->tvl2w.u1_[idx_img];
        p_data->stuffBa.tvl2w.u2_[idx_par] = stuffBa->tvl2w.u2_[idx_img];

        p_data->stuffGo.tvl2w.u1Aux[idx_par] = stuffGo->tvl2w.u1Aux[idx_img];
        p_data->stuffGo.tvl2w.u2Aux[idx_par] = stuffGo->tvl2w.u2Aux[idx_img];

        p_data->stuffBa.tvl2w.u1Aux[idx_par] = stuffBa->tvl2w.u1Aux[idx_img];
        p_data->stuffBa.tvl2w.u2Aux[idx_par] = stuffBa->tvl2w.u2Aux[idx_img];

        p_data->stuffGo.tvl2w.I1x[idx_par] = stuffGo->tvl2w.I1x[idx_img];
        p_data->stuffGo.tvl2w.I1y[idx_par] = stuffGo->tvl2w.I1y[idx_img];

        p_data->stuffBa.tvl2w.I1x[idx_par] = stuffBa->tvl2w.I1x[idx_img];
        p_data->stuffBa.tvl2w.I1y[idx_par] = stuffBa->tvl2w.I1y[idx_img];

        p_data->stuffGo.tvl2w.I1wx[idx_par] = stuffGo->tvl2w.I1wx[idx_img];
        p_data->stuffGo.tvl2w.I1wy[idx_par] = stuffGo->tvl2w.I1wy[idx_img];

        p_data->stuffBa.tvl2w.I1wx[idx_par] = stuffBa->tvl2w.I1wx[idx_img];
        p_data->stuffBa.tvl2w.I1wy[idx_par] = stuffBa->tvl2w.I1wy[idx_img];

        p_data->stuffGo.tvl2w.div_xi1[idx_par] = stuffGo->tvl2w.div_xi1[idx_img];
        p_data->stuffGo.tvl2w.div_xi2[idx_par] = stuffGo->tvl2w.div_xi2[idx_img];

        p_data->stuffBa.tvl2w.div_xi1[idx_par] = stuffBa->tvl2w.div_xi1[idx_img];
        p_data->stuffBa.tvl2w.div_xi2[idx_par] = stuffBa->tvl2w.div_xi2[idx_img];

        p_data->stuffGo.tvl2w.u_N[idx_par] = stuffGo->tvl2w.u_N[idx_img];
        p_data->stuffBa.tvl2w.u_N[idx_par] = stuffBa->tvl2w.u_N[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // Xi
        stuffGo->tvl2w.xi11[idx_img] = p_data->stuffGo.tvl2w.xi11[idx_par];
        stuffGo->tvl2w.xi12[idx_img] = p_data->stuffGo.tvl2w.xi12[idx_par];
        stuffGo->tvl2w.xi21[idx_img] = p_data->stuffGo.tvl2w.xi21[idx_par];
        stuffGo->tvl2w.xi22[idx_img] = p_data->stuffGo.tvl2w.xi22[idx_par];

        stuffBa->tvl2w.xi11[idx_img] = p_data->stuffBa.tvl2w.xi11[idx_par];
        stuffBa->tvl2w.xi12[idx_img] = p_data->stuffBa.tvl2w.xi12[idx_par];
        stuffBa->tvl2w.xi21[idx_img] = p_data->stuffBa.tvl2w.xi21[idx_par];
        stuffBa->tvl2w.xi22[idx_img] = p_data->stuffBa.tvl2w.xi22[idx_par];

        // u1, u2
        stuffGo->tvl2w.u1x[idx_img] = p_data->stuffGo.tvl2w.u1x[idx_par];
        stuffGo->tvl2w.u1y[idx_img] = p_data->stuffGo.tvl2w.u1y[idx_par];
        stuffGo->tvl2w.u2x[idx_img] = p_data->stuffGo.tvl2w.u2x[idx_par];
        stuffGo->tvl2w.u2y[idx_img] = p_data->stuffGo.tvl2w.u2y[idx_par];

        stuffBa->tvl2w.u1x[idx_img] = p_data->stuffBa.tvl2w.u1x[idx_par];
        stuffBa->tvl2w.u1y[idx_img] = p_data->stuffBa.tvl2w.u1y[idx_par];
        stuffBa->tvl2w.u2x[idx_img] = p_data->stuffBa.tvl2w.u2x[idx_par];
        stuffBa->tvl2w.u2y[idx_img] = p_data->stuffBa.tvl2w.u2y[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->tvl2w.v1[idx_img] = p_data->stuffGo.tvl2w.v1[idx_par];
        stuffGo->tvl2w.v2[idx_img] = p_data->stuffGo.tvl2w.v2[idx_par];

        stuffBa->tvl2w.v1[idx_img] = p_data->stuffBa.tvl2w.v1[idx_par];
        stuffBa->tvl2w.v2[idx_img] = p_data->stuffBa.tvl2w.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->tvl2w.rho_c[idx_img] = p_data->stuffGo.tvl2w.rho_c[idx_par];
        stuffBa->tvl2w.rho_c[idx_img] = p_data->stuffBa.tvl2w.rho_c[idx_par];

        stuffGo->tvl2w.grad[idx_img] = p_data->stuffGo.tvl2w.grad[idx_par];
        stuffBa->tvl2w.grad[idx_img] = p_data->stuffBa.tvl2w.grad[idx_par];

        stuffGo->tvl2w.u1_[idx_img] = p_data->stuffGo.tvl2w.u1_[idx_par];
        stuffGo->tvl2w.u2_[idx_img] = p_data->stuffGo.tvl2w.u2_[idx_par];

        stuffBa->tvl2w.u1_[idx_img] = p_data->stuffBa.tvl2w.u1_[idx_par];
        stuffBa->tvl2w.u2_[idx_img] = p_data->stuffBa.tvl2w.u2_[idx_par];

        stuffGo->tvl2w.u1Aux[idx_img] = p_data->stuffGo.tvl2w.u1Aux[idx_par];
        stuffGo->tvl2w.u2Aux[idx_img] = p_data->stuffGo.tvl2w.u2Aux[idx_par];

        stuffBa->tvl2w.u1Aux[idx_img] = p_data->stuffBa.tvl2w.u1Aux[idx_par];
        stuffBa->tvl2w.u2Aux[idx_img] = p_data->stuffBa.tvl2w.u2Aux[idx_par];

        stuffGo->tvl2w.I1x[idx_img] = p_data->stuffGo.tvl2w.I1x[idx_par];
        stuffGo->tvl2w.I1y[idx_img] = p_data->stuffGo.tvl2w.I1y[idx_par];

        stuffBa->tvl2w.I1x[idx_img] = p_data->stuffBa.tvl2w.I1x[idx_par];
        stuffBa->tvl2w.I1y[idx_img] = p_data->stuffBa.tvl2w.I1y[idx_par];

        stuffGo->tvl2w.I1wx[idx_img] = p_data->stuffGo.tvl2w.I1wx[idx_par];
        stuffGo->tvl2w.I1wy[idx_img] = p_data->stuffGo.tvl2w.I1wy[idx_par];

        stuffBa->tvl2w.I1wx[idx_img] = p_data->stuffBa.tvl2w.I1wx[idx_par];
        stuffBa->tvl2w.I1wy[idx_img] = p_data->stuffBa.tvl2w.I1wy[idx_par];

        stuffGo->tvl2w.div_xi1[idx_img] = p_data->stuffGo.tvl2w.div_xi1[idx_par];
        stuffGo->tvl2w.div_xi2[idx_img] = p_data->stuffGo.tvl2w.div_xi2[idx_par];

        stuffBa->tvl2w.div_xi1[idx_img] = p_data->stuffBa.tvl2w.div_xi1[idx_par];
        stuffBa->tvl2w.div_xi2[idx_img] = p_data->stuffBa.tvl2w.div_xi2[idx_par];

        stuffGo->tvl2w.u_N[idx_img] = p_data->stuffGo.tvl2w.u_N[idx_par];
        stuffBa->tvl2w.u_N[idx_img] = p_data->stuffBa.tvl2w.u_N[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the NL-TVL1 functional
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_nltvl1_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                           const int idx_img, const int idx_par, bool img_to_part)
{
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

        p_data->stuffGo.nltvl1.div_p[idx_par] = stuffGo->nltvl1.div_p[idx_img];
        p_data->stuffGo.nltvl1.div_q[idx_par] = stuffGo->nltvl1.div_q[idx_img];

        p_data->stuffBa.nltvl1.div_p[idx_par] = stuffBa->nltvl1.div_p[idx_img];
        p_data->stuffBa.nltvl1.div_q[idx_par] = stuffBa->nltvl1.div_q[idx_img];

    } else {
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


/**
 * @brief                   updates SpecificOFStuff struct for the NL-TVL1 functional with weights
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_nltvl1w_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                            const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Dual variables
        p_data->stuffGo.nltvl1w.p[idx_par] = stuffGo->nltvl1w.p[idx_img];
        p_data->stuffGo.nltvl1w.q[idx_par] = stuffGo->nltvl1w.q[idx_img];

        p_data->stuffBa.nltvl1w.p[idx_par] = stuffBa->nltvl1w.p[idx_img];
        p_data->stuffBa.nltvl1w.q[idx_par] = stuffBa->nltvl1w.q[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.nltvl1w.v1[idx_par] = stuffGo->nltvl1w.v1[idx_img];
        p_data->stuffGo.nltvl1w.v2[idx_par] = stuffGo->nltvl1w.v2[idx_img];

        p_data->stuffBa.nltvl1w.v1[idx_par] = stuffBa->nltvl1w.v1[idx_img];
        p_data->stuffBa.nltvl1w.v2[idx_par] = stuffBa->nltvl1w.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.nltvl1w.rho_c[idx_par] = stuffGo->nltvl1w.rho_c[idx_img];
        p_data->stuffBa.nltvl1w.rho_c[idx_par] = stuffBa->nltvl1w.rho_c[idx_img];

        p_data->stuffGo.nltvl1w.grad[idx_par] = stuffGo->nltvl1w.grad[idx_img];
        p_data->stuffBa.nltvl1w.grad[idx_par] = stuffBa->nltvl1w.grad[idx_img];

        p_data->stuffGo.nltvl1w.u1_[idx_par] = stuffGo->nltvl1w.u1_[idx_img];
        p_data->stuffGo.nltvl1w.u2_[idx_par] = stuffGo->nltvl1w.u2_[idx_img];

        p_data->stuffBa.nltvl1w.u1_[idx_par] = stuffBa->nltvl1w.u1_[idx_img];
        p_data->stuffBa.nltvl1w.u2_[idx_par] = stuffBa->nltvl1w.u2_[idx_img];

        p_data->stuffGo.nltvl1w.u1_tmp[idx_par] = stuffGo->nltvl1w.u1_tmp[idx_img];
        p_data->stuffGo.nltvl1w.u2_tmp[idx_par] = stuffGo->nltvl1w.u2_tmp[idx_img];

        p_data->stuffBa.nltvl1w.u1_tmp[idx_par] = stuffBa->nltvl1w.u1_tmp[idx_img];
        p_data->stuffBa.nltvl1w.u2_tmp[idx_par] = stuffBa->nltvl1w.u2_tmp[idx_img];

        p_data->stuffGo.nltvl1w.I1x[idx_par] = stuffGo->nltvl1w.I1x[idx_img];
        p_data->stuffGo.nltvl1w.I1y[idx_par] = stuffGo->nltvl1w.I1y[idx_img];

        p_data->stuffBa.nltvl1w.I1x[idx_par] = stuffBa->nltvl1w.I1x[idx_img];
        p_data->stuffBa.nltvl1w.I1y[idx_par] = stuffBa->nltvl1w.I1y[idx_img];

        p_data->stuffGo.nltvl1w.I1wx[idx_par] = stuffGo->nltvl1w.I1wx[idx_img];
        p_data->stuffGo.nltvl1w.I1wy[idx_par] = stuffGo->nltvl1w.I1wy[idx_img];

        p_data->stuffBa.nltvl1w.I1wx[idx_par] = stuffBa->nltvl1w.I1wx[idx_img];
        p_data->stuffBa.nltvl1w.I1wy[idx_par] = stuffBa->nltvl1w.I1wy[idx_img];

        p_data->stuffGo.nltvl1w.div_p[idx_par] = stuffGo->nltvl1w.div_p[idx_img];
        p_data->stuffGo.nltvl1w.div_q[idx_par] = stuffGo->nltvl1w.div_q[idx_img];

        p_data->stuffBa.nltvl1w.div_p[idx_par] = stuffBa->nltvl1w.div_p[idx_img];
        p_data->stuffBa.nltvl1w.div_q[idx_par] = stuffBa->nltvl1w.div_q[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // Dual variables
        stuffGo->nltvl1w.p[idx_img] = p_data->stuffGo.nltvl1w.p[idx_par];
        stuffGo->nltvl1w.q[idx_img] = p_data->stuffGo.nltvl1w.q[idx_par];

        stuffBa->nltvl1w.p[idx_img] = p_data->stuffBa.nltvl1w.p[idx_par];
        stuffBa->nltvl1w.q[idx_img] = p_data->stuffBa.nltvl1w.q[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->nltvl1w.v1[idx_img] = p_data->stuffGo.nltvl1w.v1[idx_par];
        stuffGo->nltvl1w.v2[idx_img] = p_data->stuffGo.nltvl1w.v2[idx_par];

        stuffBa->nltvl1w.v1[idx_img] = p_data->stuffBa.nltvl1w.v1[idx_par];
        stuffBa->nltvl1w.v2[idx_img] = p_data->stuffBa.nltvl1w.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->nltvl1w.rho_c[idx_img] = p_data->stuffGo.nltvl1w.rho_c[idx_par];
        stuffBa->nltvl1w.rho_c[idx_img] = p_data->stuffBa.nltvl1w.rho_c[idx_par];

        stuffGo->nltvl1w.grad[idx_img] = p_data->stuffGo.nltvl1w.grad[idx_par];
        stuffBa->nltvl1w.grad[idx_img] = p_data->stuffBa.nltvl1w.grad[idx_par];

        stuffGo->nltvl1w.u1_[idx_img] = p_data->stuffGo.nltvl1w.u1_[idx_par];
        stuffGo->nltvl1w.u2_[idx_img] = p_data->stuffGo.nltvl1w.u2_[idx_par];

        stuffBa->nltvl1w.u1_[idx_img] = p_data->stuffBa.nltvl1w.u1_[idx_par];
        stuffBa->nltvl1w.u2_[idx_img] = p_data->stuffBa.nltvl1w.u2_[idx_par];

        stuffGo->nltvl1w.u1_tmp[idx_img] = p_data->stuffGo.nltvl1w.u1_tmp[idx_par];
        stuffGo->nltvl1w.u2_tmp[idx_img] = p_data->stuffGo.nltvl1w.u2_tmp[idx_par];

        stuffBa->nltvl1w.u1_tmp[idx_img] = p_data->stuffBa.nltvl1w.u1_tmp[idx_par];
        stuffBa->nltvl1w.u2_tmp[idx_img] = p_data->stuffBa.nltvl1w.u2_tmp[idx_par];

        stuffGo->nltvl1w.I1x[idx_img] = p_data->stuffGo.nltvl1w.I1x[idx_par];
        stuffGo->nltvl1w.I1y[idx_img] = p_data->stuffGo.nltvl1w.I1y[idx_par];

        stuffBa->nltvl1w.I1x[idx_img] = p_data->stuffBa.nltvl1w.I1x[idx_par];
        stuffBa->nltvl1w.I1y[idx_img] = p_data->stuffBa.nltvl1w.I1y[idx_par];

        stuffGo->nltvl1w.I1wx[idx_img] = p_data->stuffGo.nltvl1w.I1wx[idx_par];
        stuffGo->nltvl1w.I1wy[idx_img] = p_data->stuffGo.nltvl1w.I1wy[idx_par];

        stuffBa->nltvl1w.I1wx[idx_img] = p_data->stuffBa.nltvl1w.I1wx[idx_par];
        stuffBa->nltvl1w.I1wy[idx_img] = p_data->stuffBa.nltvl1w.I1wy[idx_par];

        stuffGo->nltvl1w.div_p[idx_img] = p_data->stuffGo.nltvl1w.div_p[idx_par];
        stuffGo->nltvl1w.div_q[idx_img] = p_data->stuffGo.nltvl1w.div_q[idx_par];

        stuffBa->nltvl1w.div_p[idx_img] = p_data->stuffBa.nltvl1w.div_p[idx_par];
        stuffBa->nltvl1w.div_q[idx_img] = p_data->stuffBa.nltvl1w.div_q[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the TV-CSAD functional
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_tvcsad_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                           const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // PosNei (neighbours' position)
        p_data->stuffGo.tvcsad.pnei[idx_par] = stuffGo->tvcsad.pnei[idx_img];
        p_data->stuffBa.tvcsad.pnei[idx_par] = stuffBa->tvcsad.pnei[idx_img];
        // Xi
        p_data->stuffGo.tvcsad.xi11[idx_par] = stuffGo->tvcsad.xi11[idx_img];
        p_data->stuffGo.tvcsad.xi12[idx_par] = stuffGo->tvcsad.xi12[idx_img];
        p_data->stuffGo.tvcsad.xi21[idx_par] = stuffGo->tvcsad.xi21[idx_img];
        p_data->stuffGo.tvcsad.xi22[idx_par] = stuffGo->tvcsad.xi22[idx_img];

        p_data->stuffBa.tvcsad.xi11[idx_par] = stuffBa->tvcsad.xi11[idx_img];
        p_data->stuffBa.tvcsad.xi12[idx_par] = stuffBa->tvcsad.xi12[idx_img];
        p_data->stuffBa.tvcsad.xi21[idx_par] = stuffBa->tvcsad.xi21[idx_img];
        p_data->stuffBa.tvcsad.xi22[idx_par] = stuffBa->tvcsad.xi22[idx_img];

        // u1, u2
        p_data->stuffGo.tvcsad.u1x[idx_par] = stuffGo->tvcsad.u1x[idx_img];
        p_data->stuffGo.tvcsad.u1y[idx_par] = stuffGo->tvcsad.u1y[idx_img];
        p_data->stuffGo.tvcsad.u2x[idx_par] = stuffGo->tvcsad.u2x[idx_img];
        p_data->stuffGo.tvcsad.u2y[idx_par] = stuffGo->tvcsad.u2y[idx_img];

        p_data->stuffBa.tvcsad.u1x[idx_par] = stuffBa->tvcsad.u1x[idx_img];
        p_data->stuffBa.tvcsad.u1y[idx_par] = stuffBa->tvcsad.u1y[idx_img];
        p_data->stuffBa.tvcsad.u2x[idx_par] = stuffBa->tvcsad.u2x[idx_img];
        p_data->stuffBa.tvcsad.u2y[idx_par] = stuffBa->tvcsad.u2y[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.tvcsad.v1[idx_par] = stuffGo->tvcsad.v1[idx_img];
        p_data->stuffGo.tvcsad.v2[idx_par] = stuffGo->tvcsad.v2[idx_img];

        p_data->stuffBa.tvcsad.v1[idx_par] = stuffBa->tvcsad.v1[idx_img];
        p_data->stuffBa.tvcsad.v2[idx_par] = stuffBa->tvcsad.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.tvcsad.rho_c[idx_par] = stuffGo->tvcsad.rho_c[idx_img];
        p_data->stuffBa.tvcsad.rho_c[idx_par] = stuffBa->tvcsad.rho_c[idx_img];

        p_data->stuffGo.tvcsad.grad[idx_par] = stuffGo->tvcsad.grad[idx_img];
        p_data->stuffBa.tvcsad.grad[idx_par] = stuffBa->tvcsad.grad[idx_img];

        p_data->stuffGo.tvcsad.u1_[idx_par] = stuffGo->tvcsad.u1_[idx_img];
        p_data->stuffGo.tvcsad.u2_[idx_par] = stuffGo->tvcsad.u2_[idx_img];

        p_data->stuffBa.tvcsad.u1_[idx_par] = stuffBa->tvcsad.u1_[idx_img];
        p_data->stuffBa.tvcsad.u2_[idx_par] = stuffBa->tvcsad.u2_[idx_img];

        p_data->stuffGo.tvcsad.u1_tmp[idx_par] = stuffGo->tvcsad.u1_tmp[idx_img];
        p_data->stuffGo.tvcsad.u2_tmp[idx_par] = stuffGo->tvcsad.u2_tmp[idx_img];

        p_data->stuffBa.tvcsad.u1_tmp[idx_par] = stuffBa->tvcsad.u1_tmp[idx_img];
        p_data->stuffBa.tvcsad.u2_tmp[idx_par] = stuffBa->tvcsad.u2_tmp[idx_img];

        p_data->stuffGo.tvcsad.I1x[idx_par] = stuffGo->tvcsad.I1x[idx_img];
        p_data->stuffGo.tvcsad.I1y[idx_par] = stuffGo->tvcsad.I1y[idx_img];

        p_data->stuffBa.tvcsad.I1x[idx_par] = stuffBa->tvcsad.I1x[idx_img];
        p_data->stuffBa.tvcsad.I1y[idx_par] = stuffBa->tvcsad.I1y[idx_img];

        p_data->stuffGo.tvcsad.I1wx[idx_par] = stuffGo->tvcsad.I1wx[idx_img];
        p_data->stuffGo.tvcsad.I1wy[idx_par] = stuffGo->tvcsad.I1wy[idx_img];

        p_data->stuffBa.tvcsad.I1wx[idx_par] = stuffBa->tvcsad.I1wx[idx_img];
        p_data->stuffBa.tvcsad.I1wy[idx_par] = stuffBa->tvcsad.I1wy[idx_img];

        p_data->stuffGo.tvcsad.div_xi1[idx_par] = stuffGo->tvcsad.div_xi1[idx_img];
        p_data->stuffGo.tvcsad.div_xi2[idx_par] = stuffGo->tvcsad.div_xi2[idx_img];

        p_data->stuffBa.tvcsad.div_xi1[idx_par] = stuffBa->tvcsad.div_xi1[idx_img];
        p_data->stuffBa.tvcsad.div_xi2[idx_par] = stuffBa->tvcsad.div_xi2[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // PosNei
        stuffGo->tvcsad.pnei[idx_img] = p_data->stuffGo.tvcsad.pnei[idx_par];
        stuffBa->tvcsad.pnei[idx_img] = p_data->stuffBa.tvcsad.pnei[idx_par];

        // Xi
        stuffGo->tvcsad.xi11[idx_img] = p_data->stuffGo.tvcsad.xi11[idx_par];
        stuffGo->tvcsad.xi12[idx_img] = p_data->stuffGo.tvcsad.xi12[idx_par];
        stuffGo->tvcsad.xi21[idx_img] = p_data->stuffGo.tvcsad.xi21[idx_par];
        stuffGo->tvcsad.xi22[idx_img] = p_data->stuffGo.tvcsad.xi22[idx_par];

        stuffBa->tvcsad.xi11[idx_img] = p_data->stuffBa.tvcsad.xi11[idx_par];
        stuffBa->tvcsad.xi12[idx_img] = p_data->stuffBa.tvcsad.xi12[idx_par];
        stuffBa->tvcsad.xi21[idx_img] = p_data->stuffBa.tvcsad.xi21[idx_par];
        stuffBa->tvcsad.xi22[idx_img] = p_data->stuffBa.tvcsad.xi22[idx_par];

        // u1, u2
        stuffGo->tvcsad.u1x[idx_img] = p_data->stuffGo.tvcsad.u1x[idx_par];
        stuffGo->tvcsad.u1y[idx_img] = p_data->stuffGo.tvcsad.u1y[idx_par];
        stuffGo->tvcsad.u2x[idx_img] = p_data->stuffGo.tvcsad.u2x[idx_par];
        stuffGo->tvcsad.u2y[idx_img] = p_data->stuffGo.tvcsad.u2y[idx_par];

        stuffBa->tvcsad.u1x[idx_img] = p_data->stuffBa.tvcsad.u1x[idx_par];
        stuffBa->tvcsad.u1y[idx_img] = p_data->stuffBa.tvcsad.u1y[idx_par];
        stuffBa->tvcsad.u2x[idx_img] = p_data->stuffBa.tvcsad.u2x[idx_par];
        stuffBa->tvcsad.u2y[idx_img] = p_data->stuffBa.tvcsad.u2y[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->tvcsad.v1[idx_img] = p_data->stuffGo.tvcsad.v1[idx_par];
        stuffGo->tvcsad.v2[idx_img] = p_data->stuffGo.tvcsad.v2[idx_par];

        stuffBa->tvcsad.v1[idx_img] = p_data->stuffBa.tvcsad.v1[idx_par];
        stuffBa->tvcsad.v2[idx_img] = p_data->stuffBa.tvcsad.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->tvcsad.rho_c[idx_img] = p_data->stuffGo.tvcsad.rho_c[idx_par];
        stuffBa->tvcsad.rho_c[idx_img] = p_data->stuffBa.tvcsad.rho_c[idx_par];

        stuffGo->tvcsad.grad[idx_img] = p_data->stuffGo.tvcsad.grad[idx_par];
        stuffBa->tvcsad.grad[idx_img] = p_data->stuffBa.tvcsad.grad[idx_par];

        stuffGo->tvcsad.u1_[idx_img] = p_data->stuffGo.tvcsad.u1_[idx_par];
        stuffGo->tvcsad.u2_[idx_img] = p_data->stuffGo.tvcsad.u2_[idx_par];

        stuffBa->tvcsad.u1_[idx_img] = p_data->stuffBa.tvcsad.u1_[idx_par];
        stuffBa->tvcsad.u2_[idx_img] = p_data->stuffBa.tvcsad.u2_[idx_par];

        stuffGo->tvcsad.u1_tmp[idx_img] = p_data->stuffGo.tvcsad.u1_tmp[idx_par];
        stuffGo->tvcsad.u2_tmp[idx_img] = p_data->stuffGo.tvcsad.u2_tmp[idx_par];

        stuffBa->tvcsad.u1_tmp[idx_img] = p_data->stuffBa.tvcsad.u1_tmp[idx_par];
        stuffBa->tvcsad.u2_tmp[idx_img] = p_data->stuffBa.tvcsad.u2_tmp[idx_par];

        stuffGo->tvcsad.I1x[idx_img] = p_data->stuffGo.tvcsad.I1x[idx_par];
        stuffGo->tvcsad.I1y[idx_img] = p_data->stuffGo.tvcsad.I1y[idx_par];

        stuffBa->tvcsad.I1x[idx_img] = p_data->stuffBa.tvcsad.I1x[idx_par];
        stuffBa->tvcsad.I1y[idx_img] = p_data->stuffBa.tvcsad.I1y[idx_par];

        stuffGo->tvcsad.I1wx[idx_img] = p_data->stuffGo.tvcsad.I1wx[idx_par];
        stuffGo->tvcsad.I1wy[idx_img] = p_data->stuffGo.tvcsad.I1wy[idx_par];

        stuffBa->tvcsad.I1wx[idx_img] = p_data->stuffBa.tvcsad.I1wx[idx_par];
        stuffBa->tvcsad.I1wy[idx_img] = p_data->stuffBa.tvcsad.I1wy[idx_par];

        stuffGo->tvcsad.div_xi1[idx_img] = p_data->stuffGo.tvcsad.div_xi1[idx_par];
        stuffGo->tvcsad.div_xi2[idx_img] = p_data->stuffGo.tvcsad.div_xi2[idx_par];

        stuffBa->tvcsad.div_xi1[idx_img] = p_data->stuffBa.tvcsad.div_xi1[idx_par];
        stuffBa->tvcsad.div_xi2[idx_img] = p_data->stuffBa.tvcsad.div_xi2[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the TV-CSAD functional with weights
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_tvcsadw_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                            const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // PosNei (neighbours' position)
        p_data->stuffGo.tvcsadw.pnei[idx_par] = stuffGo->tvcsadw.pnei[idx_img];
        p_data->stuffBa.tvcsadw.pnei[idx_par] = stuffBa->tvcsadw.pnei[idx_img];

        // Xi
        p_data->stuffGo.tvcsadw.xi11[idx_par] = stuffGo->tvcsadw.xi11[idx_img];
        p_data->stuffGo.tvcsadw.xi12[idx_par] = stuffGo->tvcsadw.xi12[idx_img];
        p_data->stuffGo.tvcsadw.xi21[idx_par] = stuffGo->tvcsadw.xi21[idx_img];
        p_data->stuffGo.tvcsadw.xi22[idx_par] = stuffGo->tvcsadw.xi22[idx_img];

        p_data->stuffBa.tvcsadw.xi11[idx_par] = stuffBa->tvcsadw.xi11[idx_img];
        p_data->stuffBa.tvcsadw.xi12[idx_par] = stuffBa->tvcsadw.xi12[idx_img];
        p_data->stuffBa.tvcsadw.xi21[idx_par] = stuffBa->tvcsadw.xi21[idx_img];
        p_data->stuffBa.tvcsadw.xi22[idx_par] = stuffBa->tvcsadw.xi22[idx_img];

        // u1, u2
        p_data->stuffGo.tvcsadw.u1x[idx_par] = stuffGo->tvcsadw.u1x[idx_img];
        p_data->stuffGo.tvcsadw.u1y[idx_par] = stuffGo->tvcsadw.u1y[idx_img];
        p_data->stuffGo.tvcsadw.u2x[idx_par] = stuffGo->tvcsadw.u2x[idx_img];
        p_data->stuffGo.tvcsadw.u2y[idx_par] = stuffGo->tvcsadw.u2y[idx_img];

        p_data->stuffBa.tvcsadw.u1x[idx_par] = stuffBa->tvcsadw.u1x[idx_img];
        p_data->stuffBa.tvcsadw.u1y[idx_par] = stuffBa->tvcsadw.u1y[idx_img];
        p_data->stuffBa.tvcsadw.u2x[idx_par] = stuffBa->tvcsadw.u2x[idx_img];
        p_data->stuffBa.tvcsadw.u2y[idx_par] = stuffBa->tvcsadw.u2y[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.tvcsadw.v1[idx_par] = stuffGo->tvcsadw.v1[idx_img];
        p_data->stuffGo.tvcsadw.v2[idx_par] = stuffGo->tvcsadw.v2[idx_img];

        p_data->stuffBa.tvcsadw.v1[idx_par] = stuffBa->tvcsadw.v1[idx_img];
        p_data->stuffBa.tvcsadw.v2[idx_par] = stuffBa->tvcsadw.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.tvcsadw.rho_c[idx_par] = stuffGo->tvcsadw.rho_c[idx_img];
        p_data->stuffBa.tvcsadw.rho_c[idx_par] = stuffBa->tvcsadw.rho_c[idx_img];

        p_data->stuffGo.tvcsadw.grad[idx_par] = stuffGo->tvcsadw.grad[idx_img];
        p_data->stuffBa.tvcsadw.grad[idx_par] = stuffBa->tvcsadw.grad[idx_img];

        p_data->stuffGo.tvcsadw.u1_[idx_par] = stuffGo->tvcsadw.u1_[idx_img];
        p_data->stuffGo.tvcsadw.u2_[idx_par] = stuffGo->tvcsadw.u2_[idx_img];

        p_data->stuffBa.tvcsadw.u1_[idx_par] = stuffBa->tvcsadw.u1_[idx_img];
        p_data->stuffBa.tvcsadw.u2_[idx_par] = stuffBa->tvcsadw.u2_[idx_img];

        p_data->stuffGo.tvcsadw.u1_tmp[idx_par] = stuffGo->tvcsadw.u1_tmp[idx_img];
        p_data->stuffGo.tvcsadw.u2_tmp[idx_par] = stuffGo->tvcsadw.u2_tmp[idx_img];

        p_data->stuffBa.tvcsadw.u1_tmp[idx_par] = stuffBa->tvcsadw.u1_tmp[idx_img];
        p_data->stuffBa.tvcsadw.u2_tmp[idx_par] = stuffBa->tvcsadw.u2_tmp[idx_img];

        p_data->stuffGo.tvcsadw.I1x[idx_par] = stuffGo->tvcsadw.I1x[idx_img];
        p_data->stuffGo.tvcsadw.I1y[idx_par] = stuffGo->tvcsadw.I1y[idx_img];

        p_data->stuffBa.tvcsadw.I1x[idx_par] = stuffBa->tvcsadw.I1x[idx_img];
        p_data->stuffBa.tvcsadw.I1y[idx_par] = stuffBa->tvcsadw.I1y[idx_img];

        p_data->stuffGo.tvcsadw.I1wx[idx_par] = stuffGo->tvcsadw.I1wx[idx_img];
        p_data->stuffGo.tvcsadw.I1wy[idx_par] = stuffGo->tvcsadw.I1wy[idx_img];

        p_data->stuffBa.tvcsadw.I1wx[idx_par] = stuffBa->tvcsadw.I1wx[idx_img];
        p_data->stuffBa.tvcsadw.I1wy[idx_par] = stuffBa->tvcsadw.I1wy[idx_img];

        p_data->stuffGo.tvcsadw.div_xi1[idx_par] = stuffGo->tvcsadw.div_xi1[idx_img];
        p_data->stuffGo.tvcsadw.div_xi2[idx_par] = stuffGo->tvcsadw.div_xi2[idx_img];

        p_data->stuffBa.tvcsadw.div_xi1[idx_par] = stuffBa->tvcsadw.div_xi1[idx_img];
        p_data->stuffBa.tvcsadw.div_xi2[idx_par] = stuffBa->tvcsadw.div_xi2[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // PosNei
        stuffGo->tvcsadw.pnei[idx_img] = p_data->stuffGo.tvcsadw.pnei[idx_par];
        stuffBa->tvcsadw.pnei[idx_img] = p_data->stuffBa.tvcsadw.pnei[idx_par];

        // Xi
        stuffGo->tvcsadw.xi11[idx_img] = p_data->stuffGo.tvcsadw.xi11[idx_par];
        stuffGo->tvcsadw.xi12[idx_img] = p_data->stuffGo.tvcsadw.xi12[idx_par];
        stuffGo->tvcsadw.xi21[idx_img] = p_data->stuffGo.tvcsadw.xi21[idx_par];
        stuffGo->tvcsadw.xi22[idx_img] = p_data->stuffGo.tvcsadw.xi22[idx_par];

        stuffBa->tvcsadw.xi11[idx_img] = p_data->stuffBa.tvcsadw.xi11[idx_par];
        stuffBa->tvcsadw.xi12[idx_img] = p_data->stuffBa.tvcsadw.xi12[idx_par];
        stuffBa->tvcsadw.xi21[idx_img] = p_data->stuffBa.tvcsadw.xi21[idx_par];
        stuffBa->tvcsadw.xi22[idx_img] = p_data->stuffBa.tvcsadw.xi22[idx_par];

        // u1, u2
        stuffGo->tvcsadw.u1x[idx_img] = p_data->stuffGo.tvcsadw.u1x[idx_par];
        stuffGo->tvcsadw.u1y[idx_img] = p_data->stuffGo.tvcsadw.u1y[idx_par];
        stuffGo->tvcsadw.u2x[idx_img] = p_data->stuffGo.tvcsadw.u2x[idx_par];
        stuffGo->tvcsadw.u2y[idx_img] = p_data->stuffGo.tvcsadw.u2y[idx_par];

        stuffBa->tvcsadw.u1x[idx_img] = p_data->stuffBa.tvcsadw.u1x[idx_par];
        stuffBa->tvcsadw.u1y[idx_img] = p_data->stuffBa.tvcsadw.u1y[idx_par];
        stuffBa->tvcsadw.u2x[idx_img] = p_data->stuffBa.tvcsadw.u2x[idx_par];
        stuffBa->tvcsadw.u2y[idx_img] = p_data->stuffBa.tvcsadw.u2y[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->tvcsadw.v1[idx_img] = p_data->stuffGo.tvcsadw.v1[idx_par];
        stuffGo->tvcsadw.v2[idx_img] = p_data->stuffGo.tvcsadw.v2[idx_par];

        stuffBa->tvcsadw.v1[idx_img] = p_data->stuffBa.tvcsadw.v1[idx_par];
        stuffBa->tvcsadw.v2[idx_img] = p_data->stuffBa.tvcsadw.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->tvcsadw.rho_c[idx_img] = p_data->stuffGo.tvcsadw.rho_c[idx_par];
        stuffBa->tvcsadw.rho_c[idx_img] = p_data->stuffBa.tvcsadw.rho_c[idx_par];

        stuffGo->tvcsadw.grad[idx_img] = p_data->stuffGo.tvcsadw.grad[idx_par];
        stuffBa->tvcsadw.grad[idx_img] = p_data->stuffBa.tvcsadw.grad[idx_par];

        stuffGo->tvcsadw.u1_[idx_img] = p_data->stuffGo.tvcsadw.u1_[idx_par];
        stuffGo->tvcsadw.u2_[idx_img] = p_data->stuffGo.tvcsadw.u2_[idx_par];

        stuffBa->tvcsadw.u1_[idx_img] = p_data->stuffBa.tvcsadw.u1_[idx_par];
        stuffBa->tvcsadw.u2_[idx_img] = p_data->stuffBa.tvcsadw.u2_[idx_par];

        stuffGo->tvcsadw.u1_tmp[idx_img] = p_data->stuffGo.tvcsadw.u1_tmp[idx_par];
        stuffGo->tvcsadw.u2_tmp[idx_img] = p_data->stuffGo.tvcsadw.u2_tmp[idx_par];

        stuffBa->tvcsadw.u1_tmp[idx_img] = p_data->stuffBa.tvcsadw.u1_tmp[idx_par];
        stuffBa->tvcsadw.u2_tmp[idx_img] = p_data->stuffBa.tvcsadw.u2_tmp[idx_par];

        stuffGo->tvcsadw.I1x[idx_img] = p_data->stuffGo.tvcsadw.I1x[idx_par];
        stuffGo->tvcsadw.I1y[idx_img] = p_data->stuffGo.tvcsadw.I1y[idx_par];

        stuffBa->tvcsadw.I1x[idx_img] = p_data->stuffBa.tvcsadw.I1x[idx_par];
        stuffBa->tvcsadw.I1y[idx_img] = p_data->stuffBa.tvcsadw.I1y[idx_par];

        stuffGo->tvcsadw.I1wx[idx_img] = p_data->stuffGo.tvcsadw.I1wx[idx_par];
        stuffGo->tvcsadw.I1wy[idx_img] = p_data->stuffGo.tvcsadw.I1wy[idx_par];

        stuffBa->tvcsadw.I1wx[idx_img] = p_data->stuffBa.tvcsadw.I1wx[idx_par];
        stuffBa->tvcsadw.I1wy[idx_img] = p_data->stuffBa.tvcsadw.I1wy[idx_par];

        stuffGo->tvcsadw.div_xi1[idx_img] = p_data->stuffGo.tvcsadw.div_xi1[idx_par];
        stuffGo->tvcsadw.div_xi2[idx_img] = p_data->stuffGo.tvcsadw.div_xi2[idx_par];

        stuffBa->tvcsadw.div_xi1[idx_img] = p_data->stuffBa.tvcsadw.div_xi1[idx_par];
        stuffBa->tvcsadw.div_xi2[idx_img] = p_data->stuffBa.tvcsadw.div_xi2[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the NLTV-CSAD functional
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_nltvcsad_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                             const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Dual variables
        p_data->stuffGo.nltvcsad.p[idx_par] = stuffGo->nltvcsad.p[idx_img];
        p_data->stuffGo.nltvcsad.q[idx_par] = stuffGo->nltvcsad.q[idx_img];

        p_data->stuffBa.nltvcsad.p[idx_par] = stuffBa->nltvcsad.p[idx_img];
        p_data->stuffBa.nltvcsad.q[idx_par] = stuffBa->nltvcsad.q[idx_img];

        // PosNei (neighbours' position)
        p_data->stuffGo.nltvcsad.pnei[idx_par] = stuffGo->nltvcsad.pnei[idx_img];
        p_data->stuffBa.nltvcsad.pnei[idx_par] = stuffBa->nltvcsad.pnei[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.nltvcsad.v1[idx_par] = stuffGo->nltvcsad.v1[idx_img];
        p_data->stuffGo.nltvcsad.v2[idx_par] = stuffGo->nltvcsad.v2[idx_img];

        p_data->stuffBa.nltvcsad.v1[idx_par] = stuffBa->nltvcsad.v1[idx_img];
        p_data->stuffBa.nltvcsad.v2[idx_par] = stuffBa->nltvcsad.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.nltvcsad.rho_c[idx_par] = stuffGo->nltvcsad.rho_c[idx_img];
        p_data->stuffBa.nltvcsad.rho_c[idx_par] = stuffBa->nltvcsad.rho_c[idx_img];

        p_data->stuffGo.nltvcsad.grad[idx_par] = stuffGo->nltvcsad.grad[idx_img];
        p_data->stuffBa.nltvcsad.grad[idx_par] = stuffBa->nltvcsad.grad[idx_img];

        p_data->stuffGo.nltvcsad.u1_[idx_par] = stuffGo->nltvcsad.u1_[idx_img];
        p_data->stuffGo.nltvcsad.u2_[idx_par] = stuffGo->nltvcsad.u2_[idx_img];

        p_data->stuffBa.nltvcsad.u1_[idx_par] = stuffBa->nltvcsad.u1_[idx_img];
        p_data->stuffBa.nltvcsad.u2_[idx_par] = stuffBa->nltvcsad.u2_[idx_img];

        p_data->stuffGo.nltvcsad.u1_tmp[idx_par] = stuffGo->nltvcsad.u1_tmp[idx_img];
        p_data->stuffGo.nltvcsad.u2_tmp[idx_par] = stuffGo->nltvcsad.u2_tmp[idx_img];

        p_data->stuffBa.nltvcsad.u1_tmp[idx_par] = stuffBa->nltvcsad.u1_tmp[idx_img];
        p_data->stuffBa.nltvcsad.u2_tmp[idx_par] = stuffBa->nltvcsad.u2_tmp[idx_img];

        p_data->stuffGo.nltvcsad.I1x[idx_par] = stuffGo->nltvcsad.I1x[idx_img];
        p_data->stuffGo.nltvcsad.I1y[idx_par] = stuffGo->nltvcsad.I1y[idx_img];

        p_data->stuffBa.nltvcsad.I1x[idx_par] = stuffBa->nltvcsad.I1x[idx_img];
        p_data->stuffBa.nltvcsad.I1y[idx_par] = stuffBa->nltvcsad.I1y[idx_img];

        p_data->stuffGo.nltvcsad.I1wx[idx_par] = stuffGo->nltvcsad.I1wx[idx_img];
        p_data->stuffGo.nltvcsad.I1wy[idx_par] = stuffGo->nltvcsad.I1wy[idx_img];

        p_data->stuffBa.nltvcsad.I1wx[idx_par] = stuffBa->nltvcsad.I1wx[idx_img];
        p_data->stuffBa.nltvcsad.I1wy[idx_par] = stuffBa->nltvcsad.I1wy[idx_img];

        p_data->stuffGo.nltvcsad.div_p[idx_par] = stuffGo->nltvcsad.div_p[idx_img];
        p_data->stuffGo.nltvcsad.div_q[idx_par] = stuffGo->nltvcsad.div_q[idx_img];

        p_data->stuffBa.nltvcsad.div_p[idx_par] = stuffBa->nltvcsad.div_p[idx_img];
        p_data->stuffBa.nltvcsad.div_q[idx_par] = stuffBa->nltvcsad.div_q[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // Dual variables
        stuffGo->nltvcsad.p[idx_img] = p_data->stuffGo.nltvcsad.p[idx_par];
        stuffGo->nltvcsad.q[idx_img] = p_data->stuffGo.nltvcsad.q[idx_par];

        stuffBa->nltvcsad.p[idx_img] = p_data->stuffBa.nltvcsad.p[idx_par];
        stuffBa->nltvcsad.q[idx_img] = p_data->stuffBa.nltvcsad.q[idx_par];

        // PosNei
        stuffGo->tvcsadw.pnei[idx_img] = p_data->stuffGo.tvcsadw.pnei[idx_par];
        stuffBa->tvcsadw.pnei[idx_img] = p_data->stuffBa.tvcsadw.pnei[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->nltvcsad.v1[idx_img] = p_data->stuffGo.nltvcsad.v1[idx_par];
        stuffGo->nltvcsad.v2[idx_img] = p_data->stuffGo.nltvcsad.v2[idx_par];

        stuffBa->nltvcsad.v1[idx_img] = p_data->stuffBa.nltvcsad.v1[idx_par];
        stuffBa->nltvcsad.v2[idx_img] = p_data->stuffBa.nltvcsad.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->nltvcsad.rho_c[idx_img] = p_data->stuffGo.nltvcsad.rho_c[idx_par];
        stuffBa->nltvcsad.rho_c[idx_img] = p_data->stuffBa.nltvcsad.rho_c[idx_par];

        stuffGo->nltvcsad.grad[idx_img] = p_data->stuffGo.nltvcsad.grad[idx_par];
        stuffBa->nltvcsad.grad[idx_img] = p_data->stuffBa.nltvcsad.grad[idx_par];

        stuffGo->nltvcsad.u1_[idx_img] = p_data->stuffGo.nltvcsad.u1_[idx_par];
        stuffGo->nltvcsad.u2_[idx_img] = p_data->stuffGo.nltvcsad.u2_[idx_par];

        stuffBa->nltvcsad.u1_[idx_img] = p_data->stuffBa.nltvcsad.u1_[idx_par];
        stuffBa->nltvcsad.u2_[idx_img] = p_data->stuffBa.nltvcsad.u2_[idx_par];

        stuffGo->nltvcsad.u1_tmp[idx_img] = p_data->stuffGo.nltvcsad.u1_tmp[idx_par];
        stuffGo->nltvcsad.u2_tmp[idx_img] = p_data->stuffGo.nltvcsad.u2_tmp[idx_par];

        stuffBa->nltvcsad.u1_tmp[idx_img] = p_data->stuffBa.nltvcsad.u1_tmp[idx_par];
        stuffBa->nltvcsad.u2_tmp[idx_img] = p_data->stuffBa.nltvcsad.u2_tmp[idx_par];

        stuffGo->nltvcsad.I1x[idx_img] = p_data->stuffGo.nltvcsad.I1x[idx_par];
        stuffGo->nltvcsad.I1y[idx_img] = p_data->stuffGo.nltvcsad.I1y[idx_par];

        stuffBa->nltvcsad.I1x[idx_img] = p_data->stuffBa.nltvcsad.I1x[idx_par];
        stuffBa->nltvcsad.I1y[idx_img] = p_data->stuffBa.nltvcsad.I1y[idx_par];

        stuffGo->nltvcsad.I1wx[idx_img] = p_data->stuffGo.nltvcsad.I1wx[idx_par];
        stuffGo->nltvcsad.I1wy[idx_img] = p_data->stuffGo.nltvcsad.I1wy[idx_par];

        stuffBa->nltvcsad.I1wx[idx_img] = p_data->stuffBa.nltvcsad.I1wx[idx_par];
        stuffBa->nltvcsad.I1wy[idx_img] = p_data->stuffBa.nltvcsad.I1wy[idx_par];

        stuffGo->nltvcsad.div_p[idx_img] = p_data->stuffGo.nltvcsad.div_p[idx_par];
        stuffGo->nltvcsad.div_q[idx_img] = p_data->stuffGo.nltvcsad.div_q[idx_par];

        stuffBa->nltvcsad.div_p[idx_img] = p_data->stuffBa.nltvcsad.div_p[idx_par];
        stuffBa->nltvcsad.div_q[idx_img] = p_data->stuffBa.nltvcsad.div_q[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the NLTV-CSAD functional with weights
 *
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_nltvcsadw_stuffof(SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa, PartitionData *& p_data,
                              const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Copy image-wise variables to partition-specific ones
        // Dual variables
        p_data->stuffGo.nltvcsadw.p[idx_par] = stuffGo->nltvcsadw.p[idx_img];
        p_data->stuffGo.nltvcsadw.q[idx_par] = stuffGo->nltvcsadw.q[idx_img];

        p_data->stuffBa.nltvcsadw.p[idx_par] = stuffBa->nltvcsadw.p[idx_img];
        p_data->stuffBa.nltvcsadw.q[idx_par] = stuffBa->nltvcsadw.q[idx_img];

        // PosNei (neighbours' position)
        p_data->stuffGo.nltvcsadw.pnei[idx_par] = stuffGo->nltvcsadw.pnei[idx_img];
        p_data->stuffBa.nltvcsadw.pnei[idx_par] = stuffBa->nltvcsadw.pnei[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.nltvcsadw.v1[idx_par] = stuffGo->nltvcsadw.v1[idx_img];
        p_data->stuffGo.nltvcsadw.v2[idx_par] = stuffGo->nltvcsadw.v2[idx_img];

        p_data->stuffBa.nltvcsadw.v1[idx_par] = stuffBa->nltvcsadw.v1[idx_img];
        p_data->stuffBa.nltvcsadw.v2[idx_par] = stuffBa->nltvcsadw.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        p_data->stuffGo.nltvcsadw.rho_c[idx_par] = stuffGo->nltvcsadw.rho_c[idx_img];
        p_data->stuffBa.nltvcsadw.rho_c[idx_par] = stuffBa->nltvcsadw.rho_c[idx_img];

        p_data->stuffGo.nltvcsadw.grad[idx_par] = stuffGo->nltvcsadw.grad[idx_img];
        p_data->stuffBa.nltvcsadw.grad[idx_par] = stuffBa->nltvcsadw.grad[idx_img];

        p_data->stuffGo.nltvcsadw.u1_[idx_par] = stuffGo->nltvcsadw.u1_[idx_img];
        p_data->stuffGo.nltvcsadw.u2_[idx_par] = stuffGo->nltvcsadw.u2_[idx_img];

        p_data->stuffBa.nltvcsadw.u1_[idx_par] = stuffBa->nltvcsadw.u1_[idx_img];
        p_data->stuffBa.nltvcsadw.u2_[idx_par] = stuffBa->nltvcsadw.u2_[idx_img];

        p_data->stuffGo.nltvcsadw.u1_tmp[idx_par] = stuffGo->nltvcsadw.u1_tmp[idx_img];
        p_data->stuffGo.nltvcsadw.u2_tmp[idx_par] = stuffGo->nltvcsadw.u2_tmp[idx_img];

        p_data->stuffBa.nltvcsadw.u1_tmp[idx_par] = stuffBa->nltvcsadw.u1_tmp[idx_img];
        p_data->stuffBa.nltvcsadw.u2_tmp[idx_par] = stuffBa->nltvcsadw.u2_tmp[idx_img];

        p_data->stuffGo.nltvcsadw.I1x[idx_par] = stuffGo->nltvcsadw.I1x[idx_img];
        p_data->stuffGo.nltvcsadw.I1y[idx_par] = stuffGo->nltvcsadw.I1y[idx_img];

        p_data->stuffBa.nltvcsadw.I1x[idx_par] = stuffBa->nltvcsadw.I1x[idx_img];
        p_data->stuffBa.nltvcsadw.I1y[idx_par] = stuffBa->nltvcsadw.I1y[idx_img];

        p_data->stuffGo.nltvcsadw.I1wx[idx_par] = stuffGo->nltvcsadw.I1wx[idx_img];
        p_data->stuffGo.nltvcsadw.I1wy[idx_par] = stuffGo->nltvcsadw.I1wy[idx_img];

        p_data->stuffBa.nltvcsadw.I1wx[idx_par] = stuffBa->nltvcsadw.I1wx[idx_img];
        p_data->stuffBa.nltvcsadw.I1wy[idx_par] = stuffBa->nltvcsadw.I1wy[idx_img];

        p_data->stuffGo.nltvcsadw.div_p[idx_par] = stuffGo->nltvcsadw.div_p[idx_img];
        p_data->stuffGo.nltvcsadw.div_q[idx_par] = stuffGo->nltvcsadw.div_q[idx_img];

        p_data->stuffBa.nltvcsadw.div_p[idx_par] = stuffBa->nltvcsadw.div_p[idx_img];
        p_data->stuffBa.nltvcsadw.div_q[idx_par] = stuffBa->nltvcsadw.div_q[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // Dual variables
        stuffGo->nltvcsadw.p[idx_img] = p_data->stuffGo.nltvcsadw.p[idx_par];
        stuffGo->nltvcsadw.q[idx_img] = p_data->stuffGo.nltvcsadw.q[idx_par];

        stuffBa->nltvcsadw.p[idx_img] = p_data->stuffBa.nltvcsadw.p[idx_par];
        stuffBa->nltvcsadw.q[idx_img] = p_data->stuffBa.nltvcsadw.q[idx_par];

        // PosNei
        stuffGo->tvcsadw.pnei[idx_img] = p_data->stuffGo.tvcsadw.pnei[idx_par];
        stuffBa->tvcsadw.pnei[idx_img] = p_data->stuffBa.tvcsadw.pnei[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->nltvcsadw.v1[idx_img] = p_data->stuffGo.nltvcsadw.v1[idx_par];
        stuffGo->nltvcsadw.v2[idx_img] = p_data->stuffGo.nltvcsadw.v2[idx_par];

        stuffBa->nltvcsadw.v1[idx_img] = p_data->stuffBa.nltvcsadw.v1[idx_par];
        stuffBa->nltvcsadw.v2[idx_img] = p_data->stuffBa.nltvcsadw.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        stuffGo->nltvcsadw.rho_c[idx_img] = p_data->stuffGo.nltvcsadw.rho_c[idx_par];
        stuffBa->nltvcsadw.rho_c[idx_img] = p_data->stuffBa.nltvcsadw.rho_c[idx_par];

        stuffGo->nltvcsadw.grad[idx_img] = p_data->stuffGo.nltvcsadw.grad[idx_par];
        stuffBa->nltvcsadw.grad[idx_img] = p_data->stuffBa.nltvcsadw.grad[idx_par];

        stuffGo->nltvcsadw.u1_[idx_img] = p_data->stuffGo.nltvcsadw.u1_[idx_par];
        stuffGo->nltvcsadw.u2_[idx_img] = p_data->stuffGo.nltvcsadw.u2_[idx_par];

        stuffBa->nltvcsadw.u1_[idx_img] = p_data->stuffBa.nltvcsadw.u1_[idx_par];
        stuffBa->nltvcsadw.u2_[idx_img] = p_data->stuffBa.nltvcsadw.u2_[idx_par];

        stuffGo->nltvcsadw.u1_tmp[idx_img] = p_data->stuffGo.nltvcsadw.u1_tmp[idx_par];
        stuffGo->nltvcsadw.u2_tmp[idx_img] = p_data->stuffGo.nltvcsadw.u2_tmp[idx_par];

        stuffBa->nltvcsadw.u1_tmp[idx_img] = p_data->stuffBa.nltvcsadw.u1_tmp[idx_par];
        stuffBa->nltvcsadw.u2_tmp[idx_img] = p_data->stuffBa.nltvcsadw.u2_tmp[idx_par];

        stuffGo->nltvcsadw.I1x[idx_img] = p_data->stuffGo.nltvcsadw.I1x[idx_par];
        stuffGo->nltvcsadw.I1y[idx_img] = p_data->stuffGo.nltvcsadw.I1y[idx_par];

        stuffBa->nltvcsadw.I1x[idx_img] = p_data->stuffBa.nltvcsadw.I1x[idx_par];
        stuffBa->nltvcsadw.I1y[idx_img] = p_data->stuffBa.nltvcsadw.I1y[idx_par];

        stuffGo->nltvcsadw.I1wx[idx_img] = p_data->stuffGo.nltvcsadw.I1wx[idx_par];
        stuffGo->nltvcsadw.I1wy[idx_img] = p_data->stuffGo.nltvcsadw.I1wy[idx_par];

        stuffBa->nltvcsadw.I1wx[idx_img] = p_data->stuffBa.nltvcsadw.I1wx[idx_par];
        stuffBa->nltvcsadw.I1wy[idx_img] = p_data->stuffBa.nltvcsadw.I1wy[idx_par];

        stuffGo->nltvcsadw.div_p[idx_img] = p_data->stuffGo.nltvcsadw.div_p[idx_par];
        stuffGo->nltvcsadw.div_q[idx_img] = p_data->stuffGo.nltvcsadw.div_q[idx_par];

        stuffBa->nltvcsadw.div_p[idx_img] = p_data->stuffBa.nltvcsadw.div_p[idx_par];
        stuffBa->nltvcsadw.div_q[idx_img] = p_data->stuffBa.nltvcsadw.div_q[idx_par];
    }
}


/**
 * @brief                   updates SpecificOFStuff struct for the TVL1 functional with occlusions
 *
 * @param ofGo              forward flow' OpticalFlowData struct (used only to check algorithm's step: global or local)
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_of_data
 */
void update_tvl2occ_stuffof(OpticalFlowData *& ofGo, SpecificOFStuff *& stuffGo, SpecificOFStuff *& stuffBa,
                            PartitionData *& p_data, const int idx_img, const int idx_par, bool img_to_part)
{
    if (img_to_part) {
        // Chi
        p_data->stuffGo.tvl2_occ.chix[idx_par] = stuffGo->tvl2_occ.chix[idx_img];
        p_data->stuffGo.tvl2_occ.chiy[idx_par] = stuffGo->tvl2_occ.chiy[idx_img];

        p_data->stuffBa.tvl2_occ.chix[idx_par] = stuffBa->tvl2_occ.chix[idx_img];
        p_data->stuffBa.tvl2_occ.chiy[idx_par] = stuffBa->tvl2_occ.chiy[idx_img];

        // Diff u_N
        p_data->stuffGo.tvl2_occ.diff_u_N[idx_par] = stuffGo->tvl2_occ.diff_u_N[idx_img];
        p_data->stuffBa.tvl2_occ.diff_u_N[idx_par] = stuffBa->tvl2_occ.diff_u_N[idx_img];

        // g
        p_data->stuffGo.tvl2_occ.g[idx_par] = stuffGo->tvl2_occ.g[idx_img];
        p_data->stuffBa.tvl2_occ.g[idx_par] = stuffBa->tvl2_occ.g[idx_img];

        // Xi
        p_data->stuffGo.tvl2_occ.xi11[idx_par] = stuffGo->tvl2_occ.xi11[idx_img];
        p_data->stuffGo.tvl2_occ.xi12[idx_par] = stuffGo->tvl2_occ.xi12[idx_img];
        p_data->stuffGo.tvl2_occ.xi21[idx_par] = stuffGo->tvl2_occ.xi21[idx_img];
        p_data->stuffGo.tvl2_occ.xi22[idx_par] = stuffGo->tvl2_occ.xi22[idx_img];

        p_data->stuffBa.tvl2_occ.xi11[idx_par] = stuffBa->tvl2_occ.xi11[idx_img];
        p_data->stuffBa.tvl2_occ.xi12[idx_par] = stuffBa->tvl2_occ.xi12[idx_img];
        p_data->stuffBa.tvl2_occ.xi21[idx_par] = stuffBa->tvl2_occ.xi21[idx_img];
        p_data->stuffBa.tvl2_occ.xi22[idx_par] = stuffBa->tvl2_occ.xi22[idx_img];

        // u1, u2
        p_data->stuffGo.tvl2_occ.u1x[idx_par] = stuffGo->tvl2_occ.u1x[idx_img];
        p_data->stuffGo.tvl2_occ.u1y[idx_par] = stuffGo->tvl2_occ.u1y[idx_img];
        p_data->stuffGo.tvl2_occ.u2x[idx_par] = stuffGo->tvl2_occ.u2x[idx_img];
        p_data->stuffGo.tvl2_occ.u2y[idx_par] = stuffGo->tvl2_occ.u2y[idx_img];

        p_data->stuffBa.tvl2_occ.u1x[idx_par] = stuffBa->tvl2_occ.u1x[idx_img];
        p_data->stuffBa.tvl2_occ.u1y[idx_par] = stuffBa->tvl2_occ.u1y[idx_img];
        p_data->stuffBa.tvl2_occ.u2x[idx_par] = stuffBa->tvl2_occ.u2x[idx_img];
        p_data->stuffBa.tvl2_occ.u2y[idx_par] = stuffBa->tvl2_occ.u2y[idx_img];

        // v1, v2 (auxiliar minimization variables)
        p_data->stuffGo.tvl2_occ.v1[idx_par] = stuffGo->tvl2_occ.v1[idx_img];
        p_data->stuffGo.tvl2_occ.v2[idx_par] = stuffGo->tvl2_occ.v2[idx_img];

        p_data->stuffBa.tvl2_occ.v1[idx_par] = stuffBa->tvl2_occ.v1[idx_img];
        p_data->stuffBa.tvl2_occ.v2[idx_par] = stuffBa->tvl2_occ.v2[idx_img];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        // rho_c1
        p_data->stuffGo.tvl2_occ.rho_c1[idx_par] = stuffGo->tvl2_occ.rho_c1[idx_img];
        p_data->stuffBa.tvl2_occ.rho_c1[idx_par] = stuffBa->tvl2_occ.rho_c1[idx_img];

        // rho_c_1
        p_data->stuffGo.tvl2_occ.rho_c_1[idx_par] = stuffGo->tvl2_occ.rho_c_1[idx_img];
        p_data->stuffBa.tvl2_occ.rho_c_1[idx_par] = stuffBa->tvl2_occ.rho_c_1[idx_img];

        // grad_1
        p_data->stuffGo.tvl2_occ.grad_1[idx_par] = stuffGo->tvl2_occ.grad_1[idx_img];
        p_data->stuffBa.tvl2_occ.grad_1[idx_par] = stuffBa->tvl2_occ.grad_1[idx_img];

        // grad__1
        p_data->stuffGo.tvl2_occ.grad__1[idx_par] = stuffGo->tvl2_occ.grad__1[idx_img];
        p_data->stuffBa.tvl2_occ.grad__1[idx_par] = stuffBa->tvl2_occ.grad__1[idx_img];

        if (ofGo->params.step_algorithm == GLOBAL_STEP) {
            p_data->stuffGo.tvl2_occ.I0x[idx_par] = stuffGo->tvl2_occ.I0x[idx_img];
            p_data->stuffGo.tvl2_occ.I0y[idx_par] = stuffGo->tvl2_occ.I0y[idx_img];

            p_data->stuffBa.tvl2_occ.I0x[idx_par] = stuffBa->tvl2_occ.I0x[idx_img];
            p_data->stuffBa.tvl2_occ.I0y[idx_par] = stuffBa->tvl2_occ.I0y[idx_img];
        }
        // I1
        stuffGo->tvl2_occ.I1x[idx_img] = p_data->stuffGo.tvl2_occ.I1x[idx_par];
        stuffGo->tvl2_occ.I1y[idx_img] = p_data->stuffGo.tvl2_occ.I1y[idx_par];

        stuffBa->tvl2_occ.I1x[idx_img] = p_data->stuffBa.tvl2_occ.I1x[idx_par];
        stuffBa->tvl2_occ.I1y[idx_img] = p_data->stuffBa.tvl2_occ.I1y[idx_par];

        stuffGo->tvl2_occ.I1wx[idx_img] = p_data->stuffGo.tvl2_occ.I1wx[idx_par];
        stuffGo->tvl2_occ.I1wy[idx_img] = p_data->stuffGo.tvl2_occ.I1wy[idx_par];

        stuffBa->tvl2_occ.I1wx[idx_img] = p_data->stuffBa.tvl2_occ.I1wx[idx_par];
        stuffBa->tvl2_occ.I1wy[idx_img] = p_data->stuffBa.tvl2_occ.I1wy[idx_par];

        // I_1
        stuffGo->tvl2_occ.I_1x[idx_img] = p_data->stuffGo.tvl2_occ.I_1x[idx_par];
        stuffGo->tvl2_occ.I_1y[idx_img] = p_data->stuffGo.tvl2_occ.I_1y[idx_par];

        stuffBa->tvl2_occ.I_1x[idx_img] = p_data->stuffBa.tvl2_occ.I_1x[idx_par];
        stuffBa->tvl2_occ.I_1y[idx_img] = p_data->stuffBa.tvl2_occ.I_1y[idx_par];

        stuffGo->tvl2_occ.I_1wx[idx_img] = p_data->stuffGo.tvl2_occ.I_1wx[idx_par];
        stuffGo->tvl2_occ.I_1wy[idx_img] = p_data->stuffGo.tvl2_occ.I_1wy[idx_par];

        stuffBa->tvl2_occ.I_1wx[idx_img] = p_data->stuffBa.tvl2_occ.I_1wx[idx_par];
        stuffBa->tvl2_occ.I_1wy[idx_img] = p_data->stuffBa.tvl2_occ.I_1wy[idx_par];

        // vi_div1
        p_data->stuffGo.tvl2_occ.vi_div1[idx_par] = stuffGo->tvl2_occ.vi_div1[idx_img];
        p_data->stuffBa.tvl2_occ.vi_div1[idx_par] = stuffBa->tvl2_occ.vi_div1[idx_img];

        // grad_x1
        p_data->stuffGo.tvl2_occ.grad_x1[idx_par] = stuffGo->tvl2_occ.grad_x1[idx_img];
        p_data->stuffBa.tvl2_occ.grad_x1[idx_par] = stuffBa->tvl2_occ.grad_x1[idx_img];

        // grad_y1
        p_data->stuffGo.tvl2_occ.grad_y1[idx_par] = stuffGo->tvl2_occ.grad_y1[idx_img];
        p_data->stuffBa.tvl2_occ.grad_y1[idx_par] = stuffBa->tvl2_occ.grad_y1[idx_img];

        // vi_div2
        p_data->stuffGo.tvl2_occ.vi_div2[idx_par] = stuffGo->tvl2_occ.vi_div2[idx_img];
        p_data->stuffBa.tvl2_occ.vi_div2[idx_par] = stuffBa->tvl2_occ.vi_div2[idx_img];

        // grad_x2
        p_data->stuffGo.tvl2_occ.grad_x2[idx_par] = stuffGo->tvl2_occ.grad_x2[idx_img];
        p_data->stuffBa.tvl2_occ.grad_x2[idx_par] = stuffBa->tvl2_occ.grad_x2[idx_img];

        // grad_y2
        p_data->stuffGo.tvl2_occ.grad_y2[idx_par] = stuffGo->tvl2_occ.grad_y2[idx_img];
        p_data->stuffBa.tvl2_occ.grad_y2[idx_par] = stuffBa->tvl2_occ.grad_y2[idx_img];

        // g_xi
        p_data->stuffGo.tvl2_occ.g_xi11[idx_img] = stuffGo->tvl2_occ.g_xi11[idx_par];
        p_data->stuffGo.tvl2_occ.g_xi12[idx_img] = stuffGo->tvl2_occ.g_xi12[idx_par];

        p_data->stuffBa.tvl2_occ.g_xi11[idx_img] = stuffBa->tvl2_occ.g_xi11[idx_par];
        p_data->stuffBa.tvl2_occ.g_xi12[idx_img] = stuffBa->tvl2_occ.g_xi12[idx_par];

        p_data->stuffGo.tvl2_occ.g_xi21[idx_img] = stuffGo->tvl2_occ.g_xi21[idx_par];
        p_data->stuffGo.tvl2_occ.g_xi22[idx_img] = stuffGo->tvl2_occ.g_xi22[idx_par];

        p_data->stuffBa.tvl2_occ.g_xi21[idx_img] = stuffBa->tvl2_occ.g_xi21[idx_par];
        p_data->stuffBa.tvl2_occ.g_xi22[idx_img] = stuffBa->tvl2_occ.g_xi22[idx_par];

        // div_g_xi
        p_data->stuffGo.tvl2_occ.div_g_xi1[idx_img] = stuffGo->tvl2_occ.div_g_xi1[idx_par];
        p_data->stuffGo.tvl2_occ.div_g_xi2[idx_img] = stuffGo->tvl2_occ.div_g_xi2[idx_par];

        p_data->stuffBa.tvl2_occ.div_g_xi1[idx_img] = stuffBa->tvl2_occ.div_g_xi1[idx_par];
        p_data->stuffBa.tvl2_occ.div_g_xi2[idx_img] = stuffBa->tvl2_occ.div_g_xi2[idx_par];

        // eta
        p_data->stuffGo.tvl2_occ.eta1[idx_img] = stuffGo->tvl2_occ.eta1[idx_par];
        p_data->stuffGo.tvl2_occ.eta2[idx_img] = stuffGo->tvl2_occ.eta2[idx_par];

        p_data->stuffBa.tvl2_occ.eta1[idx_img] = stuffBa->tvl2_occ.eta1[idx_par];
        p_data->stuffBa.tvl2_occ.eta2[idx_img] = stuffBa->tvl2_occ.eta2[idx_par];

        // F, G
        p_data->stuffGo.tvl2_occ.F[idx_img] = stuffGo->tvl2_occ.F[idx_par];
        p_data->stuffGo.tvl2_occ.G[idx_img] = stuffGo->tvl2_occ.G[idx_par];

        p_data->stuffBa.tvl2_occ.F[idx_img] = stuffBa->tvl2_occ.F[idx_par];
        p_data->stuffBa.tvl2_occ.G[idx_img] = stuffBa->tvl2_occ.G[idx_par];

        // div_u
        p_data->stuffGo.tvl2_occ.div_u[idx_par] = stuffGo->tvl2_occ.div_u[idx_img];
        p_data->stuffBa.tvl2_occ.div_u[idx_par] = stuffBa->tvl2_occ.div_u[idx_img];

        // g_eta
        p_data->stuffGo.tvl2_occ.g_eta1[idx_img] = stuffGo->tvl2_occ.g_eta1[idx_par];
        p_data->stuffGo.tvl2_occ.g_eta2[idx_img] = stuffGo->tvl2_occ.g_eta2[idx_par];

        p_data->stuffBa.tvl2_occ.g_eta1[idx_img] = stuffBa->tvl2_occ.g_eta1[idx_par];
        p_data->stuffBa.tvl2_occ.g_eta2[idx_img] = stuffBa->tvl2_occ.g_eta2[idx_par];

        // div_g_eta
        p_data->stuffGo.tvl2_occ.div_g_eta[idx_par] = stuffGo->tvl2_occ.div_g_eta[idx_img];
        p_data->stuffBa.tvl2_occ.div_g_eta[idx_par] = stuffBa->tvl2_occ.div_g_eta[idx_img];

    } else {
        // Copy partition-wise variables to corresponding image-wise variables
        // Chi
        stuffGo->tvl2_occ.chix[idx_img] = p_data->stuffGo.tvl2_occ.chix[idx_par];
        stuffGo->tvl2_occ.chiy[idx_img] = p_data->stuffGo.tvl2_occ.chiy[idx_par];

        stuffBa->tvl2_occ.chix[idx_img] = p_data->stuffBa.tvl2_occ.chix[idx_par];
        stuffBa->tvl2_occ.chiy[idx_img] = p_data->stuffBa.tvl2_occ.chiy[idx_par];

        // Diff u_N
        stuffGo->tvl2_occ.diff_u_N[idx_img] = p_data->stuffGo.tvl2_occ.diff_u_N[idx_par];
        stuffBa->tvl2_occ.diff_u_N[idx_img] = p_data->stuffBa.tvl2_occ.diff_u_N[idx_par];

        // g
        stuffGo->tvl2_occ.g[idx_img] = p_data->stuffGo.tvl2_occ.g[idx_par];
        stuffBa->tvl2_occ.g[idx_img] = p_data->stuffBa.tvl2_occ.g[idx_par];

        // Xi
        stuffGo->tvl2_occ.xi11[idx_img] = p_data->stuffGo.tvl2_occ.xi11[idx_par];
        stuffGo->tvl2_occ.xi12[idx_img] = p_data->stuffGo.tvl2_occ.xi12[idx_par];
        stuffGo->tvl2_occ.xi21[idx_img] = p_data->stuffGo.tvl2_occ.xi21[idx_par];
        stuffGo->tvl2_occ.xi22[idx_img] = p_data->stuffGo.tvl2_occ.xi22[idx_par];

        stuffBa->tvl2_occ.xi11[idx_img] = p_data->stuffBa.tvl2_occ.xi11[idx_par];
        stuffBa->tvl2_occ.xi12[idx_img] = p_data->stuffBa.tvl2_occ.xi12[idx_par];
        stuffBa->tvl2_occ.xi21[idx_img] = p_data->stuffBa.tvl2_occ.xi21[idx_par];
        stuffBa->tvl2_occ.xi22[idx_img] = p_data->stuffBa.tvl2_occ.xi22[idx_par];

        // u1, u2
        stuffGo->tvl2_occ.u1x[idx_img] = p_data->stuffGo.tvl2_occ.u1x[idx_par];
        stuffGo->tvl2_occ.u1y[idx_img] = p_data->stuffGo.tvl2_occ.u1y[idx_par];
        stuffGo->tvl2_occ.u2x[idx_img] = p_data->stuffGo.tvl2_occ.u2x[idx_par];
        stuffGo->tvl2_occ.u2y[idx_img] = p_data->stuffGo.tvl2_occ.u2y[idx_par];

        stuffBa->tvl2_occ.u1x[idx_img] = p_data->stuffBa.tvl2_occ.u1x[idx_par];
        stuffBa->tvl2_occ.u1y[idx_img] = p_data->stuffBa.tvl2_occ.u1y[idx_par];
        stuffBa->tvl2_occ.u2x[idx_img] = p_data->stuffBa.tvl2_occ.u2x[idx_par];
        stuffBa->tvl2_occ.u2y[idx_img] = p_data->stuffBa.tvl2_occ.u2y[idx_par];

        // v1, v2 (auxiliar minimization variables)
        stuffGo->tvl2_occ.v1[idx_img] = p_data->stuffGo.tvl2_occ.v1[idx_par];
        stuffGo->tvl2_occ.v2[idx_img] = p_data->stuffGo.tvl2_occ.v2[idx_par];

        stuffBa->tvl2_occ.v1[idx_img] = p_data->stuffBa.tvl2_occ.v1[idx_par];
        stuffBa->tvl2_occ.v2[idx_img] = p_data->stuffBa.tvl2_occ.v2[idx_par];

        // Auxiliary (gradients, weighted gradients, divergence, ...)
        // rho_c1
        stuffGo->tvl2_occ.rho_c1[idx_img] = p_data->stuffGo.tvl2_occ.rho_c1[idx_par];
        stuffBa->tvl2_occ.rho_c1[idx_img] = p_data->stuffBa.tvl2_occ.rho_c1[idx_par];

        // rho_c_1
        stuffGo->tvl2_occ.rho_c_1[idx_img] = p_data->stuffGo.tvl2_occ.rho_c_1[idx_par];
        stuffBa->tvl2_occ.rho_c_1[idx_img] = p_data->stuffBa.tvl2_occ.rho_c_1[idx_par];

        // grad_1
        stuffGo->tvl2_occ.grad_1[idx_img] = p_data->stuffGo.tvl2_occ.grad_1[idx_par];
        stuffBa->tvl2_occ.grad_1[idx_img] = p_data->stuffBa.tvl2_occ.grad_1[idx_par];

        // grad__1
        stuffGo->tvl2_occ.grad__1[idx_img] = p_data->stuffGo.tvl2_occ.grad__1[idx_par];
        stuffBa->tvl2_occ.grad__1[idx_img] = p_data->stuffBa.tvl2_occ.grad__1[idx_par];

        if (ofGo->params.step_algorithm == GLOBAL_STEP) {
            stuffGo->tvl2_occ.I0x[idx_img] = p_data->stuffGo.tvl2_occ.I0x[idx_par];
            stuffGo->tvl2_occ.I0y[idx_img] = p_data->stuffGo.tvl2_occ.I0y[idx_par];

            stuffBa->tvl2_occ.I0x[idx_img] = p_data->stuffBa.tvl2_occ.I0x[idx_par];
            stuffBa->tvl2_occ.I0y[idx_img] = p_data->stuffBa.tvl2_occ.I0y[idx_par];
        }
        // I1
        stuffGo->tvl2_occ.I1x[idx_img] = p_data->stuffGo.tvl2_occ.I1x[idx_par];
        stuffGo->tvl2_occ.I1y[idx_img] = p_data->stuffGo.tvl2_occ.I1y[idx_par];

        stuffBa->tvl2_occ.I1x[idx_img] = p_data->stuffBa.tvl2_occ.I1x[idx_par];
        stuffBa->tvl2_occ.I1y[idx_img] = p_data->stuffBa.tvl2_occ.I1y[idx_par];

        stuffGo->tvl2_occ.I1wx[idx_img] = p_data->stuffGo.tvl2_occ.I1wx[idx_par];
        stuffGo->tvl2_occ.I1wy[idx_img] = p_data->stuffGo.tvl2_occ.I1wy[idx_par];

        stuffBa->tvl2_occ.I1wx[idx_img] = p_data->stuffBa.tvl2_occ.I1wx[idx_par];
        stuffBa->tvl2_occ.I1wy[idx_img] = p_data->stuffBa.tvl2_occ.I1wy[idx_par];

        // I_1
        stuffGo->tvl2_occ.I_1x[idx_img] = p_data->stuffGo.tvl2_occ.I_1x[idx_par];
        stuffGo->tvl2_occ.I_1y[idx_img] = p_data->stuffGo.tvl2_occ.I_1y[idx_par];

        stuffBa->tvl2_occ.I_1x[idx_img] = p_data->stuffBa.tvl2_occ.I_1x[idx_par];
        stuffBa->tvl2_occ.I_1y[idx_img] = p_data->stuffBa.tvl2_occ.I_1y[idx_par];

        stuffGo->tvl2_occ.I_1wx[idx_img] = p_data->stuffGo.tvl2_occ.I_1wx[idx_par];
        stuffGo->tvl2_occ.I_1wy[idx_img] = p_data->stuffGo.tvl2_occ.I_1wy[idx_par];

        stuffBa->tvl2_occ.I_1wx[idx_img] = p_data->stuffBa.tvl2_occ.I_1wx[idx_par];
        stuffBa->tvl2_occ.I_1wy[idx_img] = p_data->stuffBa.tvl2_occ.I_1wy[idx_par];

        // vi_div1
        stuffGo->tvl2_occ.vi_div1[idx_par] = p_data->stuffGo.tvl2_occ.vi_div1[idx_img];
        stuffBa->tvl2_occ.vi_div1[idx_par] = p_data->stuffBa.tvl2_occ.vi_div1[idx_img];

        // grad_x1
        stuffGo->tvl2_occ.grad_x1[idx_par] = p_data->stuffGo.tvl2_occ.grad_x1[idx_img];
        stuffBa->tvl2_occ.grad_x1[idx_par] = p_data->stuffBa.tvl2_occ.grad_x1[idx_img];

        // grad_y1
        stuffGo->tvl2_occ.grad_y1[idx_par] = p_data->stuffGo.tvl2_occ.grad_y1[idx_img];
        stuffBa->tvl2_occ.grad_y1[idx_par] = p_data->stuffBa.tvl2_occ.grad_y1[idx_img];

        // vi_div2
        stuffGo->tvl2_occ.vi_div2[idx_par] = p_data->stuffGo.tvl2_occ.vi_div2[idx_img];
        stuffBa->tvl2_occ.vi_div2[idx_par] = p_data->stuffBa.tvl2_occ.vi_div2[idx_img];

        // grad_x2
        stuffGo->tvl2_occ.grad_x2[idx_par] = p_data->stuffGo.tvl2_occ.grad_x2[idx_img];
        stuffBa->tvl2_occ.grad_x2[idx_par] = p_data->stuffBa.tvl2_occ.grad_x2[idx_img];

        // grad_y2
        stuffGo->tvl2_occ.grad_y2[idx_par] = p_data->stuffGo.tvl2_occ.grad_y2[idx_img];
        stuffBa->tvl2_occ.grad_y2[idx_par] = p_data->stuffBa.tvl2_occ.grad_y2[idx_img];

        // g_xi
        stuffGo->tvl2_occ.g_xi11[idx_img] = p_data->stuffGo.tvl2_occ.g_xi11[idx_par];
        stuffGo->tvl2_occ.g_xi12[idx_img] = p_data->stuffGo.tvl2_occ.g_xi12[idx_par];

        stuffBa->tvl2_occ.g_xi11[idx_img] = p_data->stuffBa.tvl2_occ.g_xi11[idx_par];
        stuffBa->tvl2_occ.g_xi12[idx_img] = p_data->stuffBa.tvl2_occ.g_xi12[idx_par];

        stuffGo->tvl2_occ.g_xi21[idx_img] = p_data->stuffGo.tvl2_occ.g_xi21[idx_par];
        stuffGo->tvl2_occ.g_xi22[idx_img] = p_data->stuffGo.tvl2_occ.g_xi22[idx_par];

        stuffBa->tvl2_occ.g_xi21[idx_img] = p_data->stuffBa.tvl2_occ.g_xi21[idx_par];
        stuffBa->tvl2_occ.g_xi22[idx_img] = p_data->stuffBa.tvl2_occ.g_xi22[idx_par];

        // div_g_xi
        stuffGo->tvl2_occ.div_g_xi1[idx_img] = p_data->stuffGo.tvl2_occ.div_g_xi1[idx_par];
        stuffGo->tvl2_occ.div_g_xi2[idx_img] = p_data->stuffGo.tvl2_occ.div_g_xi2[idx_par];

        stuffBa->tvl2_occ.div_g_xi1[idx_img] = p_data->stuffBa.tvl2_occ.div_g_xi1[idx_par];
        stuffBa->tvl2_occ.div_g_xi2[idx_img] = p_data->stuffBa.tvl2_occ.div_g_xi2[idx_par];

        // eta
        stuffGo->tvl2_occ.eta1[idx_img] = p_data->stuffGo.tvl2_occ.eta1[idx_par];
        stuffGo->tvl2_occ.eta2[idx_img] = p_data->stuffGo.tvl2_occ.eta2[idx_par];

        stuffBa->tvl2_occ.eta1[idx_img] = p_data->stuffBa.tvl2_occ.eta1[idx_par];
        stuffBa->tvl2_occ.eta2[idx_img] = p_data->stuffBa.tvl2_occ.eta2[idx_par];

        // F, G
        stuffGo->tvl2_occ.F[idx_img] = p_data->stuffGo.tvl2_occ.F[idx_par];
        stuffGo->tvl2_occ.G[idx_img] = p_data->stuffGo.tvl2_occ.G[idx_par];

        stuffBa->tvl2_occ.F[idx_img] = p_data->stuffBa.tvl2_occ.F[idx_par];
        stuffBa->tvl2_occ.G[idx_img] = p_data->stuffBa.tvl2_occ.G[idx_par];

        // div_u
        stuffGo->tvl2_occ.div_u[idx_par] = p_data->stuffGo.tvl2_occ.div_u[idx_img];
        stuffBa->tvl2_occ.div_u[idx_par] = p_data->stuffBa.tvl2_occ.div_u[idx_img];

        // g_eta
        stuffGo->tvl2_occ.g_eta1[idx_img] = p_data->stuffGo.tvl2_occ.g_eta1[idx_par];
        stuffGo->tvl2_occ.g_eta2[idx_img] = p_data->stuffGo.tvl2_occ.g_eta2[idx_par];

        stuffBa->tvl2_occ.g_eta1[idx_img] = p_data->stuffBa.tvl2_occ.g_eta1[idx_par];
        stuffBa->tvl2_occ.g_eta2[idx_img] = p_data->stuffBa.tvl2_occ.g_eta2[idx_par];

        // div_g_eta
        stuffGo->tvl2_occ.div_g_eta[idx_par] = p_data->stuffGo.tvl2_occ.div_g_eta[idx_img];
        stuffBa->tvl2_occ.div_g_eta[idx_par] = p_data->stuffBa.tvl2_occ.div_g_eta[idx_img];
    }
}


/**
 * @brief                   calls the update functions for the structs of the chosen functional
 *
 * @param ofGo              forward OpticalFlowData struct
 * @param ofBa              backward OpticalFlowData struct
 * @param stuffGo           forward flow' SpecificOFStuff struct with auxiliar variables
 * @param stuffBa           backward flow' SpecificOFStuff struct with auxiliar variables
 * @param p_data            struct of partition data for the current partition
 * @param idx_img           corresponding index in the image-wise 'domain'
 * @param idx_par           corresponding index in the partition-wise 'domain'
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      image_to_partitions
 */
void update_partitions_structures(OpticalFlowData *& ofGo, OpticalFlowData *& ofBa, SpecificOFStuff *& stuffGo,
                                  SpecificOFStuff *& stuffBa, PartitionData *& p_data, const int n_ch,
                                  const int idx_img, const int idx_par, bool img_to_part)
{
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
            case M_TVL1_W:          // TV-l1 coupled with weights
                update_tvl1w_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_NLTVCSAD_W:      // NLTV-CSAD with weights
                update_nltvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_NLTVL1_W:        // NLTV-L1 with weights
                update_nltvl1w_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVCSAD_W:        // TV-CSAD with weights
                update_tvcsadw_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            case M_TVL1_OCC:        // TV-l1 with occlusion
                update_tvl2occ_stuffof(ofGo, stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
            default:                // TV-l1 coupled
                update_tvl2_stuffof(stuffGo, stuffBa, p_data, idx_img, idx_par, img_to_part);
                break;
        }
    } else {
        update_of_data(ofGo, ofBa, p_data, idx_img, idx_par, n_ch, img_to_part);
    }
}


/**
 * @brief                   updates the forward and backward priority queues used in the local growing to store candidates
 *
 * @param queue_Go          forward priority queue
 * @param queue_Ba          backward priority queue
 * @param n_partitions      number of partitions
 * @param p_data            vector of partition data structs containing all the variables for all partitions
 *
 * @sa                      image_to_partitions
 */
void update_candidate_queues(pq_cand &queue_Go, pq_cand &queue_Ba, const int n_partitions,
                             std::vector<PartitionData*> *p_data)
{
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
            // Otherwise, check following partition
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


/**
 * @brief                   updates the partition structures with the latest image-wise structures info (and viceversa)
 *
 * @param oft0              image-wise forward optical flow field (dim: img_width x img_height x 2)
 * @param oft1              image-wise backward optical flow field (dim: img_width x img_height x 2)
 * @param ene_Go            image-wise forward energy (dim: img_width x img_height)
 * @param ene_Ba            image-wise backward energy (dim: img_width x img_height)
 * @param occ_Go            image-wise forward occlusion map (dim: img_width x img_height)
 * @param occ_Ba            image-wise backward occlusion map (dim: img_width x img_height)
 * @param ofGo              image-wise forward OpticalFlowData struct
 * @param ofBa              image-wise backward OpticalFlowData struct
 * @param stuffGo           image-wise forward SpecificOFStuff struct
 * @param stuffBa           image-wise backward SpecificOFStuff struct
 * @param queue_Go          image-wise forward priority queue (candidates)
 * @param queue_Ba          image-wise backward priority queue (candidates)
 * @param n_partitions      number of partitions used
 * @param w_src             image width (refered to as above: 'img_width')
 * @param h_src             image height (refered to as above: 'img_height')
 * @param p_data            vector of data structs (one per partition) that define all needed parameters (see PartitionData)
 * @param img_to_part       direction of the update. If 'true' img values are copied to partition ('false': viceversa)
 *
 * @sa                      update_partitions_structures, update_candidate_queues
 */
void image_to_partitions(float *oft0, float *oft1, float *ene_Go, float *ene_Ba, float *occ_Go, float *occ_Ba,
                         OpticalFlowData *ofGo, OpticalFlowData *ofBa, SpecificOFStuff *stuffGo,
                         SpecificOFStuff *stuffBa, pq_cand &queue_Go, pq_cand &queue_Ba, const int n_partitions,
                         const int w_src, const int h_src, std::vector<PartitionData*> *p_data, bool img_to_part)
{
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
                    // NOTE: in this code, multiple channels are stacked one after the other without alternate values
                    // So, we have: V1_C1, V2_C1, ..., VN_C1, V1_C2, V2_C2, ..., VN_C2
                    int m = j * p_data->at(p)->width + i + k * p_data->at(p)->width * p_data->at(p)->height;
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

        if (img_to_part) {
            p_data->at(p)->ene_Go = ene_Go_p;
            p_data->at(p)->ene_Ba = ene_Ba_p;
            p_data->at(p)->occ_Go = occ_Go_p;
            p_data->at(p)->occ_Ba = occ_Ba_p;
            p_data->at(p)->oft0 = oft0_p;
            p_data->at(p)->oft1 = oft1_p;
            p_data->at(p)->ofGo.params = ofGo->params;
            p_data->at(p)->ofBa.params = ofBa->params;

        } else {
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


/**
 * @brief                   checks if any of the partition-specific queues is empty
 *
 * @param p_data            vector of data structs (one per partition) that define all needed parameters (see PartitionData)
 * @param n_partitions      number of partitions used
 * @return                  'true' if any of the queues is empty, 'false' if none of them is empty
 *
 */
bool anyEmptyQueues(std::vector<PartitionData*> *p_data, const int n_partitions)
{
    int n_empty_queues = 0;
    for (unsigned n = 0; n < n_partitions; n++) {
        if (p_data->at(n)->queue_Go_size == 0) {
            n_empty_queues++;
        }

        if (p_data->at(n)->queue_Ba_size == 0) {
            n_empty_queues++;
        }

        if (n_empty_queues > 0) {
            // Enough to revert back to whole-image
            return true;
        }
        // Otherwise we keep checking other partitions
    }

    return 0 != n_empty_queues;  // Returns 'false' if 'n_emtpy_queues' = 0; 'true' otherwise
}


#endif
