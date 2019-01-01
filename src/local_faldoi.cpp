// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// Copyright (C) 2018, Ferran Pérez <fperez.gamonal@gmail.com>
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
#include "aux_partitions.h"

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


/**
 * @brief       Returns a sample of a 3D array with dimensions w x h x pd (with boundary checking)
 *
 * @param x     input array, the one to be sampled at the specified location
 * @param w     width of the input array 'x'
 * @param h     height of the input array 'x'
 * @param pd    number of channels (depth) of the input array 'x'
 * @param i     column' sample index corresponding to the first dimension
 * @param j     row' sample index corresponding to the second dimension
 * @param l     channel' sample index corresponding to the third dimension
 * @return      returns the array value sampled at the targeted coordinate
 */
static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return INFINITY;
    return x[(i + j * w) * pd + l];
}


/**
 * @brief       Checks if the patch values are too uniform to pass the consistency check
 *
 * @param a     the input array with the patch's values
 * @param tol   the threshold that defines the uniformity
 * @param i     the column index of the pixel that is being analysed
 * @param j     the row index of the pixel that is being analysed
 * @param w     width of the image where 'a' is sampled from
 * @param h     height of the image where 'a' is sampled from
 * @param pd    depth of the image where 'a' is sampled from
 * @return      returns '1' if the patch values are too uniform; '0' otherwise.
 */
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


/**
 * @brief               Checks if the areas within the frames are too uniform or not
 *
 * @param a             source frame (fwd: 'i0' at time 't', bwd: 'i1' at time 't+1')
 * @param b             second frame (fwd: 'i1' at time 't+1', bwd: 'i0' at time 't')
 * @param in0           input flow vector
 * @param trust_in0     array to be filled with 1's or 0's depending on the return value of 'too_uniform'
 * @param w             width of the input frames
 * @param h             height of the input frames
 * @param tol           threshold that defines the uniformity
 *
 * @sa                  too_uniform, pruning_method
 */
void too_uniform_areas(float *a, float *b, float *in0, int *trust_in0, int w, int h, float tol)
{
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


/**
 * @brief               Check forward-backward consistency check for optical flow |u(x) + v(x+u(x))| < eps.
 * @details             Energy map related to the pixels that do not pass the check are put to INFINITY.
 *
 * @param in0           array containing the forward flow (in time)
 * @param in1           array containing the backward flow (in time)
 * @param trust_in0     array to be filled with 1's or 0's depending on whether the flow passes the check or not
 * @param w             width of the input flow fields
 * @param h             height of the input flow fields
 * @param epsilon       threshold that defines the maximum difference between backward and forward flows (consistency)
 *
 * @sa                  pruning_method
 */
void fb_consistency_check(float *in0, float *in1, int *trust_in0, int w, int h, float epsilon)
{
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


/**
 * @brief               Calls the selected pruning method(s) ('fb-consistency' and/or 'too_uniform')
 *
 * @param i0            source frame (fwd: 'i0' at time 't', bwd: 'i1' at time 't+1')
 * @param i1            second frame (fwd: 'i1' at time 't+1', bwd: 'i0' at time 't')
 * @param w             width of the input frames
 * @param h             height of the input frames
 * @param tol           array that contains the tolerances for both pruning methods
 * @param method        defines which method(s) will be used; method[0]=1: 'fb-consistency', method[1]=1: 'too_uniform'
 * @param trust_Go      array to be filled with the pruning decisions for each pixel (forward flow)
 * @param go            forward flow field
 * @param trust_Ba      array to be filled with the pruning decisions for each pixel (backward flow)
 * @param ba            backward flow field
 *
 * @sa                  fb_consistency_check, too_uniform_areas
 */
void pruning_method(float *i0, float *i1, int w, int h, float *tol, const int *method, int *trust_Go, float *go,
                    int *trust_Ba, float *ba)
{
    int *go_fb_check   = new int[w*h];
    int *go_cons_check = new int[w*h];
    int *ba_fb_check   = new int[w*h];
    int *ba_cons_check = new int[w*h];

    for (int i = 0; i < w*h; i++)
    {
        //0 - Invalid pixel 1 - Trustable pixel.
        trust_Go[i] = 1;
        trust_Ba[i] = 1;
    }


    //FB - consistency check
    if (method[0]==1)
    {
        std::printf("FB-Consistency:%f\n",tol[0]);
        fb_consistency_check(go, ba, go_fb_check, w, h, tol[0]);
        fb_consistency_check(ba, go, ba_fb_check, w, h, tol[0]);
    }
    //Too-uniform consistency check
    if (method[1]==1)
    {
        std::printf("Too Uniform -Consistency:%f\n",tol[1]);
        too_uniform_areas(i0, i1, go, go_cons_check, w, h, tol[1]);
        too_uniform_areas(i0, i1, ba, ba_cons_check, w, h, tol[1]);
    }
    for (int i = 0; i < w*h; i++){
        if (method[0] == 1)
        {
            //FB-Consistency
            if (go_fb_check[i] == 0)
            {
                trust_Go[i] = 0;
            }
            if (ba_fb_check[i] == 0)
            {
                trust_Ba[i] = 0;
            }
        }
        //Too uniform -Consistency
        if (method[1] == 1)
        {
            if (go_cons_check[i] == 0)
            {
                trust_Go[i] = 0;
            }
            if (ba_cons_check[i] == 0)
            {
                trust_Ba[i] = 0;
            }
        }
    }

    delete [] go_fb_check;
    delete [] go_cons_check;
    delete [] ba_fb_check;
    delete [] ba_cons_check;
}


/**
 * @brief               deletes non-trustable candidates by setting their flow to NAN and its energy to INF
 * @details             uses the pruning decisions returned by 'pruning_method' to choose what candidates to delete
 *
 * @param ofD           optical flow data struct containing the flow fields and other useful information
 * @param in            output optical flow field to be updated accordingly
 * @param ene_val       energy values to be updated accordingly
 * @param w             width of the input flow field
 * @param h             height of the input flow field
 */
void delete_not_trustable_candidates(OpticalFlowData *ofD, float *in, float *ene_val, const int w, const int h)
{
    int *mask = ofD->trust_points;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float *chi = ofD->chi;

    int n = 0;
    for (int i = 0; i < w*h; i++)
    {
        if (mask[i] == 0)
        {
            //printf("%f\n", ene_val[i]);
            if (ene_val[i]==0.0)
            {
                n++;
            }
            in[i]       = NAN;
            in[i + w*h] = NAN;
            u1[i]       = NAN;
            u2[i]       = NAN;
            ene_val[i]  = INFINITY;
            // If the flow is non trustable, is considered
            // to be an occlusion
            chi[i] = 1;
        }
    }
    printf("Total_seeds%d\n", n);
}


////////////////////////////////////////////////////////////////////////////////
//////////////////LOCAL INITIALIZATION//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Poisson Interpolation
/**
 * @brief           performs the poisson interpolation for the optical flow values in the 'patch'
 *
 * @param ofD       OpticalFlowData struct containing the updated optical flow fields
 * @param patch     patch that is currently being analysed
 */
void interpolate_poisson(OpticalFlowData *ofD, const PatchIndexes &patch, const int w_src)
{
    int w = patch.ei - patch.ii;
    int h = patch.ej - patch.ij;

    int *fixed_points = ofD->fixed_points;
    int wR = w_src; //ofD->params.w;

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


/**
 * @brief           obtains undefined (not fixed) optical flow values by applying a bilateral filter on the patch
 *
 * @param ofD       OpticalFlowData struct containing the updated optical flow fields
 * @param BiFilt    BilateralFilterData struct containing the values of the filter on the patch
 * @param patch     patch that is currently being analysed
 * @param w         width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @param h         height of the optical flow data being processed (to img_height or partition_height if parallelizing)
 */
void bilateral_filter(OpticalFlowData *ofD, BilateralFilterData *BiFilt, const PatchIndexes &patch, const int w,
                      const int h)
{
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

    // For each pixel in the patch non trustable, do bilateral filtering
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


/**
 * @brief               Insert n_neigh-connected candidates into the priority queue with their energies.
 *
 * @param queue         priority queue with candidates
 * @param ene_val       array storing energy values (updated only if the new energy is better/lower)
 * @param ofD           OpticalFlowData struct containing the updated optical flow fields
 * @param i             column index of the pixel that is being currently processed
 * @param j             row index of the pixel that is being currently processed
 * @param ener_N        auxiliar variable to store the new computed energy (and compare against old one in ene_val)
 * @param w             width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @param h             height of the optical flow data being processed (to img_height or partition_height if parallelizing)
 */
void insert_candidates(pq_cand &queue, float *ene_val, OpticalFlowData *ofD, const int i, const int j,
                       const float ener_N, const int w, const int h)
{
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


/**
 * @brief           returns the relative weights corresponding to the current index i, j
 *
 * @param iiw       column weights
 * @param ijw       row weights
 * @param wr        windows radius (patch_size = 2 * wr + 1 in each direction)
 * @param i         current pixel's column index
 * @param j         current pixel's row index
 *
 * @sa              get_index_weight
 */
inline void get_relative_index_weight(int *iiw, int *ijw, const int wr, const int i, const int j)
{
    (*iiw) = (((i - wr) < 0) ? -(i - wr) : 0);
    (*ijw) = (((j - wr) < 0) ? -(j - wr) : 0);
    assert(*iiw >= 0);
    assert(*ijw >= 0);
}


/**
 * @brief           save the index weights generated by 'get_relative_index_weight' to the chosen functional struct
 *
 * @param method    functional chosen
 * @param ofS       SpecificOFStuff struct where the weights should be stored
 * @param wr        windows radius (patch_size = 2 * wr + 1 in each direction)
 * @param i         current pixel's column index
 * @param j         current pixel's row index
 *
 * @sa              get_relative_index_weight
 */
static void get_index_weight(int method, SpecificOFStuff *ofS, const int wr, int i, int j)
{
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


/**
 * @brief           Copy over ofD->u1 and ofD->u2 the presented values in out.
 *
 * @param ofD       OpticalFlowData struct containing the updated optical flow fields
 * @param out       output flow array with original values
 * @param index     indices of the patch that is currently being processed
 * @param w         width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @param h         height of the optical flow data being processed (to img_height or partition_height if parallelizing)
 */
inline void copy_fixed_coordinates(OpticalFlowData *ofD, const float *out, const PatchIndexes &index, const int w,
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


/**
 * @brief           Check if there is at least one pixel that hasn't survived to the prunning.
 *
 * @param ofD       OpticalFlowData struct containing the updated optical flow fields
 * @param index     indices of the patch that is currently being processed
 * @param w         width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @return          1's if the pixel's flow is trustable, 0's otherwise
 */
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


/**
 * @brief               adds neigbours of the pixel that is currently being processed to the queue
 * @details             computes the new energy for the patch and then inserts the candidates that lowered their energy
 *
 * @param i0            source frame at time 't'
 * @param i1            second frame at time 't+1'
 * @param i_1           previous frame at time 't-1' (used for occlusions only)
 * @param ene_val       array that stores the energy values updated (maybe) during the previous candidate's iteration
 * @param ofD           OpticalFlowData struct containing the updated optical flow fields
 * @param ofS           SpecificOFStuff struct where functional-specific variables reside
 * @param queue         priority queue that contains the candidates
 * @param i             column index of the current pixel being processed
 * @param j             row index of the current pixel being processed
 * @param iteration     iteration index of the local minimization step
 * @param out           array that stores the output flow (maybe) updated during the previous candidate's iteration
 * @param out_occ       array that stores the occlusions' map (maybe) updated during the previous candidate's iteration
 * @param BiFilt        struct that contains the indices and weights of the bilateral filter
 * @param w             width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @param h             height of the optical flow data being processed (to img_height or partition_height if parallelizing)
 */
static void add_neighbors(const float *i0, const float *i1, const float *i_1, float *ene_val, OpticalFlowData *ofD,
                          SpecificOFStuff *ofS, pq_cand *queue, const int i, const int j, const int iteration,
                          float *out, float *out_occ, BilateralFilterData *BiFilt, const int w, const int h)
{
    const int wr = ofD->params.w_radio;
    float ener_N;
    const PatchIndexes index = get_index_patch(wr, w, h, i, j, 1);
    int method = ofD->params.val_method; // used to include no occ.

    // In first iteration, Poisson interpolation
    if (iteration == 0) {
        // Interpolate by poisson on initialization
        copy_fixed_coordinates(ofD, out, index, w, h);
        // Poisson Interpolation (4wr x 4wr + 1)
        interpolate_poisson(ofD, index, w);

    } else if (check_trustable_patch(ofD, index, w) == 0) {
        // Interpolate by bilateral filtering if some points do not survive to prunning
        //if (check_trustable_patch(ofD, index, w) == 0) {

        copy_fixed_coordinates(ofD, out, index, w, h);
        interpolate_poisson(ofD, index, w);
        // TODO: fix bilateral filter as it is quite faster
        //bilateral_filter(ofD, BiFilt, index, w, h);  // yields a far worse estimation
        // check implementation details
        // }
    }

    // get index's weight (if the functional uses them)
    get_index_weight(method, ofS, wr, i, j);


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


/**
 * @brief               inserts the initial seeds to the priority queue by using initial flow derived from the sparse matches (SIFT or deepmatching)
 * @details             initialises to default values: flow to NAN, energy to INF and occlusions to 0 (if it applies)
 *
 * @param i0            source frame at time 't'
 * @param i1            second frame at time 't+1'
 * @param i_1           previous frame at time 't-1' (used for occlusions only)
 * @param in_flow       array that contains the initial flow obtained from the sparse matches
 * @param queue         priority queue where the initial candidates will be inserted
 * @param ofD           OpticalFlowData struct with default values
 * @param ofS           SpecificOFStuff struct with default values
 * @param ene_val       array that will store the energy of all pixels in the image
 * @param out_flow      array that will store the optical flow values of all pixels in the image
 * @param out_occ       array that will store the occlusion map for the image
 * @param BiFilt        struct that contains the indices and weights of the bilateral filter
 * @param w             width of the optical flow data being processed
 * @param h             height of the optical flow data being processed
 */
void insert_initial_seeds(const float *i0, const float *i1, const float *i_1, float *in_flow, pq_cand *queue,
                          OpticalFlowData *ofD, SpecificOFStuff *ofS, float *ene_val, float *out_flow, float *out_occ,
                          BilateralFilterData *BiFilt, const int w, const int h)
{
    int wr = ofD->params.w_radio;

    //Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++)
    {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out_flow[i] = NAN;
        out_flow[w*h + i] = NAN;
    }


    ofD->params.w_radio = 1;
    //Fixed the initial seeds.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            //Indicates the initial seed in the similarity map
            if (std::isfinite(in_flow[j*w +i]) && std::isfinite(in_flow[w*h + j*w +i]))
            {
                out_flow[j*w + i] = in_flow[j*w + i];
                out_flow[w*h + j*w + i] = in_flow[w*h + j*w +i];
                ofD->fixed_points[j*w + i] = 1;
                // add_neigbors 0 means that during the propagation interpolates the patch
                // based on the energy.
                add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, 0, out_flow, out_occ, BiFilt, w, h);
                out_flow[j*w + i] = NAN;
                out_flow[w*h + j*w + i] = NAN;
                ofD->fixed_points[j*w + i] = 0;
            }
        }
    ofD->params.w_radio = wr;
    //Propagate the information of the initial seeds to their neighbours.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            if (std::isfinite(in_flow[j*w +i]) && std::isfinite(in_flow[w*h + j*w +i]))
            {
                out_flow[j*w + i] = in_flow[j*w + i];
                out_flow[w*h + j*w + i] = in_flow[w*h + j*w +i];
                ofD->fixed_points[j*w + i] = 1;
                ene_val[j*w + i] = 0.0;
            }
        }
}


// Insert each pixel into the queue as possible candidate. Its related energy comes
// from the energy stored at the moment that the pixel was fixed.
/**
 * @brief           Insert each pixel into the queue as possible candidate with its energy.
 * @details         its related energy comes from the energy stored at the moment that the pixel was fixed
 *
 * @param in        array containing the flow field values
 * @param ofD       OpticalFlowData struct containing the updated optical flow fields
 * @param queue     priority queue where the new candidates will be inserted
 * @param ene_val   array that stores the energy values updated (maybe) during previous iterations
 * @param out_occ   array that stores the occlusions' map (maybe) updated during previous iterations
 * @param w         width of the optical flow data being processed
 * @param h         height of the optical flow data being processed
 */
void insert_potential_candidates(const float *in, OpticalFlowData *ofD, pq_cand &queue, float *ene_val,
                                 float *out_flow, const float *out_occ, const int w, const int h)
{
    //Fixed the initial seeds.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            //Indicates the initial seed in the similarity map
            if (std::isfinite(in[j*w +i]) && std::isfinite(in[w*h + j*w +i]))
            {
                SparseOF element;
                element.i = i; // column
                element.j = j; // row
                element.u = in[j*w +i];
                element.v = in[w*h + j*w +i];
                //Obs: Notice that ene_val contains (en)*saliency
                element.sim_node = ene_val[j*w +i];
                if (ofD->params.val_method >= 8) {
                    element.occluded = out_occ[j * w + i];
                }
                assert(std::isfinite(ene_val[j*w +i]));
                queue.push(element);
            }
        }

    //Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++)
    {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out_flow[i] = NAN;
        out_flow[w*h + i] = NAN;
    }
}


// Initialize the data to prepare everything for the region growing
/**
 * @brief               Initialize the data to prepare everything for the region growing
 * @details             energy, flow and fixed (computed or not) variables set to default values (Inf, NAN and 0)
 *
 * @param ofD           OpticalFlowData struct containing fixed/not fixed variable to be reset to default
 * @param ene_val       array that stores the energy values to be reset to default
 * @param out           output flow array to be reset to default
 * @param w             width of the optical flow data being processed
 * @param h             height of the optical flow data being processed
 */
void prepare_data_for_growing(OpticalFlowData *ofD, float *ene_val, float *out, const int w, const int h)
{
    // Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++)
    {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out[i] = NAN;
        out[w*h + i] = NAN;
    }
}


/**
 * @brief               function that manages a specific iteration of the local minimization, processing all the queue's candidates
 *
 * @param i0            source frame at time 't'
 * @param i1            second frame at time 't+1'
 * @param i_1           previous frame at time 't-1' (used for occlusions only)
 * @param queue         priority queue from/to which candidates will be obtained/inserted
 * @param ofS           SpecificOFStuff struct where functional-specific variables reside
 * @param ofD           OpticalFlowData struct containing the updated optical flow fields
 * @param iteration     index for the local minimization's current iteration
 * @param ene_val       array that stores the energy values (will be updated through the execution of this function)
 * @param out_flow      array that stores the optical flow fields (will be updated through the execution of this function)
 * @param out_occ       array that stores the occlusions map (may be updated if occlusions are estimated)
 * @param BiFilt        struct that contains the indices and weights of the bilateral filter
 * @param fwd_or_bwd    boolean that defines if we are processing a forward or backward flow (if we store partial results)
 * @param w             width of the optical flow data being processed (to img_width or partition_width if parallelizing)
 * @param h             height of the optical flow data being processed (to img_height or partition_height if parallelizing)
 */
void local_growing(const float *i0, const float *i1, const float *i_1, pq_cand *queue, SpecificOFStuff *ofS,
                   OpticalFlowData *ofD, int iteration, float *ene_val, float *out_flow, float *out_occ,
                   BilateralFilterData *BiFilt, bool fwd_or_bwd, const int w, const int h, const int part_idx)
{
    std::vector<int> percent_print = {30, 70, 80, 95, 100};
    int fixed = 0;
    const int size = w * h;
    std::printf("queue size at start = %d\n", (int) queue->size());
    while (!queue->empty()) {
        //std::printf("Fixed elements = %d\n", val);
        SparseOF element = queue->top();
        int i = element.i;
        int j = element.j;

        queue->pop();
        if (!ofD->fixed_points[j * w + i]) {
            assert(std::isfinite(element.sim_node));
            float u = element.u;
            float v = element.v;
            float energy = element.sim_node;
            float occlusion;

            if (ofD->params.val_method >= 8) {
                occlusion = element.occluded;
            } else {
                occlusion = 0.0;
            }

            if (!std::isfinite(u)) {
                std::printf("U1 = %f\n", u);
            }
            if (!std::isfinite(v)) {
                std::printf("U2 = %f\n", v);
            }

            ofD->fixed_points[j * w + i] = 1;
            fixed ++;

            out_flow[j * w + i] = u;
            out_flow[w * h + j * w + i] = v;
            ene_val[j * w + i] = energy;
            out_occ[j * w + i] = occlusion;
            // TODO: copy the values so they are taken into account in the minimization
            // ofD->u1[j*w + i] = u;
            // ofD->u2[j*w + i] = v;

            add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, iteration, out_flow, out_occ, BiFilt, w, h);

            // From here to the end of the function:
            // Code used to print partial growing results for debugging or further exploration
            // Just set add the flag '-partial_res 1' when you call any of the Python scripts (or the binary)
            float percent = 100 * fixed * 1.0 / size * 1.0;

            if (ofD->params.part_res == 1) {
                for (int k = 0; k < 4; k++) {
                    if (percent > percent_print[k] && percent < percent_print[k + 1]) {
                        string filename_flow = " ";
                        string filename_occ = " ";
                        if (fwd_or_bwd) {
                            if (ofD->params.split_img) {
                                filename_flow =
                                        "../Results/Partial_results/partial_results_fwd_" +
                                        std::to_string(percent_print[k]) +
                                        "_iter_" + std::to_string(iteration) + "_part_idx" + to_string(part_idx) + ".flo";
                                iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);

                                if (ofD->params.val_method >= 8) {
                                    filename_occ =
                                            "../Results/Partial_results/partial_results_fwd_" +
                                            std::to_string(percent_print[k]) +
                                            "_iter_" + std::to_string(iteration) + "_part_idx" + to_string(part_idx) +
                                            "_occ.png";
                                    auto *out_occ_int = new int[w * h];
                                    for (int l = 0; l < w * h; l++) {
                                        out_occ_int[l] = out_occ[l];
                                    }
                                    iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);
                                }

                            } else {
                                filename_flow =
                                        "../Results/Partial_results/partial_results_fwd_" +
                                        std::to_string(percent_print[k]) +
                                        "_iter_" + std::to_string(iteration) + ".flo";
                                iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);

                                if (ofD->params.val_method >= 8) {
                                    filename_occ =
                                            "../Results/Partial_results/partial_results_fwd_" +
                                            std::to_string(percent_print[k]) +
                                            "_iter_" + std::to_string(iteration) + "_occ.png";
                                    auto *out_occ_int = new int[w * h];
                                    for (int l = 0; l < w * h; l++) {
                                        out_occ_int[l] = out_occ[l];
                                    }
                                    iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);
                                }

                            }

                            percent_print[k] = 200;
                        }
                    }
                }
            }
        }
    }
    if (ofD->params.part_res == 1) {
        if (fwd_or_bwd) {
            string filename_flow = " ";
            string filename_occ = " ";
            if (ofD->params.split_img) {

                filename_flow =
                        "../Results/Partial_results/partial_results_fwd_100_iter_" + std::to_string(iteration)
                        + "_part_idx" + to_string(part_idx) + ".flo";
                iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);

                if (ofD->params.val_method >= 8) {
                filename_occ =
                        "../Results/Partial_results/partial_results_fwd_100_iter_" +
                        std::to_string(iteration) + "_part_idx" + to_string(part_idx) + "_occ.png";
                auto *out_occ_int = new int[w * h];
                for (int i = 0; i < w * h; i++) {
                    out_occ_int[i] = out_occ[i];
                }
                iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);
            }
            } else {
                filename_flow =
                        "../Results/Partial_results/partial_results_fwd_100_iter_" + std::to_string(iteration) + ".flo";
                iio_save_image_float_split(filename_flow.c_str(), out_flow, w, h, 2);

                if (ofD->params.val_method >= 8) {
                    filename_occ =
                            "../Results/Partial_results/partial_results_fwd_100_iter_" + std::to_string(iteration) +
                            "_occ.png";
                    auto *out_occ_int = new int[w * h];
                    for (int i = 0; i < w * h; i++) {
                        out_occ_int[i] = out_occ[i];
                    }
                    iio_save_image_int(filename_occ.c_str(), out_occ_int, w, h);
                }
            }
        }
    }


}


/**
 * @brief               manages the whole local minimization, calling 'local_growing' for each iteration and updating all variables
 * @details             every iteration involves: growing, pruning, deleting non valid candidates and updating queues
 *
 * @param go            initial forward flow field obtained from the sparse matches
 * @param ba            initial backward flow field obtained from the sparse matches
 * @param i0n           normalised (gray and smooth) source frame
 * @param i1n           normalised (gray and smooth) second frame
 * @param i_1n          normalised (gray and smooth) previous frame
 * @param i2n           normalised (gray and smooth) second frame
 * @param sal_go        forward saliency array
 * @param sal_ba        backward saliency array
 * @param params        basic parameters to check input arguments values(e.g.: functional, windows radius, etc)
 * @param ene_val       array that stores the best energy found in the local minimization (constantly updated)
 * @param out_flow      array that stores the best optical flow values in the local minimization (constantly updated)
 * @param out_occ       array that stores the occlusion map in the local minimization (may be updated if it applies)
 */
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
            insert_initial_seeds(i0n, i1n, i_1n, go, &queue_Go, &ofGo, &stuffGo, ene_Go, oft0, occ_Go, BiFilt_Go, w, h);
        } else {
            insert_initial_seeds(i1n, i0n, i2n, ba, &queue_Ba, &ofBa, &stuffBa, ene_Ba, oft1, occ_Ba, BiFilt_Ba, w, h);
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
    float tol[2] = {params.epsilon, TU_TOL};
    int p[2] = {1, 0};

    // Create partitions data structures
    std::vector<PartitionData*> p_data;
    std::vector<PartitionData*> p_data_r;

    if (params.split_img == 1) {
        auto clk_init_part = system_clock::now(); // PROFILING
        // Initialise partitions
        // Note: to avoid reinforcing discontinuities that may be caused by the partitions, we
        // flip the grid/partition to avoid 'cutting' the image twice on successive iterations on the same spot
        // (i.e.: 3x2 => 2x3)
        // Odd iterations (i == 1, 3, 5, ...) ==> h_parts x v_parts grid
        init_subimage_partitions(i0, i1, i_1, i2, i0n, i1n, i_1n, i2n, BiFilt_Go, BiFilt_Ba, sal_go, sal_ba, w, h,
                                 params.h_parts, params.v_parts, &p_data, params);

        // Even iterations (i == 2, 4, 6, ...) ==> v_parts x h_parts grid
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
                                  true, w, h, -1);
                    auto clk_fwd_grow = system_clock::now();  // PROFILING
                    duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                    cout << "(match growing) FWD local growing (it=" << i << ") took "
                         << elapsed_secs_fwd_grow.count() << endl;
                } else {
                    // BWD
                    auto clk_start_bwd = system_clock::now();
                    local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                  false, w, h, -1);
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
            insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, oft0, occ_Go, w, h);
            insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, oft1, occ_Ba, w, h);

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

        } else if ((i > 0 && i <= iter - 1) && params.split_img == 1) {
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
                                          p_data.at(m)->width, p_data.at(m)->height, m);

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
                                          p_data.at(m)->width, p_data.at(m)->height, m);

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
                    // are not empty but that would imply significant changes to the code as is
                    cout << "Reverted back to whole-image based processing due to lack of seeds (1 or more empty queues)" << endl;
#ifdef _OPENMP
#pragma omp parallel for
                    for (int k = 0; k < 2; k++) {
                        if (k == 0) {
                            // FWD
                            auto clk_start_fwd = system_clock::now();
                            local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, i, ene_Go, oft0, occ_Go, BiFilt_Go,
                                          true, w, h, -1);
                            auto clk_fwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                            cout << "(match growing) FWD local growing (it=" << i << ") took "
                                 << elapsed_secs_fwd_grow.count() << endl;
                        } else {
                            // BWD
                            auto clk_start_bwd = system_clock::now();
                            local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                          false, w, h, -1);
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
                insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, oft0, occ_Go, w, h);
                insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, oft1, occ_Ba, w, h);

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

            } else if (i % 2 == 0 && i <= iter - 1) {  // part. grid: v_parts (cols) x h_parts (rows)
                auto clk_even_start = system_clock::now();  // PROFILING

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
                                          p_data_r.at(m)->width, p_data_r.at(m)->height, m);

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
                                          p_data_r.at(m)->width, p_data_r.at(m)->height, m);

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
                                          true, w, h, -1);
                            auto clk_fwd_grow = system_clock::now();  // PROFILING
                            duration<double> elapsed_secs_fwd_grow = clk_fwd_grow - clk_start_fwd;  // PROFILING
                            cout << "(match growing) FWD local growing (it=" << i << ") took "
                                 << elapsed_secs_fwd_grow.count() << endl;
                        } else {
                            // BWD
                            auto clk_start_bwd = system_clock::now();
                            local_growing(i1n, i0n, i2n, &queue_Ba, &stuffBa, &ofBa, i, ene_Ba, oft1, occ_Ba, BiFilt_Ba,
                                          false, w, h, -1);
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
                insert_potential_candidates(oft0, &ofGo, queue_Go, ene_Go, oft0, occ_Go, w, h);
                insert_potential_candidates(oft1, &ofBa, queue_Ba, ene_Ba, oft1, occ_Ba, w, h);

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
                              p_data.at(m)->width, p_data.at(m)->height, m);

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
            local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, iter, ene_Go, oft0, occ_Go, BiFilt_Go, true, w, h, iter);
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
        local_growing(i0n, i1n, i_1n, &queue_Go, &stuffGo, &ofGo, iter, ene_Go, oft0, occ_Go, BiFilt_Go, true, w, h, iter);

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

/**
 * @brief           main function that reads the command arguments, calls 'match_growing_variational' and frees memory
 *
 * @param argc      argument count (i.e.: how many input arguments have been inputted by the user)
 * @param argv      argument variables (actual values of the inputted arguments)
 * @return          returns 0 if the execution finished without errors, other error codes otherwise
 */
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
	auto fb_threshold = pick_option(args, "fb_thresh", to_string(FB_TOL));		// Threshold for the FB pruning (if tol > thr, discard)
	auto partial_results = pick_option(args, "partial_res",
								       to_string(SAVE_RESULTS));				// Whether to store intermediate flows in "../Results/Partial_results"

    if (args.size() < 6 || args.size() > 9) {
        // Without occlusions
        fprintf(stderr, "Without occlusions (nº of params: 5 or 7 + 1 (own function name)):\n");
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff"
                        " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                        " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                        " [-split_img split_image] [-h_parts horiz_parts]"
                        " [-v_parts vert_parts] [-fb_thresh thresh] [-partial_res val]\n", args.size(), args[0].c_str());
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff sal0.tiff sal1.tiff"
                        " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                        " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                        " [-split_img split_image] [-h_parts horiz_parts]"
                        " [-v_parts vert_parts] [-fb_thresh thresh] [-partial_res val]\n", args.size(), args[0].c_str());
        fprintf(stderr, "\n");
        // With occlusions
        fprintf(stderr, "With occlusions (nº of params: 7 or 9 + 1 (own function name)):\n");
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png"
                        " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                        " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                        " [-split_img split_image] [-h_parts horiz_parts]"
                        " [-v_parts vert_parts] [-fb_thresh thresh] [-partial_res val]\n", args.size(), args[0].c_str());
        fprintf(stderr,
                "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png sal0.tiff sal1.tiff"
                " [-m method_id] [-wr windows_radio] [-p file of parameters]"
                " [-loc_it local_iters] [-max_pch_it max_iters_patch]"
                " [-split_img split_image] [-h_parts horiz_parts]"
                " [-v_parts vert_parts] [-fb_thresh thresh] [-partial_res val]\n", args.size(), args[0].c_str());
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
    //TODO: works with occlusions but we lost 'const' advantages
    string filename_occ;

    // Case w/o occ or saliency has all params already (args.size()==6)
    if (args.size() == 7 || args.size() == 9) // case with occlusions
    {
        filename_occ = args[6];
    }
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
	float fb_thresh = stof(fb_threshold);
	int partial_res = stoi(partial_results);

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
        return fprintf(stderr, "ERROR: input flow field size mismatch (global\n");

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
	params.epsilon = fb_thresh;
	params.part_res = partial_res;
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
