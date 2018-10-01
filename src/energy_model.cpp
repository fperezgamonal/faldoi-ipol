// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// All rights reserved.
#ifndef ENERGY_MODEL
#define ENERGY_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <algorithm>


extern "C" {
#include "bicubic_interpolation.h"
#include "mask.h"
}

#include "utils.h"
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "energy_model.h"

//Models
#include "tvl2_model.h"
#include "tvcsad_model.h"
#include "nltv_model.h"
#include "tvcsad_model.h"
#include "nltvcsad_model.h"

//Models with weights
#include "tvl2w_model.h"
#include "tvcsadw_model.h"
#include "nltvw_model.h"
#include "nltvcsadw_model.h"

//Models with occlusions
#include "tvl2_model_occ.h"


void rgb_to_gray(const float *in, const int w, const int h, float *out) {
    const int size = w * h;

    for (int i = 0; i < size; i++) {

        out[i] = .299 * in[i] + .587 * in[size + i] + .114 * in[2 * size + i];

    }

}

float pow2(const float f) { return f * f; }


// Non-normalized images are assumed
void rgb_to_lab(const float *in, int size, float *out) {

    const float T = 0.008856;
    const float color_attenuation = 1.5f;
    for (int i = 0; i < size; i++) {

        const float r = in[i] / 255.f;
        const float g = in[i + size] / 255.f;
        const float b = in[i + 2 * size] / 255.f;
        float X = 0.412453 * r + 0.357580 * g + 0.180423 * b;
        float Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
        float Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;
        X /= 0.950456;
        Z /= 1.088754;
        float Y3 = pow(Y, 1. / 3);
        float fX = X > T ? pow(X, 1. / 3) : 7.787 * X + 16 / 116.;
        float fY = Y > T ? Y3 : 7.787 * Y + 16 / 116.;
        float fZ = Z > T ? pow(Z, 1. / 3) : 7.787 * Z + 16 / 116.;
        float L = Y > T ? 116 * Y3 - 16.0 : 903.3 * Y;
        float A = 500 * (fX - fY);
        float B = 200 * (fY - fZ);
        // Correct L*a*b*: dark area or light area have less reliable colors
        float correct_lab = exp(-color_attenuation * pow2(pow2(L / 100) - 0.6));
        out[i] = L;
        out[i + size] = A * correct_lab;
        out[i + 2 * size] = B * correct_lab;
    }
}




//////////////////////////////////////
/////Initialize bilteral filter///////
/////////////////////////////////////


static float weight_dist(int x1, int x2,        // center of patch
                         int chi1, int chi2     // pixel in patch
) {
    float sigma = SIGMA_BILATERAL_DIST;

    float result = exp(-0.5 * (pow(x1 - chi1, 2.0) + pow(x2 - chi2, 2.0)) / pow(sigma, 2.0));
    return result;
}

static float weight_color(float color_x,
                          float color_chi) {
    float sigma = SIGMA_BILATERAL_COLOR;

    float result = exp(-0.5 * pow(std::abs(color_x - color_chi) / sigma, 2.0));
    return result;
}

BilateralFilterData *init_weights_bilateral(
        float *i0,
        int w,
        int h) {

    auto *Filter_data = new BilateralFilterData[1];
    Filter_data->weights_filtering = new Weights_Bilateral[w * h];

    Filter_data->indexes_filtering.resize(w * h);
    int wr_filter = PATCH_BILATERAL_FILTER;


    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {

            // For each pixel in the image, compute neighbor of weights
            const int ij = j * w + i;
            auto neighbor = get_index_patch(wr_filter, w, h, i, j, 1);
            Filter_data->indexes_filtering[ij] = neighbor;
            const int w_patch = neighbor.ei - neighbor.ii;
            const int h_patch = neighbor.ej - neighbor.ij;
            Filter_data->weights_filtering[ij].weight = new float[w_patch * h_patch];

            for (int idx_j = 0; idx_j < h_patch; idx_j++) {
                for (int idx_i = 0; idx_i < w_patch; idx_i++) {

                    int x = idx_i + neighbor.ii;
                    int y = idx_j + neighbor.ij;
                    int xy = y * w + x;


                    int idx_ij = idx_j * w_patch + idx_i;

                    float dist = weight_dist(i, j, x, y);
                    float color = weight_color(i0[ij], i0[xy]);

                    // Save the weight for each point in the neighborhood
                    Filter_data->weights_filtering[ij].weight[idx_ij] = color * dist;
                }
            }
        }
    }
    return Filter_data;
}


// Initialization of Optical Flow data for global method
OpticalFlowData init_Optical_Flow_Data(
        const Parameters &params
) {
    int w = params.w;
    int h = params.h;
    OpticalFlowData of{};
    of.u1 = new float[w * h * 2];
    of.u2 = of.u1 + w * h;
    of.u1_ba = new float[w * h * 2];
    of.u2_ba = of.u1_ba + w * h;
    of.chi = new float[w * h];
    of.params = params;
    return of;
}

// Initialization of Optical Flow data for local method
OpticalFlowData init_Optical_Flow_Data(
        float *saliency,
        const Parameters &params,
        const int w,
        const int h
) {
    //int w = params.w;
    //int h = params.h;
    OpticalFlowData of{};
    of.u1 = new float[w * h * 2];
    of.u2 = of.u1 + w * h;
    of.u1_ba = new float[w * h * 2];
    of.u2_ba = of.u1_ba + w * h;
    of.u1_filter = new float[w * h * 2];
    of.u2_filter = of.u1_filter + w * h;
    of.chi = new float[w * h];
    of.fixed_points = new int[w * h];
    of.trust_points = new int[w * h];
    of.saliency = saliency;
    of.params = params;
    of.params.w = w;
    of.params.h = h;
    return of;
}

void initialize_auxiliar_stuff(
        SpecificOFStuff &ofStuff,
        OpticalFlowData &ofCore,
        const int w,
        const int h
) {
    switch (ofCore.params.val_method) {

        case M_NLTVL1:          // NLTV-L1
            initialize_stuff_nltvl1(&ofStuff, &ofCore, w, h);
            break;
        case M_TVCSAD:          // TV-CSAD
            initialize_stuff_tvcsad(&ofStuff, &ofCore, w, h);
            break;
        case M_NLTVCSAD:        // NLTV-CSAD
            initialize_stuff_nltvcsad(&ofStuff, &ofCore, w, h);
            break;
        case M_TVL1_W:          // TV-l2 coupled with weights
            initialize_stuff_tvl2coupled_w(&ofStuff, &ofCore, w, h);
            break;
        case M_NLTVCSAD_W:      // NLTV-CSAD with weights
            initialize_stuff_nltvcsad_w(&ofStuff, &ofCore, w, h);
            break;
        case M_NLTVL1_W:        // NLTV-L1 with weights
            initialize_stuff_nltvl1_w(&ofStuff, &ofCore, w, h);
            break;
        case M_TVCSAD_W:        // TV-CSAD with weights
            initialize_stuff_tvcsad_w(&ofStuff, &ofCore, w, h);
            break;
        case M_TVL1_OCC:        // TV-l2 with occlusion
            initialize_stuff_tvl2coupled_occ(ofStuff, ofCore, w, h);
            break;
        default:                //TV-l2 coupled
            initialize_stuff_tvl2coupled(&ofStuff, &ofCore, w, h);
    }

}

void free_auxiliar_stuff(SpecificOFStuff *ofStuff, OpticalFlowData *ofCore) {

    const int method = ofCore->params.val_method;

    switch (method) {

        case M_NLTVL1:          // NLTVL1
            free_stuff_nltvl1(ofStuff);
            break;
        case M_TVCSAD:          // TV-CSAD
            free_stuff_tvcsad(ofStuff);
            break;
        case M_NLTVCSAD:        // NLTV-CSAD
            free_stuff_nltvcsad(ofStuff);
            break;
        case M_TVL1_W:          // TV-l2 with weights
            free_stuff_tvl2coupled_w(ofStuff);
            break;
        case M_NLTVCSAD_W:      // NLTV-CSAD with weights
            free_stuff_nltvcsad_w(ofStuff);
            break;
        case M_NLTVL1_W:        // NLTV-L1 with weights
            free_stuff_nltvl1_w(ofStuff);
            break;
        case M_TVCSAD_W:        // TV-CSAD with weights
            free_stuff_tvcsad_w(ofStuff);
            break;
        case M_TVL1_OCC:        // TV-l2 with occlusion
            free_stuff_tvl2coupled_occ(ofStuff);
            break;
        default:                // TV-l2 coupled
            free_stuff_tvl2coupled(ofStuff);
    }
}


void prepare_stuff(
        SpecificOFStuff *ofStuff1,
        OpticalFlowData *ofCore1,
        SpecificOFStuff *ofStuff2,
        OpticalFlowData *ofCore2,
        const float *i0,
        const float *i1,
        const float *i_1,
        const float *i2,
        const int pd,
        float **out_i0,
        float **out_i1,
        float **out_i_1,
        float **out_i2,
        const int w,
        const int h
) {

    const int method = ofCore1->params.val_method;

    switch(method)
    {
        case M_NLTVL1:          // NLTV-L1
        {
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];

            auto *alb = new float[w*h*pd];
            auto *blb = new float[w*h*pd];

            int n_d = NL_DUAL_VAR;
            int radius = NL_BETA;

            rgb_to_lab(i0, w*h, alb);
            rgb_to_lab(i0, w*h, blb);
            // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
            nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                    ofStuff1->nltvl1.p, ofStuff1->nltvl1.q);
            nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                    ofStuff2->nltvl1.p, ofStuff2->nltvl1.q);
            if (pd!=1)
            {
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i0, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i0,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->nltvl1.I1x, ofStuff1->nltvl1.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->nltvl1.I1x, ofStuff2->nltvl1.I1y,
                              w, h);

            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

            delete [] alb;
            delete [] blb;
        }
            break;
        case M_TVCSAD:          // TVCSAD
        {
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];
            int rdt = DT_R;
            int ndt = DT_NEI;
            std::printf("1 - Inicializado CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->tvcsad.pnei);
            std::printf("2 - Inicializado CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->tvcsad.pnei);
            if (pd!=1)
            {
                // std::printf("Numero canales:%d\n",pd);
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->tvcsad.I1x, ofStuff1->tvcsad.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->tvcsad.I1x, ofStuff2->tvcsad.I1y,
                              w, h);
            *out_i0 = a_tmp;
            *out_i1 = b_tmp;
            std::printf("Exitting CSAD\n");
        }
            break;
        case M_NLTVCSAD:        // NLTV-CSAD
        {
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];

            auto *alb = new float[w*h*pd];
            auto *blb = new float[w*h*pd];

            int n_d = NL_DUAL_VAR;
            int radius = NL_BETA;
            int rdt = DT_R;
            int ndt = DT_NEI;
            std::printf("Initializing CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->nltvcsad.pnei);
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->nltvcsad.pnei);

            rgb_to_lab(i0, w*h, alb);
            rgb_to_lab(i0, w*h, blb);
            // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
            nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                    ofStuff1->nltvcsad.p, ofStuff1->nltvcsad.q);
            nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                    ofStuff2->nltvcsad.p, ofStuff2->nltvcsad.q);
            std::printf("Initializing NLTV\n");
            if (pd!=1)
            {
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->nltvcsad.I1x, ofStuff1->nltvcsad.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->nltvcsad.I1x, ofStuff2->nltvcsad.I1y,
                              w, h);

            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

            delete [] alb;
            delete [] blb;

            std::printf("Exitting NLTV_CSAD\n");

        }
            break;
        case M_TVL1_W:      // TV-l2 coupled with weights
        {

            float *weight1 = ofStuff1->tvl2w.weight;
            float *weight2 = ofStuff2->tvl2w.weight;
            gaussian1Dweight(weight1, ofCore1->params.w_radio);
            gaussian1Dweight(weight2, ofCore2->params.w_radio);

            std::printf("Weights\n");
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];
            if (pd!=1)
            {
                // std::printf("Number of channels:%d\n",pd);
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->tvl2w.I1x, ofStuff1->tvl2w.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->tvl2w.I1x, ofStuff2->tvl2w.I1y,
                              w, h);
            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

        }
            break;
        case M_NLTVCSAD_W:      // NLTV-CSAD with weights
        {

            float *weight1 = ofStuff1->nltvcsadw.weight;
            float *weight2 = ofStuff2->nltvcsadw.weight;
            gaussian1Dweight(weight1, ofCore1->params.w_radio);
            gaussian1Dweight(weight2, ofCore2->params.w_radio);

            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];

            auto *alb = new float[w*h*pd];
            auto *blb = new float[w*h*pd];

            int n_d = NL_DUAL_VAR;
            int radius = NL_BETA;
            int rdt = DT_R;
            int ndt = DT_NEI;
            std::printf("Initializing CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->nltvcsadw.pnei);
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->nltvcsadw.pnei);

            rgb_to_lab(i0, w*h, alb);
            rgb_to_lab(i0, w*h, blb);
            // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
            nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                    ofStuff1->nltvcsadw.p, ofStuff1->nltvcsadw.q);
            nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                    ofStuff2->nltvcsadw.p, ofStuff2->nltvcsadw.q);
            std::printf("Initializing NLTV\n");
            if (pd!=1)
            {
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->nltvcsadw.I1x, ofStuff1->nltvcsadw.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->nltvcsadw.I1x, ofStuff2->nltvcsadw.I1y,
                              w, h);

            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

            delete [] alb;
            delete [] blb;

            std::printf("Exitting NLTV_CSAD with weights\n");

        }
            break;
        case M_NLTVL1_W:        // NLTV-L1 with weights
        {
            float *weight1 = ofStuff1->nltvl1w.weight;
            float *weight2 = ofStuff2->nltvl1w.weight;
            gaussian1Dweight(weight1, ofCore1->params.w_radio);
            gaussian1Dweight(weight2, ofCore2->params.w_radio);
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];

            auto *alb = new float[w*h*pd];
            auto *blb = new float[w*h*pd];

            int n_d = NL_DUAL_VAR;
            int radius = NL_BETA;

            rgb_to_lab(i0, w*h, alb);
            rgb_to_lab(i0, w*h, blb);
            // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
            nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                    ofStuff1->nltvl1w.p, ofStuff1->nltvl1w.q);
            nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                    ofStuff2->nltvl1w.p, ofStuff2->nltvl1w.q);
            if (pd!=1)
            {
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->nltvl1w.I1x, ofStuff1->nltvl1w.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->nltvl1w.I1x, ofStuff2->nltvl1w.I1y,
                              w, h);

            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

            delete [] alb;
            delete [] blb;
        }
            break;
        case M_TVCSAD_W:        // TVCSAD with weights
        {
            float *weight1 = ofStuff1->tvcsadw.weight;
            float *weight2 = ofStuff2->tvcsadw.weight;
            gaussian1Dweight(weight1, ofCore1->params.w_radio);
            gaussian1Dweight(weight2, ofCore2->params.w_radio);
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];
            int rdt = DT_R;
            int ndt = DT_NEI;
            std::printf("1 - Initializing CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->tvcsadw.pnei);
            std::printf("2 - Initializing CSAD\n");
            csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->tvcsadw.pnei);
            if (pd!=1)
            {
                // std::printf("Number of channels:%d\n",pd);
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->tvcsadw.I1x, ofStuff1->tvcsadw.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->tvcsadw.I1x, ofStuff2->tvcsadw.I1y,
                              w, h);
            *out_i0 = a_tmp;
            *out_i1 = b_tmp;
            std::printf("Exitting CSAD\n");
        }
            break;
        case M_TVL1_OCC:        // TVL1 with occlusions
        {
            auto *i1_tmp = new float[w * h];
            auto *i2_tmp = new float[w * h];
            auto *i0_tmp = new float[w * h];
            auto *i_1_tmp = new float[w * h];

            // Check if image is gray, otherwise transform
            if (pd != 1) {
                rgb_to_gray(i0, w, h, i0_tmp);
                rgb_to_gray(i1, w, h, i1_tmp);
                rgb_to_gray(i_1, w, h, i_1_tmp);
                rgb_to_gray(i2, w, h, i2_tmp);

            } else {
                memcpy(i0_tmp, i0, w * h * sizeof(float));
                memcpy(i1_tmp, i1, w * h * sizeof(float));
                memcpy(i_1_tmp, i_1, w * h * sizeof(float));
                memcpy(i2_tmp, i2, w * h * sizeof(float));
            }

            // Normalize the images between 0 and 255
            image_normalization_4(i0_tmp, i1_tmp, i_1_tmp, i2_tmp, i0_tmp, i1_tmp, i_1_tmp, i2_tmp, w * h);

            // Produce a smoothed version of the images
            gaussian(i1_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(i2_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(i0_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(i_1_tmp, w, h, PRESMOOTHING_SIGMA);

            //TODO: change the computation of derivatives
            // Check what Onofre meant
            centered_gradient(i1_tmp, ofStuff1->tvl2_occ.I1x, ofStuff1->tvl2_occ.I1y, w, h);
            centered_gradient(i_1_tmp, ofStuff1->tvl2_occ.I_1x, ofStuff1->tvl2_occ.I_1y, w, h);

            centered_gradient(i0_tmp, ofStuff2->tvl2_occ.I1x, ofStuff2->tvl2_occ.I1y, w, h);
            centered_gradient(i2_tmp, ofStuff2->tvl2_occ.I_1x, ofStuff2->tvl2_occ.I_1y, w, h);

            // Initialize g (weight)
            // The derivatives are taken from previous computation. Forward: I0, backward: I1
            init_weight(ofStuff1->tvl2_occ.g, ofStuff2->tvl2_occ.I1x, ofStuff2->tvl2_occ.I1y, w * h);
            init_weight(ofStuff2->tvl2_occ.g, ofStuff1->tvl2_occ.I1x, ofStuff1->tvl2_occ.I1y, w * h);

            // The following variables contain a gray and smooth version
            // of the corresponding image
            *out_i1 = i1_tmp;
            *out_i2 = i2_tmp;
            *out_i0 = i0_tmp;
            *out_i_1 = i_1_tmp;
        }
            break;
        default:        // TV-l2 coupled
        {
            auto *a_tmp = new float[w*h];
            auto *b_tmp = new float[w*h];
            if (pd!=1)
            {
                // std::printf("Number of channels:%d\n",pd);
                rgb_to_gray(i0, w, h, a_tmp);
                rgb_to_gray(i1, w, h, b_tmp);
            }
            else
            {
                memcpy(a_tmp,i0,w*h*sizeof(float));
                memcpy(b_tmp,i1,w*h*sizeof(float));
            }
            // normalize the images between 0 and 255
            image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
            gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
            gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
            centered_gradient(b_tmp, ofStuff1->tvl2.I1x, ofStuff1->tvl2.I1y,
                              w, h);
            centered_gradient(a_tmp, ofStuff2->tvl2.I1x, ofStuff2->tvl2.I1y,
                              w, h);
            *out_i0 = a_tmp;
            *out_i1 = b_tmp;

        }
    }
}


void of_estimation(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        float *ener_N,
        const float *i0,            // first frame
        const float *i1,            // second frame
        const float *i_1,
        const PatchIndexes index,   // Struct of indices
        const int w,
        const int h
) {

    //TODO: Each method should have its own set of parameters
    float lambda = PAR_DEFAULT_LAMBDA;
    float theta = PAR_DEFAULT_THETA;
    float tau = PAR_DEFAULT_TAU;
    float tol_OF = PAR_DEFAULT_TOL_D;
    const bool verbose = PAR_DEFAULT_VERBOSE;
    int warps = PAR_DEFAULT_NWARPS_LOCAL;

    switch (ofCore->params.val_method) {
        case M_NLTVL1:          // NLTV-L1
        {
            lambda = 2;
            theta = 0.3;
            tau = 0.1;
            // Estimate nltvl1
            guided_nltvl1(i0, i1, ofCore, &(ofStuff->nltvl1), ener_N, index,
                          lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_TVCSAD:          // TVCSAD
        {
            lambda = 0.85; //80/(49-2)
            theta = 0.3;
            tau = 0.1;
            // Estimate_tvcsad
            guided_tvcsad(i0, i1, ofCore, &(ofStuff->tvcsad), ener_N, index.ii, index.ij, index.ei, index.ej,
                          lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_NLTVCSAD:        // NLTV-CSAD
        {
            lambda = 0.85; //80/(49-2)
            theta = 0.3;
            tau = 0.1;
            // Estimate nltvcsad
            guided_nltvcsad(i0, i1, ofCore, &(ofStuff->nltvcsad), ener_N, index.ii, index.ij, index.ei, index.ej,
                            lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_TVL1_W:          // TV-l2 coupled with weights
        {
            const float central = ofStuff->tvl2w.weight[ofCore->params.w_radio + 1];
            lambda = lambda / (central * central);
            // Estimate tvl2coupled_w
            guided_tvl2coupled_w(i0, i1, ofCore, &(ofStuff->tvl2w), ener_N, index.ii, index.ij, index.ei, index.ej,
                                 lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_NLTVCSAD_W:      // NLTV-CSAD with weights
        {
            lambda = 0.85; //80/(49-2)
            const float central = ofStuff->nltvcsadw.weight[ofCore->params.w_radio + 1];
            lambda = lambda / (central * central);
            theta = 0.3;
            tau = 0.1;
            // Estimate nltvcsad_w
            guided_nltvcsad_w(i0, i1, ofCore, &(ofStuff->nltvcsadw), ener_N, index.ii, index.ij, index.ei, index.ej,
                              lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_NLTVL1_W:        // NLTV-L1 with weights
        {
            lambda = 2;
            lambda = 0.85; //80/(49-2)
            const float central = ofStuff->nltvl1w.weight[ofCore->params.w_radio + 1];
            lambda = lambda / (central * central);
            theta = 0.3;
            tau = 0.1;
            // Estimate nltvl1_w
            guided_nltvl1_w(i0, i1, ofCore, &(ofStuff->nltvl1w), ener_N, index.ii, index.ij, index.ei, index.ej,
                            lambda, theta, tau, tol_OF, warps, verbose, w, h);

        }
            break;
        case M_TVCSAD_W:        //TVCSAD with weights
        {
            lambda = 0.85; //80/(49-2)
            const float central = ofStuff->tvcsadw.weight[ofCore->params.w_radio + 1];
            lambda = lambda / (central * central);
            theta = 0.3;
            tau = 0.1;
            // Estimate tvcsad_w
            guided_tvcsad_w(i0, i1, ofCore, &(ofStuff->tvcsadw), ener_N, index.ii, index.ij, index.ei, index.ej,
                            lambda, theta, tau, tol_OF, warps, verbose, w, h);
        }
            break;
        case M_TVL1_OCC:        // TVL1 with occlusions
        {

            // Estimate_tvl2 with occlusions
            guided_tvl2coupled_occ(i0, i1, i_1, ofCore, &(ofStuff->tvl2_occ), ener_N, index, w, h);
        }
            break;
        default:                // TV-l2 coupled
            // Estimate tvl2coupled
            // Default params are those optimal for TVL1
            guided_tvl2coupled(i0, i1, ofCore, &(ofStuff->tvl2), ener_N, index.ii, index.ij, index.ei, index.ej,
                               lambda, theta, tau, tol_OF, warps, verbose, w, h);
    }
}


#endif //ENERGY_MODEL
