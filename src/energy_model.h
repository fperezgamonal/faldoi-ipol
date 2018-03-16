// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// All rights reserved.
#ifndef ENERGY_MODEL_H
#define ENERGY_MODEL_H

////INITIALIZATION OF AUXILIAR STUFF
///
///

BilateralFilterData *init_weights_bilateral(float* i0, int w, int h);

OpticalFlowData init_Optical_Flow_Data(const Parameters& params);
OpticalFlowData init_Optical_Flow_Data(float *saliency, const Parameters& params);

void initialize_auxiliar_stuff(
        SpecificOFStuff& ofStuff,
        OpticalFlowData& ofCore
        );

void free_auxiliar_stuff(SpecificOFStuff *ofStuff, OpticalFlowData *ofCore);

void prepare_stuff(SpecificOFStuff *ofStuff1,
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
                   float **out_i2);

void of_estimation(SpecificOFStuff *ofStuff,
                   OpticalFlowData *ofCore,
                   float *ener_N,
                   const float *i0,  //first frame
                   const float *i1, //second frame
                   const float *i_1,
                   const PatchIndexes index);


#endif //ENERGY_MODEL_H
