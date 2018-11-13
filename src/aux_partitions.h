#ifndef AUX_PARTITIONS_H
#define AUX_PARTITIONS_H

#include "energy_structures.h"

/// Describes the headers for all the functions involved in managing partitions in FALDOI's local growing
/// Doxygen comments included in 'aux_partitions.cpp' with more information about inputs and outputs

// Initializes subimage partitions with the correct
void  init_subimage_partitions(
        const float *i0,
        const float *i1,
        const float *i_1,
        const float *i2,
        const float *i0n,
        const float *i1n,
        const float *i_1n,
        const float *i2n,
        BilateralFilterData* BiFilt_Go,
        BilateralFilterData* BiFilt_Ba,
        float *sal_go,
        float *sal_ba,
        int w_src,
        int h_src,
        int h_parts,
        int v_parts,
        std::vector<PartitionData*> *p_data,
        Parameters params
        );

// Updates the OpticalFlowData structure for any functional (common)
void update_of_data(
        OpticalFlowData *& ofGo,
        OpticalFlowData *& ofBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        int n_ch,
        bool img_to_part
        );

/// Functions to update each SpecificOFStuff struct (functional-dependant)

// tvl1-l2
void update_tvl2_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// tvl1-l2 with weights
void update_tvl1w_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// nltvl1
void update_nltvl1_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// nltvl1 with weights
void update_nltvl1w_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// tvcsad
void update_tvcsad_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// tvcsad with weights
void update_tvcsadw_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// nltvcsad
void update_nltvcsad_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// nltvcsad with weights
void update_nltvcsadw_stuffof(
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// tvl1-l2 with occlusions
void update_tvl2occ_stuffof(
        OpticalFlowData *& ofGo,
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int idx_img,
        int idx_par,
        bool img_to_part
        );


// Calls the updater functions specific to the functional selected by the user
void update_partitions_structures(
        OpticalFlowData *& ofGo,
        OpticalFlowData *& ofBa,
        SpecificOFStuff *& stuffGo,
        SpecificOFStuff *& stuffBa,
        PartitionData *& p_data,
        int n_ch,
        int idx_img,
        int idx_par,
        bool img_to_part
        );

// Updates the candidates queues for all partitions
void update_candidate_queues(
        pq_cand &queue_Go,
        pq_cand &queue_Ba,
        int n_partitions,
        std::vector<PartitionData*> *p_data
        );

// Copies information between image-wise variables and partitions-specific ones
// Note: it employs all the functions defined above
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
        int n_partitions,
        int w_src,
        int h_src,
        std::vector<PartitionData*> *p_data,
        bool img_to_part
        );

// Checks if any of the partitions queues is empty
bool anyEmptyQueues(
        std::vector<PartitionData*> *p_data,
        int n_partitions
        );

#endif //FALDOI_QTCREATOR_AUX_PARTITIONS_H
