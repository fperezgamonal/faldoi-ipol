#!/bin/bash

EPSILON_MIN=0.5
EPSILON_MAX=5
EPSILON_STEP=0.25

#EVAL_DIR='../../Evaluation_code/AEE_EPE_error/'

#for EPSILON in `seq $EPSILON_MIN $EPSILON_STEP $EPSILON_MAX` 
for EPSILON in {6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9,9.25,9.5,9.75,10}
#for EPSILON in {4.75}
do
#EPSILON=4.75
# Avoid overwriting (if the first directory exists, we assume the rest do (for simplicity))
	#if [ -d "../Results/robustness_epsilon/dm_${EPSILON}/easy_sift/" ]; then
	#	echo "Path already exists, exiting without overwriting"
	#else
		# DM
		time ./faldoi_deep.py ../../more_example_data/final/sintel_one_frame_easy.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/easy/

		time ./faldoi_deep.py ../../more_example_data/final/sintel_one_frame_medium.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/medium/

		time ./faldoi_deep.py ../../more_example_data/final/sintel_one_frame_hard.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/hard/

		time ./faldoi_deep.py ../../more_example_data/dm_vs_sift/easy_sift.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/easy_sift/

                time ./faldoi_deep.py ../../more_example_data/dm_vs_sift/difficult_sift.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/hard_sift/

                time ./faldoi_deep.py ../../more_example_data/iterated_ambush6/iter_ambush6.txt -nt 4 -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/dm_${EPSILON}/ambush_6/


		# SIFT
		time ./faldoi_sift.py ../../more_example_data/final/sintel_one_frame_easy.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/easy/

		time ./faldoi_sift.py ../../more_example_data/final/sintel_one_frame_medium.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/medium/

		time ./faldoi_sift.py ../../more_example_data/final/sintel_one_frame_hard.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/hard/
		time ./faldoi_sift.py ../../more_example_data/dm_vs_sift/easy_sift.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/easy_sift/

                time ./faldoi_sift.py ../../more_example_data/dm_vs_sift/difficult_sift.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/hard_sift/

                time ./faldoi_sift.py ../../more_example_data/iterated_ambush6/iter_ambush6.txt -vm 0 -wr 4 -fb_thresh ${EPSILON} -res_path ../Results/robustness_epsilon/sift_${EPSILON}/ambush_6/


	#fi
done

# EVALUATION (done in MATLAB)

#matlab -nodesktop -nosplash -r "computeAEE_EPE(path_estim, path_gt, path_occ, path_invalid, 1);

