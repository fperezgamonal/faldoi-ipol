#!/bin/bash

time ./faldoi_deep.py ../../more_example_data/show_outliers/ambush4_out.txt -nt 4 -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/ambush4/dm/


time ./faldoi_deep.py ../../more_example_data/show_outliers/temple3_out.txt -nt 4 -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/temple3/dm/


time ./faldoi_deep.py ../../more_example_data/show_outliers/cave2_out.txt -nt 4 -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/cave2/dm/

time ./faldoi_sift.py ../../more_example_data/show_outliers/ambush4_out.txt -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/ambush4/sift/


time ./faldoi_sift.py ../../more_example_data/show_outliers/temple3_out.txt -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/temple3/sift/


time ./faldoi_sift.py ../../more_example_data/show_outliers/cave2_out.txt -vm 0 -wr 4 -local_iter 6 -res_path ../Results/show_outliers/cave2/sift/

