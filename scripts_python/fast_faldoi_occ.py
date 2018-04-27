#! /usr/bin/env python3
"""
Litle script for faldoi to execute the data from deep matches w.occlusions, if(method >= 8).

"""

import argparse
import os
import subprocess
import shlex
import math
import time	# for profiling
from rescore_prunning import  confidence_values as confi
from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete

# Init timer
init_occ = time.time()

# Set the arguments to compute the images
parser = argparse.ArgumentParser(description = 'Faldoi Minimization')
parser.add_argument("file_images", help = "File with images")

method = 8
matchings = True
sparse_flow = True
local_of = True
global_of = True

print("Code blocks activation value:\n" +\
	"\tmatchings =\t" + str(matchings) + "\n" +\
	"\tsparse_flow =\t" + str(sparse_flow) + "\n" +\
	"\tlocal_of =\t" + str(local_of) + "\n" +\
	"\tglobal_of =\t" + str(global_of) + "\n")

# Energy model
parser.add_argument("-vm", default = str(method),
                    help = "Variational Method "
                    "(tv-l2 coupled: 0, tv-l2 coupled_w: 1, NLTVL1:2, NLTVL1:3...")
# M_TVL1       0
# M_TVL1_W     1
# M_NLTVL1     2 
# M_NLTVL1_W   3 
# M_TVCSAD     4
# M_TVCSAD_W   5
# M_NLTVCSAD   6
# M_NLTVCSAD_W 7
# M_TVL1_OCC   8       


# Local Wise Minimization 
parser.add_argument("-wr", default = '5',
                    help = "Windows Radio Local patch"
                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
# Global Mininization
parser.add_argument("-warps", default = '7',
                    help = "Number of warps finest scale")
# Threshold for Deep Flow
parser.add_argument("-th", default = '0.45',
                    help = "Threshold to discard outliers from DeepFlow")

# Results "sub"path (e.g.: /Results/experiment1/iter3/)
parser.add_argument("-res_path", default='../Results/',
		    help="Subfolder under '../Results/' where data is stored");

args = parser.parse_args()
with open(args.file_images, 'r') as file:
	# read a list of lines into data
	data = file.readlines()
for i in range(len(data)):
	data[i] = data[i][:-1]

sequence = data[0].split('.')[2].split('/')[-2]
core_name1 = data[0].split('.')[2].split('/')[-1]
core_name2 = data[1].split('.')[2].split('/')[-1]
var_m = args.vm
warps = args.warps
windows_radio = args.wr
threshold = args.th
r_path = args.res_path;

# C++ program names
match_comparison = "../build/deepmatching"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"
evaluation = "../build/evaluation"	# Not used(in example_data there is no ground truth!)
					# Evaluation computed end-point-error (EPE), etc. but is not included


# Set the main directory that contains all the stuff
root_path = os.getcwd()
binary_path = '../build/'

f_path = r_path
if not os.path.exists(f_path):
	os.makedirs(f_path)


# Set the images input.
im_name0 = os.path.abspath(data[0])
im_name1 = os.path.abspath(data[1])
# Get the image size
from PIL import Image
with open(im_name1, 'rb') as f:
    image = Image.open(f)
    width_im = image.size[0]
    height_im = image.size[1]

os.chdir(binary_path)

match_name_1 = '%s%s_exp_mt_1.txt'%(f_path, core_name1)
sparse_name_1 = '%s%s_exp_mt_1.flo'%(f_path, core_name1)
sparse_in1 = '%s%s_exp_mt_1_saliency_out_cut.txt'%(f_path, core_name1)

match_name_2 = '%s%s_exp_mt_2.txt'%(f_path, core_name2)
sparse_name_2 = '%s%s_exp_mt_2.flo'%(f_path, core_name2)
sparse_in2 = '%s%s_exp_mt_2_saliency_out_cut.txt'%(f_path, core_name2)

region_growing = '%s%s_rg.flo'%(f_path, core_name1)
sim_value = '%s%s_exp_sim.tiff'%(f_path, core_name1)
var_flow = '%s%s_exp_var.flo'%(f_path, core_name1)

occlusions_rg = '%s%s_rg_occ.png'%(f_path, core_name1)
occlusions_var = '%s%s_var_occ.png'%(f_path, core_name1)

# Elapsed time (loadings)
load_timer = time.time()
print("Loading everything took " + str(load_timer - init_occ) + " secs.")

# Obtain the matches' list for both (I0-I1 and I1-I0)
if matchings:
	print('Obtaining list of matches from DeepMatching')
	max_scale = math.sqrt(2)

	# I0-I1
	param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name0, im_name1, max_scale, match_name_1)
	command_line = '%s %s\n'%(match_comparison, param)
	print "Computing matches I0-I1..."	
	os.system(command_line)

	# I1-I0
	param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name1, im_name0, max_scale, match_name_2)
	command_line = '%s %s\n'%(match_comparison, param)
	print "Computing matches I1-I0..."
	os.system(command_line)
	# Elapsed time (matches)
	matches_timer = time.time()
	print "Computing matches btw. I0 and I1 took " + str(matches_timer - load_timer) + " secs."	
else:
	# Need the timer anyway to compute the rest of relative values!
	matches_timer = time.time()


# Filter matches (prunning)
cut(delete(confi(im_name0, im_name1, match_name_1, f_path), threshold))
cut(delete(confi(im_name1, im_name0, match_name_2, f_path), threshold))

# Filter matches timer
# Elapsed time (prunning)
prun_timer = time.time()
print "Filtering (prunning) matches took " + str(prun_timer - matches_timer) + " secs."

# Create a sparse flow from the deepmatching matches.
if sparse_flow:
	param = '%s %s %s %s\n'%(sparse_in1, width_im, height_im, sparse_name_1)
	command_line = '%s %s\n'%(sparse_flow, param)
	print('Creating sparse from matches I0-I1...')
	os.system(command_line)
	param = '%s %s %s %s\n'%(sparse_in2, width_im, height_im, sparse_name_2)
	command_line = '%s %s\n'%(sparse_flow, param)
	os.system(command_line)
	print('Creating sparse from matches I1-I0...')

	# Elapsed time (sparse flow)
	sparse_timer = time.time()
	print "Computing a sparse flow from matches took " + str(sparse_timer - prun_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	sparse_timer = time.time()

# Create a dense flow from a sparse set of initial seeds
if local_of:
	options = '-m %s -wr %s'%(var_m, windows_radio)
	param = '%s %s %s %s %s %s %s\n'%(args.file_images, sparse_name_1, sparse_name_2, 
		                    region_growing, sim_value, occlusions_rg, options)

	command_line = '%s %s\n'%(match_propagation, param)
	print('Computing local faldoi...')
	os.system(command_line)
	
	dense_timer = time.time()
	print "Computing a dense flow from a set of initial seeds (local faldoi) took " + str(dense_timer - sparse_timer) + " secs."

# Put the dense flow as input for a variational method
if global_of:
	# Tv-l2 coupled 0 Du 1
	options = '-m %s -w %s'%(var_m, warps)
	param = '%s %s %s %s %s %s\n'%(args.file_images,
		                    region_growing, var_flow, occlusions_rg, occlusions_var, options)
	command_line = '%s %s\n'%(of_var, param)

	print('Computing global faldoi')
	os.system(command_line)
	
	inVA_timer = time.time()
	print "Putting the dense flow as input for a variational method (global faldoi) took " + str(inVA_timer - dense_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	inVA_timer = time.time()

end_faldoi_occ = time.time()
print "Everything computed in a total of " + str(end_faldoi_occ - init_dm_occ) + " secs."


