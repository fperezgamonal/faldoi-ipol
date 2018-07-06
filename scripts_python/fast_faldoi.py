#! /usr/bin/env python3
"""
 Litle script for faldoy to execute the data from deep matches.

"""
import argparse
import math
import os
import shlex
import subprocess
import time  # added for 'profiling'

from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import delete_outliers as delete
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from rescore_prunning import confidence_values as confi

# Start global timer
init_faldoi = time.time()
# Set the arguments to compute the images
parser = argparse.ArgumentParser(description='Faldoy Minimization')
#parser.add_argument("i0", help="first image")
#parser.add_argument("i1", help="second image")
parser.add_argument("file_images", help = "File with images path")

# Default values
# 	Deep Matching
matchings = True
def_downscale = 2
def_scale = math.sqrt(2)
def_num_threads = 4
def_rot_range_plus = 45
def_rot_range_minus = 45
def_threshold = 0.45

#	Sparse flow
sparse_flow = True

#	Local minimisation
local_of = True
def_method = 0
def_winsize = 5
def_local_iter = 3
def_patch_iter = 4
def_split_img = 0
def_hor_parts = 3
def_ver_parts = 2

#	Global minimisation
global_of = True
def_global_iter = 400
def_global_warps = 5

print("Code blocks activation value:\n" + \
	  "\tmatchings =\t" + str(matchings) + "\n" + \
	  "\tsparse_flow =\t" + str(sparse_flow) + "\n" + \
	  "\tlocal_of =\t" + str(local_of) + "\n" + \
	  "\tglobal_of =\t" + str(global_of) + "\n")
#Energy model
parser.add_argument("-vm", default=str(def_method),
					help="Variational Method "
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
parser.add_argument("-wr", default=str(def_winsize),
					help="Windows Radio Local patch"
						 "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)

#       Number of local faldoi iterations
parser.add_argument("-local_iter", default=str(def_local_iter),
					help="Number of iterations of the local minimisation (def.=3)")

#       Number of iterations per patch (for each local iteration)
parser.add_argument("-patch_iter", default=str(def_patch_iter),
					help="Number of iterations per patch (in each local minimisation iteration) (def.=4)")

# 	Whether to split the image into partitions or not
parser.add_argument("-split_img", default=str(def_split_img),
					help="Enable local minimization w. subpartions instead of whole image"
						 "1 - enabled, othewise - disabled.")

# 	Number of horizontal splits
parser.add_argument("-h_parts", default=str(def_hor_parts),
					help="Number of horizontal parts"
						 "An integer (>0). Default is 3")

#	Number of vertical splits
parser.add_argument("-v_parts", default=str(def_ver_parts),
					help="Number of vertical parts"
						 "An integer (>0). Default is 2")
# Global Mininization
parser.add_argument("-warps", default='5',
					help="Number of warps finest scale")

#       Number of global faldoi iterations
parser.add_argument("-glob_iter", default=str(def_global_iter),
					help="Number of iterations of the global minimisation (def.=400)")

# Threshold for Deep Flow
parser.add_argument("-th", default=str(def_threshold),
					help="Threshold to discard outliers from DeepFlow")
# Number of threads for Deep Matching
parser.add_argument("-nt", default=str(def_num_threads),
					help="Number of threads to use for Deep Matching"
						 "Recommended values: up to 8 threads (improvement with more cpus is almost negligible)")

# Results "sub"path (e.g.: /Results/experiment1/iter3/)
parser.add_argument("-res_path", default='../Results/',
					help="Subfolder under '../Results/' where data is stored");


#args = parser.parse_args()
#core_name1 = args.i0.split('.')[0].split('/')[-1]
#core_name2 = args.i1.split('.')[0].split('/')[-1]
args = parser.parse_args()
with open(args.file_images, 'r') as file:
	# read a list of lines into data
	data = file.readlines()
for i in range(len(data)):
	data[i] = data[i][:-1]

sequence = data[0].split('.')[-2].split('/')[-2]  # not used
core_name1 = data[0].split('.')[-2].split('/')[-1]
core_name2 = data[1].split('.')[-2].split('/')[-1]

var_m = args.vm
warps = args.warps
windows_radio = args.wr
loc_iter = args.local_iter
pch_iter = args.patch_iter
split_image = args.split_img
hor_parts = args.h_parts
ver_parts = args.v_parts
glb_iter = args.glob_iter
threshold = args.th
r_path = args.res_path
num_threads = args.nt

# Deep matching (at least in IPOL's server) is not faster when using 16 cpu's or more
if (int(num_threads) > 16):
    num_threads = '16'
    # Otherwise, use as much cpu's as possible or the number inputted by the user

# C++ program names
match_comparison = "../build/deepmatching"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"


# Set the main directory that contains all the stuff
root_path = "{}/".format(os.getcwd())
binary_path = '../build/'
f_path = r_path
if not os.path.exists(f_path):
	os.makedirs(f_path)
# Set the folder where the binaries are.
# Set the images input.
#im_name1 = os.path.abspath(args.i0)
#im_name1 = os.path.abspath(args.i1)
#im_name0 = os.path.abspath(data[0])
#im_name1 = os.path.abspath(data[1])
im_name0 = os.path.expanduser(data[0])  # does nothing if no "~/folder..."
im_name1 = os.path.expanduser(data[1])


# Get the image size
from PIL import Image
with open(im_name1, 'rb') as f:
	image = Image.open(f)
	width_im = image.size[0]
	height_im = image.size[1]

os.chdir(binary_path)
match_name_1 = "{}{}_dm_mt_1.txt".format(f_path, core_name1)
sparse_name_1 = "{}{}_dm_mt_1.flo".format(f_path, core_name1)

match_name_2 = "{}{}_dm_mt_2.txt".format(f_path, core_name2)
sparse_name_2 = "{}{}_dm_mt_2.flo".format(f_path, core_name2)

region_growing = "{}{}_dm_rg.flo".format(f_path, core_name1)
sim_value = "{}{}_dm_sim.tiff".format(f_path, core_name1)
var_flow = "{}{}_dm_var.flo".format(f_path, core_name1)

# Elapsed time (loadings)
load_timer = time.time()
print("Loading everything took {} secs.".format(load_timer - init_faldoi))

# Obtain the matches' list for both (I0-I1 and I1-I0)
if matchings:
	# Added -nt 0 to avoid process hanging issues
	param = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name0, im_name1, num_threads, def_downscale, def_scale, def_rot_range_plus, def_rot_range_minus, match_name_1)
	command_line = "{} {}\n".format(match_comparison, param)
	os.system(command_line)
	#I1-I0
	param = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name1, im_name0, num_threads, def_downscale, def_scale, def_rot_range_plus, def_rot_range_minus, match_name_2)
	command_line = "{} {}\n".format(match_comparison, param)
	os.system(command_line)
	# Elapsed time (deep matches)
	matches_timer = time.time()
	print("Computing matches btw. I0 and I1 ('./deepmatching') took {} secs.".format(matches_timer - load_timer))

else:
	# Need the timer anyway to compute the rest of relative values!
	matches_timer = time.time()

#Create a sparse flow from the deep matches.
if sparse_flow:
	param = "{} {} {} {}\n".format(cut(delete(confi(im_name0, im_name1, match_name_1, f_path), threshold)), width_im, height_im, sparse_name_1)
	command_line = "{} {}\n".format(sparse_flow, param)
	os.system(command_line)
	param = "{} {} {} {}\n".format(cut(delete(confi(im_name1, im_name0, match_name_2, f_path),threshold)), width_im, height_im, sparse_name_2)
	command_line = "{} {}\n".format(sparse_flow, param)
	os.system(command_line)
	# Elapsed time (sparse flow from matches)
	sparse_timer = time.time()
	print("Computing sparse flow from matches ('./sparse_flow') took {} secs.".format(sparse_timer - matches_timer))

else:
	# Need the timer anyway to compute the rest of relative values!
	sparse_timer = time.time()

if local_of:
	#Create a dense flow from a sparse set of initial seeds
	options = "-m {} -wr {} -loc_it {} -max_pch_it {} -split_img {} -h_parts {} -v_parts {}".format(
		var_m,windows_radio, loc_iter, pch_iter, split_image, hor_parts, ver_parts)
	param = "{} {} {} {} {} {}\n".format(args.file_images, sparse_name_1, sparse_name_2,
										 region_growing, sim_value, options)
	command_line = "{} {}\n".format(match_propagation, param)
	os.system(command_line)
	# Elapsed time (dense flow)
	dense_timer = time.time()
	print("Computing dense flow from set of initial seeds ('./local_faldoi') took {} secs.".format(
		dense_timer - sparse_timer))

else:
	# Need the timer anyway to compute the rest of relative values!
	dense_timer = time.time()

if global_of:
	# Put the dense flow as input for a variational method
	# Tv-l2 coupled 0 Du 1
	options = "-m {} -w {} -glb_iters {}".format(var_m, warps, glb_iter)
	param = "{} {} {} {}\n".format(args.file_images,
								   region_growing, var_flow, options)
	command_line = "{} {}\n".format(of_var, param)
	os.system(command_line)
	# Elapsed time (put the dense flow as input for a variational method)
	dense_variational_timer = time.time()
	print("Putting dense flow as input for a variational method ('./global_faldoi') took {} secs.".format(
		dense_variational_timer - dense_timer))

else:
	# Need the timer anyway to compute the rest of relative values!
	dense_variational_timer = time.time()

# Global timer (ending)
end_faldoi = time.time()
print("Everything computed in {} secs.".format(end_faldoi - init_faldoi))
