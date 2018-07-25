#! /usr/bin/env python3
"""
Litle script for faldoi to execute the data from deep matches w.occlusions, if(method >= 8).

"""

import argparse
import os
import subprocess
import shlex
import multiprocessing
import math
import time	# for profiling
from rescore_prunning import  confidence_values as confi
from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete

# Init timer
init_faldoi_occ = time.time()

# Set the arguments to compute the images
parser = argparse.ArgumentParser(description = 'Faldoi Minimization')
parser.add_argument("file_images", help = "File with images")

# Default values
# 	Deep Matching
matchings = True
def_downscale = 2
def_max_scale = math.sqrt(2)
def_num_threads = 4
def_rot_plus = 45
def_rot_minus = 45
def_threshold = 0.45

#	Sparse flow
sparse_flow = True

#	Local minimisation
local_of = True
def_method = 8
def_winsize = 5

#	Global minimisation
global_of = True
def_global_warps = 5

print("Code blocks activation value:\n" +\
	"\tmatchings =\t" + str(matchings) + "\n" +\
	"\tsparse_flow =\t" + str(sparse_flow) + "\n" +\
	"\tlocal_of =\t" + str(local_of) + "\n" +\
	"\tglobal_of =\t" + str(global_of) + "\n")

# Energy model
parser.add_argument("-vm", default = str(def_method),
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
#TODO: add occlusions to the rest of functionals listed above       


# Local Wise Minimization 
parser.add_argument("-wr", default=str(def_winsize),
                    help="Windows Radio Local patch"
                         "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
# Global Minimisation
parser.add_argument("-warps", default='5',
                    help="Number of warps finest scale")

# Threshold for deepmatching's outliers
parser.add_argument("-th", default=str(def_threshold),
                    help="Threshold to discard outliers from deepmatching")
# Number of threads for Deep Matching
parser.add_argument("-nt", default=str(def_num_threads),
                    help="Number of threads to use for Deep Matching"
                         "Recommended values: up to 8 threads (improvement with more cpus is almost negligible)")

# Subsample factor
parser.add_argument("-downscale", default=str(def_downscale),
                    help="Subsample factor to reduce images dimensions (def=2)")

# Maximum scale
parser.add_argument("-max_scale", default=str(def_max_scale),
                    help="Maximum scale of deep matching (def=sqrt(2)")

# Rotation angle (positive) 
parser.add_argument("-rot_plus", default=str(def_rot_plus),
                    help="Maximum rotation angle (positive)")

# Rotation angle (negative) 
parser.add_argument("-rot_minus", default=str(def_rot_minus),
                    help="Minimum rotation angle (negative)")


# Results "sub"path (e.g.: /Results/experiment1/iter3/)
parser.add_argument("-res_path", default='../Results/',
                    help="Subfolder under '../Results/' where data is stored");

args = parser.parse_args()
with open(args.file_images, 'r') as file:
	# read a list of lines into data
	data = file.readlines()
for i in range(len(data)):
	data[i] = data[i][:-1]

sequence = data[0].split('.')[-2].split('/')[-2]
core_name1 = data[0].split('.')[-2].split('/')[-1]
core_name2 = data[1].split('.')[-2].split('/')[-1]

var_m = args.vm
warps = args.warps
windows_radio = args.wr
threshold = args.th
r_path = args.res_path
num_threads = args.nt
downscale = args.downscale
rot_plus = args.rot_plus
rot_minus = args.rot_minus
max_scale = args.max_scale

# Deep matching (at least in IPOL's server) is not faster when using 16 cpu's or more
# if (int(num_threads) > 16):
#    num_threads = '16'
# Otherwise, use as much cpu's as possible or the number inputted by the user

# Auxiliar function to handle multiprocessing (parallel calls to functions in native Python)
def run_process(process):
    os.system("{}".format(process))

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
# im_name1 = os.path.abspath(args.i0)
# im_name1 = os.path.abspath(args.i1)
# im_name0 = os.path.abspath(data[0])
# im_name1 = os.path.abspath(data[1])
im_name0 = os.path.expanduser(data[0])  # does nothing if no "~/folder..."
im_name1 = os.path.expanduser(data[1])


# Get the image size
from PIL import Image
with open(im_name1, 'rb') as f:
    image = Image.open(f)
    width_im = image.size[0]
    height_im = image.size[1]

#===============================================================================
# IF YOU DO NOT WANT/HAVE PILLOW, UNCOMMENT 3 LINES BELOW AND COMMENT 4 ABOVE)
# Using imageMagick to get width and height (PIL is not in the IPOL server)
# cmd = 'identify -ping -format "%w %h" ' + im_name1
# tmp_out = subprocess.check_output(cmd, shell=True, universal_newlines=True)
# width_im, height_im = tmp_out.split(' ')
#===============================================================================

os.chdir(binary_path)

match_name_1 = "{}{}_dm_mt_1.txt".format(f_path, core_name1)
sparse_name_1 = "{}{}_dm_mt_1.flo".format(f_path, core_name1)
sparse_in1 = "{}{}_dm_mt_1_saliency_out_cut.txt".format(f_path, core_name1)

match_name_2 = "{}{}_dm_mt_2.txt".format(f_path, core_name2)
sparse_name_2 = "{}{}_dm_mt_2.flo".format(f_path, core_name2)
sparse_in2 = "{}{}_dm_mt_2_saliency_out_cut.txt".format(f_path, core_name2)

region_growing = "{}{}_dm_rg.flo".format(f_path, core_name1)
sim_value = "{}{}_dm_sim.tiff".format(f_path, core_name1)
var_flow = "{}{}_dm_var.flo".format(f_path, core_name1)

occlusions_rg = "{}{}_dm_rg_occ.png".format(f_path, core_name1)
occlusions_var = "{}{}_dm_var_occ.png".format(f_path, core_name1)

# Elapsed time (loadings)
load_timer = time.time()
print("Loading everything took {} secs.".format(load_timer - init_faldoi_occ))

# Obtain the matches' list for both (I0-I1 and I1-I0)
if matchings:
    if int(num_threads) <= 3:
        # Execute dm calls sequentially
        nt_bwd = str(1)
        nt_fwd = str(1)
        # I0-I1
        param_fwd = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name0, im_name1, nt_fwd, downscale, max_scale, rot_plus, rot_minus, match_name_1)
        command_line_fwd = "{} {}\n".format(match_comparison, param_fwd)

        # I1-I0
        param_bwd = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name1, im_name0, nt_bwd, downscale, max_scale, rot_plus, rot_minus, match_name_2)
        command_line_bwd = "{} {}\n".format(match_comparison, param_bwd)

        os.system(command_line_fwd)
        os.system(command_line_bwd)
    else:
        if int(num_threads) > 18:  # more cpus do not improve the timing
            num_threads = '18'

        # Compute number of threads available for each matching call
        # print("Available n_cpus = {}".format(multiprocessing.cpu_count()))
        nt_bwd = int(int(num_threads)/2)
        nt_fwd = int(num_threads) - nt_bwd
        nt_bwd = str(nt_bwd)
        nt_fwd = str(nt_fwd)

        # Added -nt 0 to avoid process hanging issues
        param_fwd = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name0, im_name1, nt_fwd, downscale, max_scale, rot_plus, rot_minus, match_name_1)
        command_line_fwd = "{} {}\n".format(match_comparison, param_fwd)

        # I1-I0
        param_bwd = "{} {} -nt {} -downscale {} -max_scale {} -rot_range -{} +{} > {}".format(im_name1, im_name0, nt_bwd, downscale, max_scale, rot_plus, rot_minus, match_name_2)
        command_line_bwd = "{} {}\n".format(match_comparison, param_bwd)

        # Execute in parallel
        # Define processes to be run in parallel
        commands = (command_line_fwd, command_line_bwd)

        # Create pool of processes to be executed and map them to a thread
        pool = multiprocessing.Pool(processes=int(num_threads))
        pool.map(run_process, commands)

    # os.system(command_line)

    # Elapsed time (deep matches)
    matches_timer = time.time()
    print("Computing matches btw. I0 and I1 ('./deepmatching') took {} secs.".format(matches_timer - load_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    matches_timer = time.time()

# Create a sparse flow from the deep matches.
if sparse_flow:
    param_fwd = "{} {} {} {}\n".format(cut(delete(confi(im_name0, im_name1, match_name_1, f_path), threshold)), width_im, height_im, sparse_name_1)
    command_line_fwd = "{} {}\n".format(sparse_flow, param_fwd)
    # os.system(command_line)
    param_bwd = "{} {} {} {}\n".format(cut(delete(confi(im_name1, im_name0, match_name_2, f_path),threshold)), width_im, height_im, sparse_name_2)
    command_line_bwd = "{} {}\n".format(sparse_flow, param_bwd)
    # os.system(command_line)
    # Execute in parallel
    # Define processes to be run in parallel
    commands = (command_line_fwd, command_line_bwd)

    # Create pool of processes to be executed and map them to a thread
    pool = multiprocessing.Pool(processes=2)
    pool.map(run_process, commands)

    # Elapsed time (sparse flow from matches)
    sparse_timer = time.time()
    print("Computing sparse flow from matches ('./sparse_flow') took {} secs.".format(sparse_timer - matches_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    sparse_timer = time.time()

if local_of:
    # Create a dense flow from a sparse set of initial seeds
    options = "-m {} -wr {}".format(
        var_m,windows_radio)
    param = "{} {} {} {} {} {} {}\n".format(args.file_images, sparse_name_1, sparse_name_2,
                                         region_growing, sim_value, occlusions_rg, options)
    command_line = "{} {}\n".format(match_propagation, param)
    print(command_line)
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
    options = "-m {} -w {}".format(var_m, warps)
    param = "{} {} {} {} {} {}\n".format(args.file_images, region_growing, var_flow,
 					occlusions_rg, occlusions_var, options)
    command_line = "{} {}\n".format(of_var, param)
    print(command_line)
    os.system(command_line)
    # Elapsed time (put the dense flow as input for a variational method)
    dense_variational_timer = time.time()
    print("Putting dense flow as input for a variational method ('./global_faldoi') took {} secs.".format(dense_variational_timer - dense_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    dense_variational_timer = time.time()


# Global timer (ending)
end_faldoi_occ = time.time()
print("Everything computed in {} secs.".format(end_faldoi_occ - init_faldoi_occ))


