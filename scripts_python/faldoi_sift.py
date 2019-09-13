#! /usr/bin/env python3
"""
 Litle script for faldoi to execute the data from sift matches.

"""
import argparse
import os
import shlex
import subprocess
import sys
import time  # added for 'profiling'
import multiprocessing

from auxiliar_faldoi_functions import cut_matching_list as cut

# Start global timer
init_sift = time.time()
# Set the arguments to compute the images
parser = argparse.ArgumentParser(description='Faldoi Minimization')
parser.add_argument("file_images", help="File with images paths")

# Default values
#   SIFT
descriptors = True
matchings = True
def_num_scales_octave = 15

#   Sparse flow
sparse_flow_val = True

#   Local minimisation
local_of = True
def_method = 0
def_energy_params = ""
def_winsize = 5
def_local_iter = 3
def_patch_iter = 4
def_split_img = 0
def_hor_parts = 3
def_ver_parts = 2
def_fb_thresh = 2
def_partial_results = 0
partial_location = '../Results/Partial_results/'

#   Global minimisation
global_of = True
def_global_iter = 400
def_global_warps = 5

print('''Code blocks activation value:
        descriptors =   {}
        matchings =     {}
        sparse_flow =   {}
        local_of =      {}
        global_of =     {}
'''.format(descriptors, matchings, sparse_flow_val, local_of, global_of))

# Energy model
parser.add_argument("-vm", default=str(def_method),
                    help="Variational Method "
                         "(tv-l2 coupled: 0, ||Du+Du'||: 1, NLTVL1:3")
# M_TVL1       0
# M_TVL1_W     1
# M_NLTVL1     2 
# M_NLTVL1_W   3 
# M_TVCSAD     4
# M_TVCSAD_W   5
# M_NLTVCSAD   6
# M_NLTVCSAD_W 7
# M_TVL1_OCC   8       

# Energy parameters file (.txt)
parser.add_argument("-energy_params", default=str(def_method),
                    help="Optical Flow energy parameters (txt file)."
                         "The format is one parameter value per line, as follows: "
                         "lambda_value\ntheta_value\ntau_value\nbeta_value\nalpha_value\ntau_u_value\n"
                         "tau_eta_value\ntau_chi_value")

# Local Wise Minimization
# 	Window's radius
parser.add_argument("-wr", default=str(def_winsize),
                    help="Windows Radio Local patch"
                         "1 -  3x3, 2 - 5x5,...")  # (2*r +1) x (2*r+1)
#       Number of local faldoi iterations
parser.add_argument("-local_iter", default=str(def_local_iter),
                    help="Number of iterations of the local minimisation (def.=3)")

#       Number of iterations per patch (for each local iteration)
parser.add_argument("-patch_iter", default=str(def_patch_iter),
                    help="Number of iterations per patch (in each local minimisation iteration) (def.=4)")

# 	Whether to split the image into partitions or not
parser.add_argument("-split_img", default=str(def_split_img),
                    help="Enable local minimization w. subpartions instead of whole image"
                         "1 - enabled, otherwise - disabled.")

# 	Number of horizontal splits
parser.add_argument("-h_parts", default=str(def_hor_parts),
                    help="Number of horizontal parts"
                         "An integer (>0). Default is 3")

#	Number of vertical splits
parser.add_argument("-v_parts", default=str(def_ver_parts),
                    help="Number of vertical parts"
                         "An integer (>0). Default is 2")

#	FB consistency check threshold (epsilon)
parser.add_argument("-fb_thresh", default=str(def_fb_thresh),
                    help="Threshold for FB consistency check (greater ==> more permissive)"
                         "A real number (>0). Default is 2")

#	Whether to save partial results (aside from last local iteration and final flow)
#   This is usually used for debugging purposes or to show in detail the evolution of the flow field across iterations.
parser.add_argument("-partial_res", default=str(def_partial_results),
                    help="Whether to save intermediate iteration results or not"
                         "0(false) or 1(true). Default is 0")

# Global Mininization
parser.add_argument("-warps", default=str(def_global_warps),
                    help="Number of warps finest scale")

#       Number of global faldoi iterations
parser.add_argument("-glob_iter", default=str(def_global_iter),
                    help="Number of iterations of the global minimisation (def.=400)")

# Initial seeds (SIFT parameters)
parser.add_argument("-nsp", default=str(def_num_scales_octave),
                    help="Increase the sift matches")

# Gaussian weight over the data term (not used anymore ==> legacy)
parser.add_argument("-m", 		default='0',
                    help="It uses the Gaussian weight over the Data Term")

# Results "sub"path (e.g.: /Results/experiment1/iter3/)
parser.add_argument("-res_path", default='../Results/',
                    help="Subfolder under '../Results/' where data is stored")

# Verbose
parser.add_argument("-verbose", default='0', help="Toggle verbosity on/off")

args = parser.parse_args()
with open(args.file_images, 'r') as file:
    # read a list of lines into data
    data = file.readlines()
for i in range(len(data)):
    data[i] = data[i][:-1]

core_name1 = data[0][:-4].split('/')[-1]  # changed to allow for dots in the filenames
core_name2 = data[1][:-4].split('/')[-1]

var_m = args.vm
energy_params = args.energy_params
warps = args.warps
windows_radio = args.wr
loc_iter = args.local_iter
pch_iter = args.patch_iter
split_image = args.split_img
hor_parts = args.h_parts
ver_parts = args.v_parts
fb_thresh = args.fb_thresh
partial_res = args.partial_res
glb_iter = args.glob_iter
gauss = args.m
nsp = args.nsp
r_path = args.res_path
verbose = args.verbose

param_sif = "-ss_nspo {}".format(nsp)

# If the user wants to store partial results, create destination folder (if it does not exist)
if int(partial_res) == 1:
    if not os.path.exists(partial_location):
        os.makedirs(partial_location)


# Auxiliary function to parallelise calls to sift functions
def run_process(process):
    os.system("{}".format(process))


feature_descriptor = "../build/sift_cli "
match_comparison = "../build/match_cli"
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

im_name0 = os.path.expanduser(data[0])  # does nothing if no "~/folder..."
im_name1 = os.path.expanduser(data[1])

# To avoid doing the same preprocessing inside
# Get the image size
from PIL import Image
with open(im_name1, 'rb') as f:
    image = Image.open(f)
    width_im = image.size[0]
    height_im = image.size[1]

# ===============================================================================
# IF YOU DO NOT WANT/HAVE PILLOW, UNCOMMENT 3 LINES BELOW AND COMMENT 4 ABOVE)
# Using imageMagick to get width and height (PIL is not in the IPOL server)
# cmd = 'identify -ping -format "%w %h" ' + im_name1
# tmp_out = subprocess.check_output(cmd, shell=True, universal_newlines=True)
# width_im, height_im = tmp_out.split(' ')
# ===============================================================================

# os.chdir(binary_path)
desc_name_1 = "{}_sift_desc_1.txt".format(os.path.join(f_path, core_name1))
desc_name_2 = "{}_sift_desc_2.txt".format(os.path.join(f_path, core_name2))

match_name_1 = "{}_sift_mt_1.txt".format(os.path.join(f_path, core_name1))
match_name_2 = "{}_sift_mt_2.txt".format(os.path.join(f_path, core_name2))

sparse_name_1 = "{}_sift_mt_1.flo".format(os.path.join(f_path, core_name1))
sparse_name_2 = "{}_sift_mt_2.flo".format(os.path.join(f_path, core_name2))

region_growing = "{}_sift_rg.flo".format(os.path.join(f_path, core_name1))
sim_value = "{}_sift_sim.tiff".format(os.path.join(f_path, core_name1))
var_flow = "{}_sift_var.flo".format(os.path.join(f_path, core_name1))

# Elapsed time (loadings)
load_timer = time.time()
print("Loading everything took {} secs.".format(load_timer - init_sift))

# Obtain the matches' list for both (I0-I1 and I1-I0)
# Initial seeds (SIFT descriptors)
if descriptors:
    command_line_fwd = "{} {} {} > {}".format(feature_descriptor, im_name0, param_sif, desc_name_1)
    command_line_bwd = "{} {} {} > {}".format(feature_descriptor, im_name1, param_sif, desc_name_2)
    
    # Define processes to be run in parallel
    commands = (command_line_fwd, command_line_bwd)
    # Create pool of processes to be executed and map them to a thread
    pool = multiprocessing.Pool(processes=2)
    pool.map(run_process, commands)

    # Elapsed time (computing SIFT descriptors)
    desc_timer = time.time()
    print("Computing the SIFT descriptors both I0 & I1 ('./sift_cli') took {} secs.".format(desc_timer - load_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    desc_timer = time.time()  # it will add nothing to the previous clock (only time to check 'if')

# Obtain the matches' list
if matchings:
    command_line_fwd = "{} {} {} > {}".format(match_comparison, desc_name_1, desc_name_2, match_name_1)
    command_line_bwd = "{} {} {} > {}".format(match_comparison, desc_name_2, desc_name_1, match_name_2)

    # Define processes to be run in parallel
    commands = (command_line_fwd, command_line_bwd)
    # Create pool of processes to be executed and map them to a thread
    pool = multiprocessing.Pool(processes=2)
    pool.map(run_process, commands)

    # Elapsed time (matches)
    matches_timer = time.time()
    print("Computing matches btw. I0 & I1 ('./match_cli') took {}".format(matches_timer - desc_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    matches_timer = time.time()

# Create a sparse flow from the sift matches.
if sparse_flow_val:
    param = "{} {} {} {}".format(cut(match_name_1), width_im, height_im, sparse_name_1)
    
# If we use custom matches that are already filtered, cut needs to be avoided
    # param = "{} {} {} {}".format(match_name_1, width_im, height_im, sparse_name_1)
    command_line = "{} {}".format(sparse_flow, param)
    os.system(command_line)
    # Create a sparse flow from the sift matches (img I1).
    param = "{} {} {} {}".format(cut(match_name_2), width_im, height_im, sparse_name_2)
    # param = "{} {} {} {}".format(match_name_2, width_im, height_im, sparse_name_2)
    command_line = "{} {}".format(sparse_flow, param)
    os.system(command_line)
    # Elapsed time (create sparse flow from SIFT matches)
    sparse_timer = time.time()
    print("Computing sparse flow from SIFT matches ('./sparse_flow') took {}".format(sparse_timer - matches_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    sparse_timer = time.time()

# Create a dense flow from a sparse set of initial seeds
if local_of:

    options = "-m {} -p {} -wr {} -loc_it {} -max_pch_it {} -split_img {} -h_parts {} -v_parts {} -fb_thresh {}" \
              " -partial_res {} -verbose {}".format(var_m, energy_params, windows_radio, loc_iter, pch_iter,
                                                    split_image, hor_parts, ver_parts, fb_thresh, partial_res, verbose)
    param = "{} {} {} {} {} {}\n".format(args.file_images, sparse_name_1, sparse_name_2, region_growing, sim_value,
                                         options)
    command_line = "{} {}\n".format(match_propagation, param)
    os.system(command_line)
    # Elapsed time (dense flow from sparse set of initial seeds)
    dense_timer = time.time()
    print("Computing dense flow from a sparse set of initial seeds ('./local_faldoi') took {}".format(
        dense_timer - sparse_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    dense_timer = time.time()

# Put the dense flow as input for a variational method
# Tv-l2 coupled 0 Du 1
if global_of:
    options = "-m {} -p {} -w {} -glb_iters {} -verbose {}".format(var_m, energy_params, warps, glb_iter, verbose)
    param = "{} {} {} {}\n".format(args.file_images, region_growing, var_flow, options)
    command_line = "{} {}\n".format(of_var, param)
    os.system(command_line)
    # Elapsed time (Put the dense flow as input for a variational method)
    denseInputVM_timer = time.time()
    print("Putting dense flow as input for a variational method ('./global_faldoi') took {}".format(
        denseInputVM_timer - dense_timer))

else:
    # Need the timer anyway to compute the rest of relative values!
    denseInputVM_timer = time.time()

# Elapsed time (whole script)
end_sift = time.time()
print("Computing everything took {} secs.".format(end_sift - init_sift))
