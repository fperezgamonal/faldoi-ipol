#! /usr/bin/env python3
"""
 Litle script for faldoy to execute the data from deep matches.

"""
import argparse
import os
import subprocess
import shlex
import math
import time     # added for profiling
from rescore_prunning import  confidence_values as confi
from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete
# Start global timer
init_faldoi = time.time()
# Set the arguments to compute the images
parser = argparse.ArgumentParser(description='Faldoy Minimization')
#parser.add_argument("i0", help="first image")
#parser.add_argument("i1", help="second image")
parser.add_argument("file_images", help = "File with images path")

method = 0
matchings = True
sparse_flow = True
local_of = True
global_of = True

print("Code blocks activation value:\n" +\
	"\tmatchings =\t" + str(matchings) + "\n" +\
	"\tsparse_flow =\t" + str(sparse_flow) + "\n" +\
	"\tlocal_of =\t" + str(local_of) + "\n" +\
	"\tglobal_of =\t" + str(global_of) + "\n")
#Energy model
parser.add_argument("-vm", default=str(method),
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
parser.add_argument("-wr", default='5',
                    help="Windows Radio Local patch"
                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
# Global Mininization
parser.add_argument("-warps",default='5',
                    help="Number of warps finest scale")
# Threshold for Deep Flow
parser.add_argument("-th",default='0.45',
                    help="Threshold to discard outliers from DeepFlow")
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
windows_radio = args.wr;
threshold = args.th;
r_path = args.res_path;

# C++ program names
match_comparison = "../build/deepmatching"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"


# Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())
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
match_name_1 = '%s%s_exp_mt_1.txt'%(f_path, core_name1)
sparse_name_1 = '%s%s_exp_mt_1.flo'%(f_path, core_name1)

match_name_2 = '%s%s_exp_mt_2.txt'%(f_path, core_name2)
sparse_name_2 = '%s%s_exp_mt_2.flo'%(f_path, core_name2)

region_growing = '%s%s_rg.flo'%(f_path, core_name1)
sim_value = '%s%s_exp_sim.tiff'%(f_path, core_name1)
var_flow = '%s%s_exp_var.flo'%(f_path, core_name1)

# Elapsed time (loadings)
load_timer = time.time()
print "Loading everything took " + str(load_timer - init_faldoi) + " secs."

# Obtain the matches' list for both (I0-I1 and I1-I0)
if matchings:
	max_scale = math.sqrt(2)
	# Added -nt 0 to avoid process hanging issues
	param = '%s %s -nt 0 -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name0, im_name1,max_scale,match_name_1)
	command_line = '%s %s\n'%(match_comparison, param)
	os.system(command_line)
	#I1-I0
	param = '%s %s -nt 0 -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name1, im_name0,max_scale,match_name_2)
	command_line = '%s %s\n'%(match_comparison, param)
	os.system(command_line)
	# Elapsed time (deep matches)
	matches_timer = time.time()
	print "Computing matches btw. I0 and I1 ('./deepmatching') took " + str(matches_timer - load_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	matches_timer = time.time()

#Create a sparse flow from the deep matches.
	if sparse_flow:
	param = '%s %s %s %s\n'%(cut(delete(confi(im_name0,im_name1,match_name_1,f_path),threshold)),width_im, height_im, sparse_name_1)
	command_line = '%s %s\n'%(sparse_flow, param)
	os.system(command_line)
	param = '%s %s %s %s\n'%(cut(delete(confi(im_name1,im_name0,match_name_2,f_path),threshold)),width_im, height_im, sparse_name_2)
	command_line = '%s %s\n'%(sparse_flow, param)
	os.system(command_line)
	# Elapsed time (sparse flow from matches)
	sparse_timer = time.time()
	print "Computing sparse flow from matches ('./sparse_flow') took " + str(sparse_timer - matches_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	sparse_timer = time.time()

if local_of:
	#Create a dense flow from a sparse set of initial seeds
	options = '-m %s -wr %s'%(var_m, windows_radio)
	param = '%s %s %s %s %s %s %s\n'%(im_name0, im_name1, sparse_name_1,sparse_name_2, 
		                    region_growing, sim_value, options)
	#print param
	command_line = '%s %s\n'%(match_propagation, param)
	os.system(command_line)
	# Elapsed time (dense flow)
	dense_timer = time.time()
	print "Computing dense flow from set of initial seeds ('./local_faldoi') took " + str(dense_timer - sparse_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	dense_timer = time.time()

if global_of:
	# Put the dense flow as input for a variational method
	# Tv-l2 coupled 0 Du 1
	options = '-m %s -w %s'%(var_m, warps)
	param = '%s %s %s %s %s\n'%(im_name0, im_name2,
		                    region_growing, var_flow, options)
	command_line = '%s %s\n'%(of_var, param)
	os.system(command_line)
	# Elapsed time (put the dense flow as input for a variational method)
	dense_variational_timer = time.time()
	print "Putting dense flow as input for a variational method ('./global_faldoi') took " + str(dense_variational_timer - dense_timer) + " secs."

else:
	# Need the timer anyway to compute the rest of relative values!
	dense_variational_timer = time.time()

# Global timer (ending)
end_faldoi = time.time()
print "Everything computed in " + str(end_faldoi - init_faldoi) + " secs."
