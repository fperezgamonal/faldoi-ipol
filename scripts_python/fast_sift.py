#! /usr/bin/python3
"""
 Litle script for faldoy to execute the data from sift matches.

"""
import argparse
import os
import subprocess
import shlex
import time  # added for 'profiling'
from auxiliar_faldoi_functions import cut_matching_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete

# Start global timer
init_sift = time.time()
# Set the arguments to compute the images
parser = argparse.ArgumentParser(description='Faldoi Minimization')
# Need to change this to fully integrate the original faldoi's functionality and
# the TFM's into one project (currently there are some bugs that do not allow this)
#parser.add_argument("i0", help="first image")
#parser.add_argument("i1", help="second image")
parser.add_argument("file_images", help = "File with images path")

method = 0
descriptors = True
matchings = True
sparse_flow = True
local_of = True
global_of = True

print("Code blocks activation value:\n" +\
	"\tdescriptors =\t" + str(descriptors) + "\n" +\
	"\tmatchings =\t" + str(matchings) + "\n" +\
	"\tsparse_flow =\t" + str(sparse_flow) + "\n" +\
	"\tlocal_of =\t" + str(local_of) + "\n" +\
	"\tglobal_of =\t" + str(global_of) + "\n")
	

# Energy model
parser.add_argument("-vm", default=str(method),
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


# Local Wise Minimization
parser.add_argument("-wr", default='5',
                    help="Windows Radio Local patch"
                         "1 -  3x3, 2 - 5x5,...")  # (2*r +1) x (2*r+1)
# Global Mininization
parser.add_argument("-warps", default='5',
                    help="Number of warps finest scale")
# Initial seeds (SIFT parameters)
parser.add_argument("-nsp", default='15',
                    help="Increase the sift matches")

parser.add_argument("-m", default='0',
                    help="It uses the Gaussian weight over the Data Term");
# Results "sub"path (e.g.: /Results/experiment1/iter3/)
parser.add_argument("-res_path", default='../Results/',
		    help="Subfolder under '../Results/' where data is stored");

#args = parser.parse_args()
#core_name1 = args.i0.split('.')[-2].split('/')[-1]
#core_name2 = args.i1.split('.')[-2].split('/')[-1]

args = parser.parse_args()
with open(args.file_images, 'r') as file:
	# read a list of lines into data
	data = file.readlines()
for i in range(len(data)):
	data[i] = data[i][:-1]
	#tmp = data[i][:-1]
	#data[i] = os.path.expanduser(tmp)
	#print("data[i]= " + data[i])

# Save to tmp file
tmp_filename = "tmp_absPaths.txt"
#with open(tmp_filename, 'w') as tmp_file:
	# Store the modified image paths' file
	#tmp_file.writelines(data)
#	tmp_file.writelines(["%s\n" % item  for item in data])
sequence = data[0].split('.')[-2].split('/')[-2]  # not used
core_name1 = data[0].split('.')[-2].split('/')[-1]
core_name2 = data[1].split('.')[-2].split('/')[-1]

var_m = args.vm
warps = args.warps
windows_radio = args.wr;
gauss = args.m;
nsp = args.nsp;
r_path = args.res_path;

param_sif = '-ss_nspo %s' % (nsp)

feature_descriptor = "../build/sift_cli "
match_comparison = "../build/match_cli"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"


#Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())
#print(root_path)
#binary_path = root_path + "bin/"
binary_path = '../build/'
#f_path = root_path + "Results/"
f_path = r_path;
if not os.path.exists(f_path):
    os.makedirs(f_path)
# Set the folder where the binaries are.
# Set the images input.
#im_name1 = os.path.abspath(args.i0)
#im_name2 = os.path.abspath(args.i1)
# Get the image size
#from PIL import Image

#with open(im_name1, 'r') as f:
#    image = Image.open(f)
#    width_im = image.size[0]
#    height_im = image.size[1]

#Set the images input.
# Included changes to make the reading not user-dependent even if
# we use paths such as "~/Folder/subfolder/subsubfolder/.../i0.png"
# relative paths w.r.t current dir work also fine (e.g.: "../example/i0.png")
#im_name0 = os.path.abspath(data[0])
#im_name1 = os.path.abspath(data[1])
im_name0 = os.path.expanduser(data[0])  # does nothing if no "~/folder..."
im_name1 = os.path.expanduser(data[1])

# To avoid doing the same preprocessing inside
#Get the image size
from PIL import Image
with open(im_name1, 'rb') as f:
    image = Image.open(f)
    width_im = image.size[0]
    height_im = image.size[1]

#os.chdir(binary_path)
desc_name_1 = '%s%s_sift_desc_1.txt' % (f_path, core_name1)
desc_name_2 = '%s%s_sift_desc_2.txt' % (f_path, core_name2)

match_name_1 = '%s%s_sift_mt_1.txt' % (f_path, core_name1)
match_name_2 = '%s%s_sift_mt_2.txt' % (f_path, core_name2)

sparse_name_1 = '%s%s_sift_mt_1.flo' % (f_path, core_name1)
sparse_name_2 = '%s%s_sift_mt_2.flo' % (f_path, core_name2)

region_growing = '%s%s_sift_rg.flo' % (f_path, core_name1)
sim_value = '%s%s_sift_sim.tiff' % (f_path, core_name1)
var_flow = '%s%s_sift_var.flo' % (f_path, core_name1)

# Elapsed time (loadings)
load_timer = time.time()
print("Loading everything took " + str(load_timer - init_sift) + " secs.")

# Obtain the matches' list for both (I0-I1 and I1-I0)
# Initial seeds (SIFT descriptors)
if descriptors:
	command_line = '%s %s %s\n' % (feature_descriptor, im_name0, param_sif)
	exe_prog(command_line, desc_name_1)
	command_line = '%s %s %s\n' % (feature_descriptor, im_name1, param_sif)
	exe_prog(command_line, desc_name_2)
	# Elapsed time (computing SIFT descriptors)
	desc_timer = time.time()
	print("Computing the SIFT descriptors both I0 & I1 ('./sift_cli') took " + str(desc_timer - load_timer) + " secs.")

else:
	# Need the timer anyway to compute the rest of relative values!
	desc_timer = time.time()	# it will add nothing to the previous clock (only time to check 'if')

# Obtain the matches' list
if matchings:
	command_line = '%s %s %s\n' % (match_comparison, desc_name_1, desc_name_2)
	exe_prog(command_line, match_name_1)
	command_line = '%s %s %s\n' % (match_comparison, desc_name_2, desc_name_1)
	exe_prog(command_line, match_name_2)
	# Elapsed time (matches)
	matches_timer = time.time()
	print("Computing matches btw. I0 & I1 ('./match_cli') took " + str(matches_timer - desc_timer) + " secs.")

else:
	# Need the timer anyway to compute the rest of relative values!
	matches_timer = time.time()

# Create a sparse flow from the sift matches.
if sparse_flow:
	param = '%s %s %s %s\n' % (cut(match_name_1), width_im, height_im, sparse_name_1)
	command_line = '%s %s\n' % (sparse_flow, param)
	os.system(command_line)
	# Create a sparse flow from the sift matches (img I1).
	param = '%s %s %s %s\n' % (cut(match_name_2), width_im, height_im, sparse_name_2)
	command_line = '%s %s\n' % (sparse_flow, param)
	os.system(command_line)
	# Elapsed time (create sparse flow from SIFT matches)
	sparse_timer = time.time()
	print("Computing sparse flow from SIFT matches ('./sparse_flow') took " + str(sparse_timer - matches_timer) + " secs.")

else:
	# Need the timer anyway to compute the rest of relative values!
	sparse_timer = time.time()

# Load absolute paths (hacky*, maybe there is a better way!)
# * converted original path files to absolute equivalents and store it (since loc_faldoi only accepts .txt input)
#print("we are at: " + os.getcwd())
#print("List of files at current directory: ")
#print(os.listdir('.'))
#print(tmp_filename)
#print(os.path.join(os.getcwd(), tmp_filename))
#print("is file (directly file): " + str(os.path.isfile(tmp_filename)))
#print("is file (abs path): " + str(os.path.isfile(os.path.join(os.getcwd(), tmp_filename))))
#print(os.path.abspath(tmp_filename))
#print(os.path.join(os.path.abspath(__file__), tmp_filename))
#os.chdir(root_path)
#with open(r"/home/fperezgamonal/Documents/Papers_code/Faldoi_tfm-master/scripts_python/tmp_absPaths.txt", 'r') as tmp_file: # Closed at the end

# Create a dense flow from a sparse set of initial seeds
if local_of:
	options = '-m %s -wr %s' % (var_m, windows_radio)
	param = '%s %s %s %s %s %s\n' % (args.file_images, sparse_name_1, sparse_name_2,
			                    region_growing, sim_value, options)
	# print param
	command_line = '%s %s\n' % (match_propagation, param)
	#print("l_of cmd:\n{}\n".format(command_line))
	os.system(command_line)
	# Elapsed time (dense flow from sparse set of initial seeds)
	dense_timer = time.time()
	print("Computing dense flow from a sparse set of initial seeds ('./local_faldoi') took " + str(dense_timer - sparse_timer) + " secs.")

else:
	# Need the timer anyway to compute the rest of relative values!
	dense_timer = time.time()

	# Put the dense flow as input for a variational method
	# Tv-l2 coupled 0 Du 1
if global_of:
	options = '-m %s -w %s' % (var_m, warps)
	param = '%s %s %s %s\n' % (args.file_images,
			              region_growing, var_flow, options)
	command_line = '%s %s\n' % (of_var, param)
	os.system(command_line)
	print(command_line)
	# Elapsed time (Put the dense flow as input for a variational method)
	denseInputVM_timer = time.time()
	print("Putting dense flow as input for a variational method ('./global_faldoi') took " + str(denseInputVM_timer - dense_timer) + " secs.")

else:
	# Need the timer anyway to compute the rest of relative values!
	denseInputVM_timer = time.time()
# File is closed here
# Delete temporal absolute paths file
#os.remove(tmp_filename)
# Elapsed time (whole script)
end_sift = time.time()
print("Computing everything took " + str(end_sift - init_sift) + " secs.")
