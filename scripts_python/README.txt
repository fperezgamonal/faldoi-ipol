========================= scripts_python/ =========================
This folder contains all the python scripts used by FALDOI.

======== TYPES ========
There are two groups of python scripts in the folder: auxiliary and execution scripts.
	- 'auxiliary': those that are used to complement the functionality of execution scripts and
 		       are called from them. List:
		- 'auxiliar_faldoi_functions.py': mainly related to pre/post processing of SIFT/
		       DeepMatching matches.

		- 'rescore_prunning.py': functions to prune DeepMatching matches according to their
 					 confidence value.

		- 'utils.py': 		(for now) only contains a function to list all the images in
					a target dataset (Middlebury' or 'MPI-Sintel').


	- 'execution': scripts that call the needed binary executables and pass the necessary +
 		       auxiliary parameters.
		- 'faldoi_sift.py': 	script that manages the execution of FALDOI with SIFT matches
 			            	(see 'usage' for a list of parameters).

		- 'faldoi_deep.py': 	script that manages the execution of FALDOI with DeepMatching
                                    	matches (see 'usage' for a list of parameters).

		- 'faldoi_deep_occ.py': script that manages the execution of FALDOI with DeepMatching
 					matches for functionals that estimate occlusions (for now
 					only tvl1). The difference with the former script is that
 					this manages the input and output files for occlusions
					(masks).


======== USAGE ========
All 3 execution scripts are used similarly. To execute FALDOI w. SIFT with the default parameters, just execute the following command in this directory ('scripts_python/'):

	./faldoi_sift.py ../example_data/final/sintel_one_frame_easy.txt

You can do the same with FALDOI+DeepMatching by calling './faldoi_deep.py' instead of './faldoi_sift.py'. For './faldoi_deep_occ.py' include '-vm 8' as it is the only functionals that models occlusions (as of August 2018), otherwise, the algorithm will revert to using the TVL1 functions without occlusions (number '0' for reference).

======== PARAMETERS ========
As shown above, the scripts only have one mandatory parameter: the text file defining the route to the frames to be processed. Aside from that, each script has several optional parameters which are defined below:

	1. Common (optional) parameters (those related to local and global minimization):
		-vm		energy functional (defined by a number):
				# M_TVL1       0
				# M_TVL1_W     1
				# M_NLTVL1     2 
				# M_NLTVL1_W   3 
				# M_TVCSAD     4
				# M_TVCSAD_W   5
				# M_NLTVCSAD   6
				# M_NLTVCSAD_W 7
				# M_TVL1_OCC   8
				Def. value = 0 (TVL1) 

		-wr		windows radius (patch_size = 2*wr + 1; wr=4 ==> 2*4+1==> 9x9 patch)
				Def. value = 5 (11x11 patch/window).

		-local_iter	number of iterations of the local minimization (iterated faldoi).
				Def. value = 3.

		-patch_iter	number of iterations per patch (for each 'local_iter').
				Def. value = 4.

		-split_img	splits input image into partitions to boost algorithm's speed by
				exploiting parallelization (more threads can be active at a time).
				It significantly reduces the execution time (w.SIFT and DM 1.71x 					faster(**)) and achieves nearly equal end-point error (maximum 0.1%
 				worse).
				Def. value = 1.
		
		-h_parts	number of horizontal parts into which the image will be split(*)
				Def. value = 3.

		-v_parts	number of vartical parts into which the image will be split(*)
				Def. value = 2.

		-warps		number of warpings performed during the final global minimization.
				Def. value = 5.

		-glob_iters	number of iterations of the global minimization (for each warping).
				Def. value = 400.

		-res_path	relative path (w.r.t. 'scripts_python') to store the results.
				Def.value = '../Results/'. Other: '../Results/sift/middlebury/test1/'


	2. Specific parameters (those related to the matching algorithm):
		2.1. SIFT (faldoi_sift.py)
			-nsp	number of scales per octave (see http://www.ipol.im/pub/art/2014/82/ 
				for more info).
				Def. value = 15.

		2.2. DeepMatching (faldoi_deep.py and faldoi_deep_occ.py)
			-th		threshold to discard outliers for deepmatching.
					Def. value = 0.45.
			-downscale	subsample factor to reduce original's image dimensions
					Def. value = 2.

			-max_scale	maximum scale for deep matching
					Def. value = sqrt(2).

			-rot_plus	maximum rotation angle (positive)
					Def. value = 45 (+45).

			-rot_minus	maximum rotation angle (negative)
					Def. value = 45 (-45).
	
			-nt		number of threads to use. Be careful setting this value
 					depending on your PC's architecture.
					Def. value = 4.

	(*)Note: if this value is very high (i.e.: > 4), specially with small images, the probability
 	of no seed falling in one partition gets bigger. If this happens, the algorithm loses
 	parallelization (as more threads sat idle). We recommend using the default values as they
	have been tested in a big cluster with up to 32 cores and produced the best overall results.

	(**) this means getting down to 1 minute for 8 cores vs nearly taking 2 minutes. For more
 	cpus (in the IPOL cluster), this difference gets increased to about 2x (down to 55 secs of
 	execution for sift).

More examples of usage may be found in the README.md of the root folder 'faldoi-ipol_[VERSION_NUMBER]'	

