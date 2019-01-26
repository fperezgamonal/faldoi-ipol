[![Python 2.7](https://img.shields.io/badge/python-2.7-green.svg)](https://www.python.org/downloads/release/python-270/)
[![Python 3.5](https://img.shields.io/badge/python-3.5-green.svg)](https://www.python.org/downloads/release/python-350/)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# FALDOI-IPOL
Stems from the basic [FALDOI: A New Minimization Strategy for Large Displacement Variational Optical Flow][1] by Roberto P. Palomares, Enric Meinhardt-Lopis, Coloma Ballester and Gloria Haro algorithm and aims to add occlusion estimation to several energy functionals and optimise the code to be published on the [IPOL journal](http://www.ipol.im/) with an interactive demo.

## Paper(s) and citation
If you use FALDOI, please cite *any* of the following papers:

[Springer original paper (November 2016)][1]:

```
@Article{Palomares2017,
author="Palomares, Roberto P.
and Meinhardt-Llopis, Enric
and Ballester, Coloma
and Haro, Gloria",
title="FALDOI: A New Minimization Strategy for Large Displacement Variational Optical Flow",
journal="Journal of Mathematical Imaging and Vision",
year="2017",
month="May",
day="01",
volume="58",
number="1",
pages="27--46",
issn="1573-7683",
doi="10.1007/s10851-016-0688-y",
url="https://doi.org/10.1007/s10851-016-0688-y"
}
```

the [arXiv paper (September 2016)][2]:

```
@article{DBLP:journals/corr/PalomaresHBM16,
  author    = {Roberto P. Palomares and
               Gloria Haro and
               Coloma Ballester and
               Enric Meinhardt{-}Llopis},
  title     = {{FALDOI:} Large Displacement Optical Flow by Astute Initialization},
  journal   = {CoRR},
  volume    = {abs/1602.08960},
  year      = {2016},
  url       = {http://arxiv.org/abs/1602.08960},
  archivePrefix = {arXiv},
  eprint    = {1602.08960},
  timestamp = {Wed, 07 Jun 2017 14:42:11 +0200},
  biburl    = {https://dblp.org/rec/bib/journals/corr/PalomaresHBM16},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

You can reference this implementation with: [IPOL article, demo and supplementary material (December 2018)][3]:
```
bibTeX for IPOL
``` 
___
## Table of contents
* [Getting Started](#getting-started)
  * [Pre-requisites](#pre-requisites)
  * [Compilation](#compilation)
  * [Algorithm's steps](#algorithm's-steps)
  * [From algorithms to actual code](#from-algorithms-to-actual-code)
* [Execution](#execution)
  * [C++ executables](#c-executables)
  * [Python scripts](#python-scripts)
    * [Example of usage](#example-of-usage)
* [Bugs, issues and questions](#bugs,-issues-and-questions)
* [Developers](#developers)
* [License and copyright](#license-and-copyright)
___

## Getting Started
The code is written in C/C++ and includes Python scripts to unify and simplify all the tasks the algorithms performs.

### Pre-requisites
The software needs the following programs to be installed in order to function:

- Python 3.5 (also works with Python 2.7)
- Pillow for Python 3 (if you do not want to use Pillow, in the python scripts there is a commented section above the PIL import that uses [imagemagick's identify](https://www.imagemagick.org/script/identify.php) instead). Nevertheless, notice that Pillow will be usually faster (since it is built-in python).
- lipng, libjpeg, libtiff: by default, we use your system's default version (although we provide the original version of libpng which the code was created with). Default system versions of Ubuntu 16.04 LTS have been checked as well.

- OpenMP (Optional but recommended). Should be included with your compiler if you have a relatively new version (e.g.: gcc supports it since version 4.2). For a gentle introduction to OpenMP, please read [this post](https://helloacm.com/simple-tutorial-with-openmp-how-to-use-parallel-block-in-cc-using-openmp/) or visit the official website.

### Compilation
The code needs to be compiled in order to generate the needed executables. Nevertheless, if you do not want to compile the code before testing it once, we added already-compiled executables that are already linked (see [_Execution_](#Execution)). These executables have been compiled in a machine with 4 cores running Ubuntu 16.04 LTS.

To compile the code, we have added a bash script that runs all the needed commands to remove the old _build_ directory, create a new one and re-compile the whole project, integrating all the changes done to any C/C++ file. To run the script, just do the following (from the project's root directory):
```bash
source recompile_changes.sh
```
or
```bash
. recompile_changes.sh
```

'source' or '.' is a bash shell built-in command that executes the content of the file passed as argument, in the current shell.

If you run into any problems with the binaries for the matchers located in 'ext_bin', you may try to download their sources and compile it under your machine. To do so, you can follow the steps explained in the ['Compute matches' section](#compute-matches) below.

### Algorithm's steps
In order to obtain the final optical flow estimation, the algorithm follows these steps:

1. Reads a text file (.txt) containing the paths to the input frames to be processed (see [Python scripts/NOTE](#python-scripts) for details).
2. Computes the descriptors and matches between frames I0 and I1 (with SIFT or Deepmatching).
3. Computes a sparse flow from the matches.
4. Computes a dense flow from the initial sparse flow (local step).
5. Computes global minimization taking as input the flow of the previous step (global step).

### From algorithms to actual code
The first time you download the code, you will see that the *'src'* folder contains several source files. Trying to understand what each of those files does at once will be too time-consuming and hard. For those reasons, we **strongly** encourage you to read one of the articles linked at the top of this file. More concretely, if you want a thorough explanation of the theoretical background of the algorithm, read the [original Springer paper][1]; if you want a summarised version containing more implementations details, read the IPOL article instead.

In the IPOL article, we have broken down the algorithm into 8 different pseudo-codes which are implemented in code inside the *'src'* folder. Note that we will refer to the algorithms with the same notation used in the IPOL article. The algorithms are the following:

- **Algorithm 1 (main algorithm (end-to-end))**: it is implemented in the python scripts ['faldoi_sift.py'](scripts_python/faldoi_sift.py) and ['faldoi_deep.py'](scripts_python/faldoi_deep.py), depending on the chosen matching method.
- **Algorithm 2 (compute-matches)**: it is implemented in lines [#201 to #240](scripts_python/faldoi_sift.py#L201) of *'faldoi_sift.py'* and lines [#223 to #273](scripts_python/faldoi_deep.py#L223) of *'faldoi_deep.py'*.
- **Algorithm 3 (saliency-matches)**: this algorithm is only run if the selected matcher is DeepMatching. It is implemented in the file ['rescore_prunning.py'](scripts_python/rescore_prunning.py) (compute confidence scores) and *'auxiliar_faldoi_functions.py'*, lines [#47 to #64](scripts_python/auxiliar_faldoi_functions.py) (delete outliers) and lines [#32 to #44](scripts_python/auxiliar_faldoi_functions.py#L32) (cut matches list, removing scores).
- **Algorithm 4 (sparse-flow)**: it is implemented in the file ['sparse_flow.cpp'](src/sparse_flow.cpp).
- **Algorithm 5 (dense-flow)**: it is implemented in the file *'local_faldoi.cpp'*, in the function *'match_growing_variational'* (lines [#1174 to #1530](src/local_faldoi.cpp#L1174)).
- **Algorithm 6 (basic-faldoi-growing)**: it is implemented in the file *'local_faldoi.cpp'*, in the function *'local_growing'* (lines [#890 to #1039](src/local_faldoi.cpp#L890)).
- **Algorithm 7 (forwardBackward-pruning)**: it is implemented in the file *'local_faldoi.cpp'*, in the function *'fb_consistency_check'* (lines [#166 to #189](src/local_faldoi.cpp#L166)).
- **Algorithm 8 (flow-refinement)**: this algorithm has two different implementations following the pseudo-code written in the article (one for the local step, over a patch; one for the global step, over the whole image domain). <br> <br> 
As a consequence, it is implemented in ['global_faldoi.cpp'](src/global_faldoi.cpp), inside the *main()* function, using different implementations for each type of functional that FALDOI supports. For example, the *TVL2-L1* functional, the global step is performed in function *'tvl2OF'* (lines [#555 to #843](src/global_faldoi.cpp#L555)). <br>
For the local step, the implementation resides in functional-specific functions defined alongside the energy model in its own source + header (cpp + h) files. For instance, for the *TVL2-L1* functional, the pseudo-code is implemented in the file *'tvl2_model.cpp'* (lines [#249 to #425](src/local_faldoi.cpp#L249)).

For more information about all these algorithms, please, refer to the papers (section *'Proposed minimization strategy'*).
Additionally, we provide a list of the relation between the naming convention for the energy functionals used in the IPOL paper and in the code.
- Paper: "TVL2-L1", code: "TVL1". Implemented in ['tvl1_model.cpp'](src/tvl1_model.cpp)
- Paper: "NLTV-CSAD", code "NLTVCSAD". Implemented in ['nltvcsad_model.cpp'](src/nltvcsad_model.cpp)
- Paper: "TVL2-CSAD", code "TVCSAD". Implemented in ['tvcsad_model.cpp'](src/tvcsad_model.cpp)

Not in the paper, but combination of different data term (CSAD, L2) and regularization term (L1, NLTV):
- "NLTV1" combines NLTV as regularizer + L1 as data term. Implemented in ['nltv_model.cpp'](src/nltv_model.cpp)

## Execution
Once you have successfully compiled the program, you can execute the algorithm in two ways: using the C++ executables or running the Python scripts (**recommended**). The advantage of choosing the latter is that all the algorithm's steps are run at once sequentially, while with C++ executables, you will have to manually call each block of the algorithm, including the necessary input files and parameters (this may be needed if you want to combine your own matching algorithm, with steps 3, 4 and 5 of FALDOI). 
Moreover, in both cases you can decide to run only some of the algorithm's steps (the python scripts set boolean variables for each step). This can be useful to avoid recomputing matches several times if you plan to run the minimization with different parameters.

### C++ executables
Given a text file with the input images paths (e.g.: 'sintel_one_frame_easy.txt' in [example_data](example_data/final/)) you can obtain the output flow by following the [Algorithm's steps](#Algorithm's-steps) and calling the executables as follows:

#### Compute matches
* With SIFT:

Computing descriptors (once per image: i0 + i1)
```bash
./sift_cli im_name0.png -nspo 15 > im_name0_descriptors.txt
```
Computing matches (once forward i0 => i1, once backward i0 <= i1)
```bash
./match_cli im_name0_descriptors.txt im_name1_descriptors.txt > im_name0_matches.txt
```
* With DeepMatching:

Computing matches (once forward i0 => i1, once backward i0 <= i1)
```bash
./deepmatching im_name0.png im_name1.png -nt 4 -downscale 2
```
To see a more in detail usage of the *sift_cli* , *match_cli* (for SIFT) and *deepmatching* (DeepMatching) executables, visit the [SIFT anatomy][4] and [DeepMatching][5] pages or/and check their source code README.md files. In order to generate the needed binaries, you MUST download the source code from the links above and compile them as follows:

* SIFT: download and unzip the compressed file in the desired directory. Change to the directory where the *Makefile* is located in a terminal. Simply run the command *'make'* which should yield the *sift_cli* , *match_cli* binaries.
* DeepMatching: download and unzip the compressed file in the desired directory. Open a terminal and change to the directory where the *Makefile* is located. In this case, we needed to change some lines of the original *Makefile* from the source code in order to successfully compile it under Linux. Just comment out the lines inside the *'ifeq ($(OS_NAME),Linux)'* statement for *'LAPACKLDFLAGS=-latlas -lblas'*. If this does not work, you may try the original lines as it may work in your case.


Please note that in both cases, you may need to install the required libraries (see README files from both sources for details).
Also in both cases, you will need to copy the binaries to the *'build'* directory (created after recompiling the project with *'./recompile_changes.sh'* in the root folder). <br> 
To avoid having to repeat this step every time you change the code, you could just add a line to the bash script above of the type: *'cp ../precompiled_binaries/* .'*  to copy them to the build directory after each new compilation. 

For a simplified usage of the algorithm, you can take a look at any of the Python scripts on the following section to see some usage examples with FALDOI (only some default parameters of the matching algorithms are tweaked and most of these are fixed in the code).

NOTE: *'nspo'* means the number of scales per octave; *'nt'* means the number of threads, *'downscale'* is the downscaling factor to apply to the original image.

#### Sparse flow
```bash
./sparse_flow list_matches.txt img_width img_height out_sparse_n.flo
```

#### Local faldoi
- **Without** occlusions:
```bash
./local_faldoi file_paths_to_images.txt out_sparse_1.flo out_sparse_2.flo out_local.flo sim_map.tiff [options...]
```

or (if you have input saliency files for both images)

```bash
./local_faldoi file_paths_to_images.txt out_sparse_1.flo out_sparse_2.flo out_local.flo sim_map.tiff sal0.tiff sal1.tiff [options...]
```
- **With** occlusions:
```bash
./local_faldoi file_paths_to_images.txt out_sparse_1.flo out_sparse_2.flo out_local.flo sim_map.tiff occ_loc.png [options...]
```

or (if you have input saliency files for both images)

```bash
./local_faldoi file_paths_to_images.txt out_sparse_1.flo out_sparse_2.flo out_local.flo sim_map.tiff occ_loc.png sal0.tiff sal1.tiff [options...]
```

options (python scripts have equivalent ones with similar names and longer explanation):
+ `-m (0)` &emsp;&emsp;&emsp;chooses the functional out of the following:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1       &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;0<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1_W     &emsp;&ensp;&emsp;&emsp;&nbsp;&nbsp;1<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVL1     &emsp;&ensp;&nbsp;&emsp;&emsp;&nbsp;&nbsp;2<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVL1_W   &nbsp;&emsp;&ensp;&nbsp;&nbsp;&nbsp; 3<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVCSAD     &emsp;&emsp;&emsp;&ensp;&nbsp;4<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVCSAD_W   &nbsp;&emsp;&ensp;&nbsp;&nbsp;5<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVCSAD   &nbsp;&emsp;&ensp;&nbsp;&nbsp; 6<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVCSAD_W &ensp;&nbsp; 7<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1_OCC   &emsp;&emsp;&ensp;8
+ `-wr (5)`     	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;radius value wr (5) - patch 2\*wr + 1 x 2\*wr +1 (11 x 11).
+ `-p (None)`   	&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;file of parameters (see function _init_params_ in [utils_preprocess.cpp](src/utils_preprocess.cpp) for more details).
+ `-loc_it (3)` 	&emsp;&emsp;&emsp;&emsp;&ensp;number of iterations for the local minimization.
+ `-max_pch_it (3)` 	&emsp;&emsp;&nbsp; number of iterations per patch (for each 'loc_it')
+ `-split_img (0)`     	&emsp;&emsp;&ensp;&nbsp; whether to split image into parts to boost speed.
+ `-h_parts (3)`     	&emsp;&emsp;&ensp;&emsp;&ensp;number of horizontal parts.
+ `-v_parts (2)`     	&emsp;&emsp;&emsp;&ensp;&ensp;number of vertical parts.
+ `-fb_thresh (*)`      &emsp;&emsp;&ensp;&ensp;forward-backward consistency check threshold(*) Def. 0.45 for SIFT, 13 for DeepMatching.
+ `-partial_res (0)`    &emsp;&ensp;&ensp;whether to save all intermediate flows or not.

#### Global faldoi
- **Without** occlusions:
```bash
./global_faldoi file_paths_to_images.txt out_local.flo out_final.flo [options...]
```
- **With** occlusions:
```bash
./global_faldoi file_paths_to_images.txt out_local.flo out_final.flo occ_loc.png occ_final.png [options...]
```

options:

+ `-m (0)`      	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;changes the functional (check aux_energy_model.h).
+ `-w (5)`      	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;number of warpings.
+ `-p (None)`   	&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp; file of parameters (see function _init_params_ in
[utils_preprocess.cpp](src/utils_preprocess.cpp) for more details).
+ `-glb_iters (400)`    &emsp;&nbsp;&nbsp;&nbsp; number of iterations for the global minimization. 

### Python scripts - Usage
As you saw above, calling each binary with the correct parameters and keeping track of all output files to pass them to the following step, etc. can be very convoluted. For that reason, we suggest that you try using the python scripts to simplify the process. In the directory [scripts_python](scripts_python), you will find three main scripts that execute the whole algorithm following all the steps detailed in the [Algorithm's steps section](#algorithm's-steps).

#### faldoi_sift.py
Given a text file containing the paths to the input frames, computes the optical flow estimation based on the FALDOI algorithm's energy minimisation. Matches are extracted with SIFT. Usage:
```bash
./faldoi_sift.py path2imgs.txt [options]
```
options: 

+ `-vm (0)`    &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&emsp;&ensp;changes the functional (check aux_energy_model.h).
+ `-wr (5)`    &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&emsp;&ensp;windows radius or patch size (2\*wr + 1 x 2\*wr + 1). E.g.: wr=5 means a 11x11 patch size.
+ `-local_iter` &emsp;&emsp;&ensp;&emsp;&emsp;number of iterations for the local minimization.
+ `-patch_iter` &emsp;&emsp;&ensp;&emsp;&emsp;number of iterations per patch (for each 'local_iter')
+ `-split_img (0)`     	&emsp;&emsp;&ensp;&nbsp; whether to split image into parts to boost speed.
+ `-h_parts (3)`     	&emsp;&emsp;&ensp;&emsp;&ensp;number of horizontal parts.
+ `-v_parts (2)`     	&emsp;&emsp;&emsp;&ensp;&ensp;number of vertical parts.
+ `-fb_thresh (*)`      &emsp;&emsp;&ensp;&ensp;forward-backward consistency check threshold(*) Def. 0.45 for SIFT, 13 for DeepMatching.
+ `-partial_res (0)`    &emsp;&ensp;&ensp;whether to save all intermediate flows or not.
+ `-warps (5)` &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;numbers of warps at the finest scale (global minimisation).
+ `-glob_iter (400)` &emsp;&ensp;&ensp;number of iterations for the global minimization. 
+ `-nsp (15)`  &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&ensp;number of scales per octave to be computed by the SIFT algorithm.
+ `-res_path`  &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&ensp;path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.

#### faldoi_deep.py
Does the same as the above script but the matches are extracted with the DeepMatching algorithm instead of SIFT. Usage:
```bash
./faldoi_deep.py path2imgs.txt [options]
```
options:

+ `-vm (0)`     &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;changes the functional (check aux_energy_model.h).
+ `-wr (5)`     &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;windows radius or patch size (2\*wr + 1 x 2\*wr + 1). E.g.: wr=5 means a 11x11 patch size.
+ `-local_iter` &emsp;&emsp;&ensp;&emsp;&emsp;number of iterations for the local minimization.
+ `-patch_iter` &emsp;&emsp;&ensp;&emsp;&emsp;number of iterations per patch (for each 'local_iter')
+ `-split_img (0)`     	&emsp;&emsp;&ensp;&nbsp; whether to split image into parts to boost speed.
+ `-h_parts (3)`     	&emsp;&emsp;&ensp;&emsp;&ensp;number of horizontal parts.
+ `-v_parts (2)`     	&emsp;&emsp;&emsp;&ensp;&ensp;number of vertical parts.
+ `-fb_thresh (*)`      &emsp;&emsp;&ensp;&ensp;forward-backward consistency check threshold(*) Def. 0.45 for SIFT, 13 for DeepMatching.
+ `-partial_res (0)`    &emsp;&ensp;&ensp;whether to save all intermediate flows or not.
+ `-warps (5)` &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;numbers of warps at the finest scale (global minimisation).
+ `-glob_iter (400)` &emsp;&ensp;&ensp;number of iterations for the global minimization. 
+ `-warps (5)`  &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;numbers of warps at the finest scale (global minimisation).
+ `-th (0.45)`  &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;threshold to discard outliers from DeepMatching.
+ `-nt (4)`     &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;number of CPUs used to compute DeepMatching.
+ `-downscale (2)` &emsp;&emsp;&emsp;related to the scale at which DM will work (by default, half original resolution).
+ `-max_scale (sqrt(2))` maximum scale of DeepMatching.
+ `-rot_plus (45)` &emsp;&emsp;&emsp;positive rotation angle for DeepMatching.
+ `-rot_minus (45)` &emsp;&emsp;&ensp;negative rotation angle for DeepMatching.
+ `-res_path`   &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.


#### NOTE
#### faldoi_deep_occ.py
Includes the optional parameters to model occlusions (only available with the TVL1 energy functional right now). Matches are computed with Deep Matching. Usage
```bash
./fast_faldoi_occ.py file_paths_to_images.txt [options]
```
options: same as fast_faldoi.py (see above).

#### REMARKS
1\. Out of the list of available functionals listed above (see [C++ Executables/Local faldoi](#local-faldoi)), in the code published in the IPOL journal all of them are available with the exception of the TVL1+occlusions. If you wish to test that particular function, please refer to the Github repository linked below. For more details about how the occlusions are model for FALDOI, please read the base article where the [occlustion estimation model](https://www.dtic.upf.edu/~cballester/DAGM2012.pdf) was presented.

2\. The ability of working with partitions is currently **only compatible with the TVL2-L1 functional** (the one that benefits more from the speed-up as commented in the IPOL paper).


3\. The format of the **path2imgs.txt** file, is the following:
```bash
path/to/frame_0001.png
path/to/frame_0002.png
```
For now, the paths should be absolute or relative to the current path (using the 'tilde' character does not work at the moment).

Notice that the order in the text file *is* important. The paths should be entered (one per line) as follows: first line: I0, second line: I1. Chronologically, they are ordered as: I0 (t) --> I1 (t+1).
Since we are not modeling occlusions (see GitHub link [below](#Bugs,-issues-or-questions) for the complete version), we can pass only 2 paths to frames I0 and I1. Nevertheless, specifying the paths to the 4 consecutive frames as if you were using occlusions will enable you to use the same input files if you desire to start modeling occlusions at some point.

#### Example of usage
With the source code, you can already run the algorithm with some sample frames stored in folders inside the [example_data](/example_data) folder. Most of these frames are part of a bigger dataset known as [MPI Sintel][6]. In the folder, two directories have been created, each containing the first 4 frames of the same sequence ('alley_1'). 

The 'clean' directory contains the synthetic frames without any defussion or distortion applied to them. The second directory, 'final', applies some transformations to the original frames which makes the optical flow estimation more challenging.

If you want to run the algorithm with SIFT matches and specify your own results path, you just need to navigate to the [scripts_python](scripts_python/) folder and execute the following line in your terminal:
```bash
./faldoi_sift.py path2imgs.txt -res_path ../tmp_faldoi/Experiment1/Results/
```
Remember to add the final slash '/' so the files are created *inside* the child folder (in the example, 'Results') and not in its parent directory. Finally, you may run the '*faldoi_deep.py*' script in a similar fashion.

The other folders contain special cases tested throughout the algorithm's development and optimization. All the sequences with known ground truth include a subfolder called *'gt'* with extra subfolders for occlusions and invalid pixel masks (so one can compute error metrics). We follow the MPI-Sintel naming nomenclature for the folders (you can find more information in [Sintel's downloads' website](http://sintel.is.tue.mpg.de/downloads)). More images can be found by downloading the supplementary material submitted to the IPOL article (see link at the top of this file).

## Bugs, issues or questions
If you encounter any bugs, issues or have any questions about this source code or the implemented algorithm, please, do not hesitate to open an issue in this repository. Alternatively, you may want to send an e-mail to the [last active developer](mailto:fperez.gamonal@gmail.com). The complete code can be found on [GitHub](https://github.com/fperezgamonal/faldoi-ipol/) including some extra files left out of the IPOL publication for the sake of simplicity and coherence with the accompanying article.


## Developers
- Roberto Palomares, 2014 (main dev)
- Onofre Martorell, 2017 (MsC Thesis: migrated most of the code to C++, added occlusions estimation for the TVL2-L1 functional and fixed bugs)
- Ferran Pérez, 2018-2019 (optimise the code for IPOL, fixed some bugs with some functionals. All the work is summarised in the [following article](http://www.ipol.im/))

## License and copyright
This software is licensed under the BSD 3-Clause license. For details, see [LICENSE.md](LICENSE.md/LICENSE.md)<br>
Copyright &copy; 2014, Roberto P.Palomares  &emsp; *r.perezpalomares@gmail.com*<br>
Copyright &copy; 2017, Onofre Martorell  &emsp; &emsp; &emsp;*onofremartorelln@gmail.com*<br>
Copyright &copy; 2018-2019, Ferran Pérez &emsp;&emsp; [*fperez.gamonal@gmail.com*](mailto:fperez.gamonal@gmail.com)<br>
All rights reserved.


[1]: https://link.springer.com/content/pdf/10.1007%2Fs10851-016-0688-y.pdf "original Springer paper"
[2]: https://arxiv.org/abs/1602.08960v3 "arXiv paper"
[3]: http://www.ipol.im/ "IPOL paper"
[4]: http://www.ipol.im/pub/art/2014/82/ "SIFT implementation"
[5]: http://lear.inrialpes.fr/src/deepmatching/ "DeepMatching"
[6]: http://sintel.is.tue.mpg.de/ "MPI Sintel Database"

