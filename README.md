[![Python 2.7](https://img.shields.io/badge/python-2.7-green.svg)](https://www.python.org/downloads/release/python-270/)
[![Python 3.5](https://img.shields.io/badge/python-3.5-green.svg)](https://www.python.org/downloads/release/python-350/)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# FALDOI-IPOL
Stems from the basic [FALDOI: A New Minimization Strategy for Large Displacement Variational Optical Flow](https://link.springer.com/content/pdf/10.1007%2Fs10851-016-0688-y.pdf) by Roberto P. Palomares, Enric Meinhardt-Lopis, Coloma Ballester and Gloria Haro algorithm and aims to add occlusion estimation to several energy functionals and optimise the code to be published on the [IPOL journal](http://www.ipol.im/) with an interactive demo.

## Paper(s) and citation
If you use FALDOI, please cite _any_ of the following papers:

[Springer original paper (November 2016)](https://link.springer.com/content/pdf/10.1007%2Fs10851-016-0688-y.pdf):
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
the [arXiv paper (September 2016)](https://arxiv.org/abs/1602.08960v3):
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
You can reference this implementation with: [IPOL article, demo and supplementary material (August 2018)](_doi_to_ipol_to__be__updated_):
```
bibTeX for IPOL
``` 
___
## Table of contents
* [Getting Started](#getting-started)
  * [Pre-requisites](#pre-requisites)
  * [Compilation](#compilation)
  * [Algorithm's steps](#algorithm's-steps)
* [Execution](#execution)
  * [C++ executables - Usage](#c-executables---usage)
  * [Python scripts - Usage](#python-scripts---usage)
    * [Example of usage](#example-of-usage)
* [Bugs, issues and questions](#bugs,-issues-and-questions)
* [Developers](#developers)
* [License and copyright](#license-and-copyright)
___
## Getting Started
The code is written in C/C++ and includes Python scripts to unify and simplify all the tasks the algorithms performs.
### Pre-requisites
The software needs the following programs to be installed in order to function:
- Python 3.5
- Pillow for Python 3 (if you do not want to use Pillow, in the python script there is a commented section above the PIL import that uses [_imagemagick_](https://www.imagemagick.org/script/identify.php) instead). Pillow is faster (since it is built-in python).
- libpng (included)
- OpenMP (Optional but recommended). Should be included with your compiler if you have a relatively new version (e.g.: gcc supports it since version 4.2).
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
### Algorithm's steps
In order to obtain the final optical flow estimation, the algorithm follows these steps:
1. Reads a _.txt_ file containing the paths to the input frames to be processed (see [Python scripts - Usage/NOTE](#python-scripts---usage) for details)
2. Computes the descriptors and matches between frames I0 and I1 (with SIFT or Deepmatching).
3. Computes a sparse flow from the matches.
4. Computes a dense flow from the initial sparse flow (local step).
5. Computes global minimization taking as input the flow of the previous step (global step).

## Execution
Once you have successfully compiled the program, you can execute the algorithm in two ways: using the C++ executables or running the Python scripts (*suggested*). The advantage of choosing the latter is that all the algorithm's steps are run at once sequentially, while with C++ executables, you will have to manually call each block of the algorithm, including the necessary input files (this may be needed if you want to combine your own matching algorithm, with steps 3, 4 and 5 of FALDOI). 
In both cases, the execution varies if you want to include occlusions or not. Moreover, in both cases you can decide to run only some of the algorithm's steps (the python scripts set boolean variables for each step). This can be useful to avoid recomputing matches several times if you plan to run the minimization with different parameters.

### C++ executables - Usage
Given a text file with the input images paths (e.g.: 'sintel_one_frame_easy.txt' in [example_data](example_data/final/)) you can obtain the output flow by following the [Algorithm's steps](#Algorithm's-steps) and calling the executables as follows:
#### Compute matches
- With SIFT\
Computing descriptors (once per image: i0 + i1)
```bash
./sift_cli im_name0.png -nspo 15 > im_name0_descriptors.txt
```
Computing matches (once forward i0=>i1, once backward i0<==i1)
```bash
./match_cli im_name0_descriptors.txt im_name1_descriptors.txt > im_name0_matches.txt
```
- With DeepMatching:
```bash
./deepmatching im_name0.png im_name1.png -nt 4 -downscale 2
```
To see a more in detail usage of the _sift_cli_ , _match_cli_ (for SIFT) and _deepmatching_ (DeepMatching) executables, visit the [SIFT anatomy](http://www.ipol.im/pub/art/2014/82/) and [DeepMatching](http://lear.inrialpes.fr/src/deepmatching/) pages or/and check their source code's README.md files.
Alternatively, you can take a look at any of the Python scripts on the following section to see some usage examples with FALDOI (only some default parameters of the matching algorithms are tweaked and most of these are fixed in the code).

NOTE: 'nspo' means the number of scales per octave; 'nt' means the number of threads, 'downscale' is the downscaling factor to apply to the original image.

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
+ `-m (0)`      &emsp;&emsp;&emsp;chooses the functional out of the following:\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1       &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;0\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1_W     &emsp;&ensp;&emsp;&emsp;&nbsp;&nbsp;1\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVL1     &emsp;&ensp;&nbsp;&emsp;&emsp;&nbsp;&nbsp;2\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVL1_W   &nbsp;&emsp;&ensp;&nbsp;&nbsp;&nbsp; 3\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVCSAD     &emsp;&emsp;&emsp;&ensp;&nbsp;4\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVCSAD_W   &nbsp;&emsp;&ensp;&nbsp;&nbsp;5\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVCSAD   &nbsp;&emsp;&ensp;&nbsp;&nbsp; 6\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_NLTVCSAD_W &ensp;&nbsp; 7\
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;M_TVL1_OCC   &emsp;&emsp;&ensp;8
+ `-wr (5)`     	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;radius value wr (5) - patch 2\*wr + 1 x 2\*wr +1 (11 x 11).
+ `-p (None)`   	&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;file of parameters (see function _init_params_ in [utils_preprocess.cpp](src/utils_preprocess.cpp) for more details).
+ `-loc_it (3)` 	&emsp;&emsp;&emsp;&emsp;&ensp;number of iterations for the local minimization.
+ `-max_pch_it (3)` 	&emsp;&emsp;&nbsp; number of iterations per patch (for each 'loc_it')
+ `-split_img (1)`     	&emsp;&emsp;&ensp;&nbsp; whether to split image into parts to boost speed.
+ `-h_parts (3)`     	&emsp;&emsp;&ensp;&emsp;&ensp;number of horizontal parts.
+ `-v_parts (2)`     	&emsp;&emsp;&emsp;&ensp;&ensp;number of vertical parts.

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
./faldoi_sift.py file_paths_to_images.txt [options]
```
options: 
+ `-vm (0)`    &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&emsp;&ensp;changes the functional (check aux_energy_model.h).
+ `-wr (5)`    &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&emsp;&ensp;windows radius or patch size (2\*wr + 1 x 2\*wr + 1). For instance, wr=5 means a 11x11 patch size.
+ `-warps (5)` &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;numbers of warps at the finest scale (global minimisation).
+ `-nsp (15)`  &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&ensp;number of scales per octave to be computed by the SIFT algorithm.
+ `-res_path`  &emsp;&emsp;&ensp;&ensp;&emsp;&emsp;&ensp;path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.

#### faldoi_deep.py
Does the same as the above script but the matches are extracted with the DeepMatching algorithm instead of SIFT. Usage:
```bash
./faldoi_deep.py file_paths_to_images.txt [options]
```
options:
+ `-vm (0)`     &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;changes the functional (check aux_energy_model.h).
+ `-wr (5)`     &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;windows radius or patch size (2*wr + 1 x 2*wr + 1). For instance, wr=5 means a 11x11 patch size.
+ `-warps (5)`  &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;numbers of warps at the finest scale (global minimisation).
+ `-th (0.45)`  &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;threshold to discard outliers from DeepMatching.
+ `-res_path`   &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.

#### faldoi_deep_occ.py
Includes the optional parameters to model occlusions (only available with the TVL1 energy functional right now). Matches are computed with Deep Matching. Usage
```bash
./fast_faldoi_occ.py file_paths_to_images.txt [options]
```
options: same as fast_faldoi.py (see above).

#### NOTE
The format of the **file_paths_to_images.txt** file, is the following:
```bash
path/to/frame_0002.png
path/to/frame_0003.png
path/to/frame_0001.png
path/to/frame_0004.png
```
For now, the paths should be absolute or relative to the current path (using the 'tilde' character does not work at the moment).

Notice that the order in the text file *is* important. The paths should be entered (one per line) as follows: first line: I0, second line: I1, third line: I_1 and fourth line: I2. Chronologically, they are ordered as: I_1 (t-1) --> I0 (t) --> I1 (t+1) --> I2 (t+2).
If you are not modeling occlusions, you can pass only 2 paths to frames I0 and I1. Nevertheless, specifying the 4 paths as if you were using occlusions will enable you to use the same input files if you desire to start modeling occlusions at some point.

#### Example of usage
With the source code, you can already run the algorithm with some sample frames stored in folders inside the [example_data](/example_data) folder. Most of these frames are part of a bigger dataset known as [MPI Sintel](http://sintel.is.tue.mpg.de/). In the folder, two directories have been created, each containing the first 4 frames of the same sequence ('alley_1'). 

The 'clean' directory contains the synthetic frames without any defussion or distortion applied to them. The second directory, 'final', applies some transformations to the original frames which makes the optical flow estimation more challenging.

If you want to run the algorithm with SIFT matches and specify your own results path, you just need to navigate to the [scripts_python](scripts_python/) folder and execute the following line in your terminal:
```bash
./fast_sift.py file_paths_to_images.txt -vm 0 -wr 5 -res_path ../tmp_faldoi/Experiment1/Results/
```
Remember to add the final slash '/' so the files are created _inside_ the child folder (in the example 'Results') and not in its parent directory. Finally, you may run the '*faldoi_deep.py*' and '*faldoi_deep_occ*' scripts in a similar fashion (the last has not been tested thoroughly so some bugs may be present).

The other folders contain special cases tested throughout the algorithm's development and optimization. All the sequences with known ground truth include a subfolder called _'gt'_ with extra subfolders for occlusions and invalid pixel masks (so one can compute error metrics). We follow the MPI-Sintel naming nomenclature for the folders (you can find more information in [Sintel's website](http://sintel.is.tue.mpg.de/downloads)).

## Bugs, issues or questions
If you encounter any bugs, issues or have any questions about this source code or the implemented algorithm, please, do not hesitate to open an issue in this repository. Alternatively, you may want to send an e-mail to the [last active developer](mailto:fperez.gamonal@gmail.com).


## Developers
- Roberto Palomares, 2014 (main dev)
- Onofre Martorell, 2017 (MsC Thesis: migrated most of the code to C++, added occlusions estimation for the TVL1 functional and fixed bugs)
- Ferran Pérez, 2018 (optimise the code for IPOL, fixed some bugs with some functionals. All the work is summarised in the [following article](http://www.ipol.im/))

## License and copyright
This software is licensed under the BSD 3-Clause license. For details, see [LICENSE.md](LICENSE.md/LICENSE.md)\
Copyright &copy; 2014, Roberto P.Palomares _r.perezpalomares@gmail.com_\
Copyright &copy; 2017, Onofre Martorell _onofremartorelln@gmail.com_\
Copyright &copy; 2018, Ferran Pérez _fperez.gamonal@gmail.com_\
All rights reserved.

