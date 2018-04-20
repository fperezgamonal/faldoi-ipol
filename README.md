# FALDOI-IPOL
Stems from the basic FALDOI algorithm and aims to add occlusion estimation to several energy functionals and optimise the code to be published on the [IPOL journal](http://www.ipol.im/) with an interactive demo.\
Paper: [FALDOI: A New Minimization Strategy for Large Displacement Variational Optical Flow](https://arxiv.org/abs/1602.08960) by Roberto P. Palomares, Enric Meinhardt-Lopis, Coloma Ballester and Glòria Haro.
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
* [Developers](#developers)
* [License](#license)
___
## Getting Started
The code is written in C/C++ and includes Python scripts to unify and simplify all the tasks the algorithms performs.
### Pre-requisites
The software needs the following programs to be installed in order to function:
- Python 3.5
- libpng (included)
- OpenMP (should be included with your compiler if you have a relatively new version). If you run int
### Compilation
The code needs to be compiled in order to generate the needed executables. Nevertheless, if you do not want to compile the code before testing it once, we added already-compiled executables that are already linked (see [_Execution_](#Execution)).

To compile the code, we have added a new bash script that runs all the needed commands to remove the old _build_ directory, create a new one and re-compile the whole project, integrating all the changes done to any C/C++ file. To run such script, just do the following (from the project's root directory):
```bash
source recompile_changes.sh
```
or
```bash
. recompile_changes.sh
```
### Algorithm's steps
In order to obtain the final optical flow estimation, the algorithm goes through the next steps:
1. Reads a _.txt_ file containing the paths to the input frames.
2. Computes the descriptors and matches between frames I0 and I1 (with SIFT or Deepmatching).
3. Computes a sparse flow from the matches.
4. Computes a dense flow from the initial sparse flow (local step).
5. Minimizes globally the flow computed in the previous step (global step).

## Execution
Once you have successfully compiled the program, you can execute the algorithm in two ways: using the C++ executables or running the Python scripts. The advantage of choosing the latter is that all the algorithm's steps are run at once sequentially, while with C++ executables, you will have to manually call each block of the algorithm, including the necessary input files. In both cases, the execution varies if you want to include occlusions or not.

### C++ executables - Usage
Given the text file '' with the input images paths in [example_data](example_data), you can obtain the output flow by following the [Algorithm's steps](#Algorithm's-steps) and calling the executables as follows:
#### Compute matches
To see a more in detail usage of the _sift_cli_ , _match_cli_ (for SIFT) and _deepmatching_ (Deep Matching) executables, visit the [SIFT anatomy](http://www.ipol.im/pub/art/2014/82/) and [DeepMatching](http://lear.inrialpes.fr/src/deepmatching/) pages or/and check their source code's README.md files.
Alternatively, you can take a look at any of the Python scripts on the following section to see some usage examples with FALDOI (only some default parameters of the matching algorithms are tweaked and most of these are fixed in the code).

#### Sparse flow
```bash
./sparse_flow list_matches.txt num_cols num_rows out_sparse_n.flo
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
options:
+ `-m (0)`      &emsp;&emsp;&emsp;Changes the functional (check [aux_energy_model.h](src/aux_energy_model.h)).
+ `-wr (5)`     &emsp;&emsp;&ensp;Radius value wr (5) - patch 2\*wr + 1 x 2\*wr +1 (11 x 11). 

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
+ `-m (0)`      &emsp;&emsp;&emsp;Changes the functional (check aux_energy_model.h).
+ `-w (5)`      &emsp;&emsp;&emsp;Number of warpings.  

### Python scripts - Usage
In the directory [scripts_python](scripts_python), you will find three main scripts that execute the whole algorithm following all the steps detailed in the [Algorithm's steps section](#algorithm's-steps).
#### fast_sift.py
Given a text file containing the paths to the input frames, computes the optical flow estimation based on the FALDOI algorithm's energy minimisation. Matches are extracted with SIFT. Usage:
```bash
./fast_sift.py file_paths_to_images.txt [options]
```
options: 
+ `-vm (0)`    &emsp;&emsp;&ensp;Changes the functional (check aux_energy_model.h).
+ `-wr (5)`    &emsp;&emsp;&ensp;Windows radius or patch size (2\*wr + 1 x 2\*wr + 1). For instance, wr=5 means a 11x11 patch size.
+ `-warps (5)` &emsp;Numbers of warps at the finest scale (global minimisation).
+ `-nsp (15)`  &emsp;&ensp;Number of scales per octave to be computed by the SIFT algorithm.
+ `-m (0)`     &emsp;&emsp;&emsp;Whether to use gaussian weighting over the data term of Faldoi.
+ `-res_path`  &emsp;&ensp;Path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.

#### fast_faldoi.py
Does the same as the above script but the matches are extracted with the DeepMatching algorithm instead of SIFT. Usage:
```bash
./fast_faldoi.py file_paths_to_images.txt [options]
```
options:
+ `-vm (0)`     &emsp;&emsp;&ensp;Changes the functional (check aux_energy_model.h).
+ `-wr (5)`     &emsp;&emsp;&ensp;Windows radius or patch size (2*wr + 1 x 2*wr + 1). For instance, wr=5 means a 11x11 patch size.
+ `-warps (5)`  &emsp;Numbers of warps at the finest scale (global minimisation).
+ `-th (0.45)`  &emsp;Threshold to discard outliers from DeepMatching.
+ `-res_path`   &emsp;&ensp;Path where the output files will be saved (partial and final flow, descriptors and matches). If "None", the results are stored in the [Results](Results) folder.

#### fast_faldoi_occ.py
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
For now, the paths should be absolute or relative to the current path (only those with the tilde '~' character do not work ATM, e.g.: '~/Data/Results/Experiment1/'). We are trying to fix this as soon as possible.

Notice that the order is important. The paths should be entered (one per line) as follows: first line: I0, second line: I1, third line: I_1 and fourth line: I2. Chronologically, they are ordered as: I_1 (t-1) --> I0 (t) --> I1 (t+1) --> I2 (t+2).
If you are not modeling occlusions, you can pass only 2 paths to frames I0 and I1. Nevertheless, specifying the 4 paths as if you were using occlusions will enable you to use the same input files if you desire to start modeling occlusions at some point.

#### Example of usage
With the source code, you can already run the algorithm with some sample frames stored in the [example_data](/example_data) folder. These frames are part of a bigger dataset known as [MPI Sintel](http://sintel.is.tue.mpg.de/). In the folder, two directories have been created, each containing the first 4 frames of the same sequence ('alley_1'). 

The 'clean' directory contains the synthetic frames without any defussion or distortion applied to them. The second directory, 'final', applies some transformations to the original frames which makes the optical flow estimation more challenging.

If you want to run the algorithm with SIFT matches and specify your own results path, you just need to navigate to the [scripts_python](scripts_python/) folder and execute the following line in your terminal:
```bash
./fast_sift.py file_paths_to_images.txt -vm 0 -wr 5 -res_path ~/Documents/tmp_faldoi/Experiment1/Results/
```
Remember to add the final slash '/' so the files are created _inside_ the child folder (in the example 'Results') and not in its parent directory. You may run the '*faldoi_deep.py*' and '*faldoi_deep_occ*' in a similar fashion (we have not run several tests with them so some bugs may be present).

## Developers
- Roberto Palomares, 2014 (main dev)
- Onofre Martorell, 2017 (MsC Thesis: migrated most of the code to C++, added occlusions estimation for the TVL1 functional and fixed bugs)
- Ferran Pérez, 2018 (_ONGOING_: optimise the code for IPOL and later, model occlusions with other energy functionals)

## License
This software is licensed under the BSD 3-Clause license. For details, see [LICENSE.md](LICENSE.md/LICENSE.md)
