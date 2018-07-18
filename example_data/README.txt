========================= example_data/ =========================
This folder contains some of the data used to test the FALDOI algorithm.

======== FOLDERS ========
Currently, as of August 2018, there are 4 folders, which are:

	- 'clean': named following the nomenclature of MPI-Sintel(at http://sintel.is.tue.mpg.de/downloads), this folder contains three subfolders: 'easy', 'medium' and 'hard' where three sequences of 4 frames are presented. The first one, 'easy' is extracted from MPI-Sintel clean alley_1 sequence and contains frames_0001-0004.png; the second one, 'medium' is extracted from clean cave_2 sequence and contains frames_0001-0004.png; finally, the 'hard' folder contains frames_0001-0004.png of the market_5 sequence.

	- 'final': also named following the nomenclature of MPI-Sintel, this folder contains the same frames as 'clean' but extracted from the final subset of sintel (i.e.: including motion blur, noise and other effects).

	- 'few_seeds': this folder was created to analyse the robustness of the algorithm to a very small number of seeds. Inside it you will find some frames from MPI-Sintel that yield very few seeds for SIFT or equivalently, a lot of outliers for DeepMatching.

	- 'large_displacements': this folder was created to analyse the algorithm's ability to detect large displacements of relatively small objects. Inside it you will find two sequences: one composite image and one real image (from the paper of LDOF at: https://www.cs.cmu.edu/~katef/LDOF.html).
	- Other folders may have been added since the last modification of this README file.

Note: The first sequence was already analysed in the first publication of FALDOI in Springer at https://link.springer.com/article/10.1007/s10851-016-0688-y, but we could not find the original ground truth flo file to report any analytical metrics on the sequence (you can find the png version of the ground truth flow, faldoi's and also LDOF' for a visual interpretation). The second sequence does not have a ground truth flow either.

Note: you can get further information about a specific sequence by opening the README.txt file inside one specific folder.

======== FOLDERS CONTENTS ========
The first 3 folders listed above contain in each subfolder (corresponding to one sequence) follow this structure:
	-----FOLDER----/SUBF1/SUBF2/			CONTENTS
	$SEQUENCE_FOLDER/				contains the original png frames to analyse
	$SEQUENCE_FOLDER/gt/				" the ground truth files (.flo + .png) " "
	$SEQUENCE_FOLDER/gt/occlusions/			" the occlusions mask to compute EPE-unmatched
	$SEQUENCE_FOLDER/gt/invalid/			" the invalid mask that defines pixels w. invalid flow

======== USAGE ========
To test any of the former sequences, you only need to pass the corresponding '.txt' file located at the parent folder of the desired sequence (e.g.: for final, medium sequence, the txt file is at 'example_data/final/', one level above 'example_data/final/medium/') as the first argument to faldoi.

For example, if you want to analyse the 'hard' sequence of the final pass of MPI-Sintel using FALDOI with Deep Matching seeds, a window size of 9x9 (corresponding to a radius=4) and the TVL1 energy functional (number '0'), you can use the command (from the 'scripts_python' folder):

	./faldoi_deep.py ../example_data/final/sintel_one_frame_hard.txt -wr 4 -vm 0 [extra params]

For more information about the parameters of the scripts 'faldoi_sift.py', 'faldoi_deep.py' and 'faldoi_deep_occ.py', check the README.txt located at the 'scripts_python' subfolder. In addition, more examples of usage are showcased in the README.md on the root directory 'faldoi-ipol_[version_number]'.
	
