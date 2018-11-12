========================= example_data/ =========================
This folder contains data provided to test the FALDOI algorithm.

======== FOLDERS ========
For simplicity, we only provide two sequences of images for testing. Both of them belong to the MPI Sintel synthetic dataset, one off the clean pass and one from the final pass.

In more detail:
	- 'clean': named following the nomenclature of MPI-Sintel(at http://sintel.is.tue.mpg.de/downloads), this folder contains three subfolders: 'easy', 'medium' and 'hard' where three sequences of 4 frames are presented. The first one, 'easy' is extracted from MPI-Sintel clean alley_1 sequence and contains frames_0001-0004.png; the second one, 'medium' is extracted from clean cave_2 sequence and contains frames_0001-0004.png; finally, the 'hard' folder contains frames_0001-0004.png of the market_5 sequence.

	- 'final': also named following the nomenclature of MPI-Sintel, this folder contains the same frames as 'clean' but extracted from the final subset of sintel (i.e.: including motion blur, noise and other effects).

You can get more sample frames by downloading the supplementary material provided alongside the article published in ipol.im (search for 'FALDOI').

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

For example, if you want to analyse the 'final' sequence of the clean pass of MPI-Sintel using FALDOI with Deep Matching seeds, a window size of 9x9 (corresponding to a radius=4) and the TVL1 energy functional (number '0'), you can use the command (from the 'scripts_python' folder):

	./faldoi_deep.py ../example_data/clean/sintel_one_frame_easy.txt -wr 4 -vm 0 [extra params]

For more information about the parameters of the scripts 'faldoi_sift.py', 'faldoi_deep.py' and 'faldoi_deep_occ.py', check the README.txt located at the 'scripts_python' subfolder. In addition, more examples of usage are showcased in the README.md on the root directory 'faldoi-ipol_[version_number]'.
	
