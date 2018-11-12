* Auxiliar binaries for the matchers *

In this directory, three files are included, belonging to two different matchers:
	SIFT
		- 'sift_cli':		computes descriptors of input image following the SIFT algorithm.
		- 'match_cli':		takes the descriptors generated above for two frames and computes the matches.
	The implementation used is the one from:


	DeepMatching
		- 'deepmatching':	computes DeepMatching matches between two input frames.
	The implementation used is the one from: 


	NOTES:
	If the binaries included with the source files do not work, you can download the newest version online from the link above and try again.
	Notice that the error may due to your system being different from the one the binaries (only for SIFT) were compiled in.
	In that case, after you download and compile the matchers code for your architecture, you should run '. recompile_changes.sh' from the root directory to apply the changes and be able to execute the program again.
