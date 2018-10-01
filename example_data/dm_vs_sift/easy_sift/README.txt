# Frame(s) structure
# One frames is analysed for the MPI-Sintel 'bandage_1' training/final sequence

# first frame: frame_0035.png
	* frame_0035.png ==> i0 (source frame, the one we estimate the optical flow from at time 't')
	* frame_0036.png ==> i1 (second frame, the one after i0 at time 't+1', used to estimate the flow w.r.t i0)
	* frame_0034.png ==> i_1 (previous frame to i0, at time 't-1', used to compute occlusions)
	* frame_0037.png ==> i_2 (following frame to i1, at time 't+2', used to compute occlusions)

# Extra frames (frame_0034.png and frame_0037.png) are used to estimate forward and backward occlusions (if the functional models them).

