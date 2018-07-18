# Frame(s) structure
# Two frames are analysed for the MPI-Sintel 'ambush_4' training/final sequence (that is why you see more than 4 frames in this folder)

# a) first frame: frame_0004.png

	* frame_0004.png ==> i0 (source frame, the one we estimate the optical flow from at time 't')
	* frame_0005.png ==> i1 (second frame, the one after i0 at time 't+1', used to estimate the flow w.r.t i0)
	* frame_0003.png ==> i_1 (previous frame to i0, at time 't-1', used to compute occlusions)
	* frame_0006.png ==> i_2 (following frame to i1, at time 't+2', used to compute occlusions)

	# Extra frames (frame_0003.png and frame_0006.png) are used to estimate forward and backward occlusions (if the functional models them).

# b) second frame: frame_0006.png

	* frame_0006.png ==> i0 (source frame, the one we estimate the optical flow from at time 't')
	* frame_0007.png ==> i1 (second frame, the one after i0 at time 't+1', used to estimate the flow w.r.t i0)
	* frame_0005.png ==> i_1 (previous frame to i0, at time 't-1', used to compute occlusions)
	* frame_0008.png ==> i_2 (following frame to i1, at time 't+2', used to compute occlusions)

	# Extra frames (frame_0005.png and frame_0008.png) are used to estimate forward and backward occlusions (if the functional models them).

