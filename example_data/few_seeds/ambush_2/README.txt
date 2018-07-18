# Frame(s) structure
	* frame_0020.png ==> i0 (source frame, the one we estimate the optical flow from at time 't')
	* frame_0021.png ==> i1 (second frame, the one after i0 at time 't+1', used to estimate the flow w.r.t i0)

# There are no frame i_1 (frame at 't-1') nor i_2 (frame at 't+2') since the sequence has no more frames available.

# Extra frames are used to estimate forward and backward occlusions.

