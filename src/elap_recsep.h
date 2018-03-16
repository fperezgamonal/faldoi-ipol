#ifndef ELAP_RECSEP_H
#define ELAP_RECSEP_H

void elap_recursive_separable(float *out, float *in, int w, int h, int pd,
		float timestep, int niter, int scale);

#endif// ELAP_RECSEP_H