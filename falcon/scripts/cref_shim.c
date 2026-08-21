/* Expose the static ffSampling_fft with a deterministic sampler, so the
   recursion's arithmetic can be compared without the PRNG in the way.
   The "sampler" is plain floor(): trivially reproducible on the Julia side. */
#include "sign.c"

static int
shim_floor(void *ctx, fpr mu, fpr isigma)
{
	(void)ctx; (void)isigma;
	return (int)fpr_floor(mu);
}

void
shim_ffsampling(fpr *z0, fpr *z1, const fpr *tree,
	const fpr *t0, const fpr *t1, unsigned logn, fpr *tmp)
{
	static fpr tt0[1024], tt1[1024];
	size_t n = (size_t)1 << logn;
	memcpy(tt0, t0, n * sizeof *t0);
	memcpy(tt1, t1, n * sizeof *t1);
	ffSampling_fft(shim_floor, NULL, z0, z1, tree, tt0, tt1, logn, tmp);
}
