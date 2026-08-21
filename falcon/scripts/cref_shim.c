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

/* --- expose the real sampler with a PRNG we seed directly ---------------- */
static sampler_context g_spc;

void
shim_sampler_init(const uint8_t *state56, unsigned logn)
{
	memcpy(g_spc.p.state.d, state56, 56);
	Zf(prng_refill)(&g_spc.p);
	g_spc.sigma_min = fpr_sigma_min[logn];
}

int
shim_sampler(double mu, double isigma)
{
	return Zf(sampler)(&g_spc, FPR(mu), FPR(isigma));
}

double
shim_sigma_min(unsigned logn)
{
	return fpr_sigma_min[logn].v;
}

/* how many bytes of the prng buffer have been consumed so far */
size_t
shim_sampler_ptr(void)
{
	return g_spc.p.ptr;
}

/* ffSampling_fft driven by the REAL sampler off the shim's prng, so the
   order and count of randomness consumption is part of what is compared. */
void
shim_ffsampling_real(fpr *z0, fpr *z1, const fpr *tree,
	const fpr *t0, const fpr *t1, unsigned logn, fpr *tmp)
{
	static fpr tt0[1024], tt1[1024];
	size_t n = (size_t)1 << logn;
	memcpy(tt0, t0, n * sizeof *t0);
	memcpy(tt1, t1, n * sizeof *t1);
	ffSampling_fft(Zf(sampler), &g_spc, z0, z1, tree, tt0, tt1, logn, tmp);
}

/* --- end to end: sign with an expanded key and a PRNG we seed ------------ */
/* Mirrors Zf(sign_tree) but takes the 56-byte ChaCha20 state directly, so the
   Julia side can drive the identical stream. */
int
shim_sign_tree(int16_t *sig, const uint8_t *state56,
	const fpr *expanded_key, const uint16_t *hm, unsigned logn, uint8_t *tmp)
{
	sampler_context spc;
	fpr *ftmp = (fpr *)tmp;
	spc.sigma_min = fpr_sigma_min[logn];
	memcpy(spc.p.state.d, state56, 56);
	Zf(prng_refill)(&spc.p);
	for (;;) {
		if (do_sign_tree(Zf(sampler), &spc, sig, expanded_key, hm, logn, ftmp))
			return 1;
	}
}

/* Compress s2 exactly as the reference's signature encoder does. */
size_t
shim_comp_encode(void *out, size_t max_out_len, const int16_t *x, unsigned logn)
{
	return Zf(comp_encode)(out, max_out_len, x, logn);
}
