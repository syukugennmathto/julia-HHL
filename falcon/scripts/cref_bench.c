/*
 * Time the C reference implementation's keygen / sign / verify, so that the
 * Julia implementation has something honest to be compared against.
 *
 * Usage:  cref_bench <logn> <keygen_iters> <sign_iters> <verify_iters>
 * Output: one "op n median_ms mean_ms min_ms iters" line per operation, plus
 *         a header line recording which floating-point back end was compiled
 *         in -- FPEMU and FPNATIVE are not the same program and it would be
 *         dishonest to report a number without saying which one ran.
 *
 * Build (from the directory holding the reference .c files):
 *
 *     cc -O2 -o cref_bench cref_bench.c \
 *        codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c \
 *        sign.c vrfy.c -lm
 *
 * Notes on making the comparison fair:
 *
 *   - `falcon_sign_dyn` is used, not `falcon_sign_tree`, because our
 *     `falcon_sign` rebuilds the Falcon tree on every call.  Comparing
 *     sign_tree against our sign would be comparing an expanded key against an
 *     unexpanded one.  `falcon_expand_privkey` + `falcon_sign_tree` is timed
 *     separately, since that is the split our `expand_privkey` also offers.
 *
 *   - Every operation gets a fresh RNG draw but the same key, so the numbers
 *     are per-call and not amortised over one lucky key.
 *
 *   - Keygen timings are wildly variable by nature: `ntru_solve` retries until
 *     the equation is solvable and the Gram-Schmidt norm is small enough.  The
 *     median is the number to quote; the mean is printed too, and if the two
 *     are far apart, that spread *is* the result.
 *
 *   - Timing uses CLOCK_MONOTONIC.  This measures wall time, so a loaded
 *     machine inflates it; run it on an idle one.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "falcon.h"
#include "inner.h"

static double
now_ms(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

static int
cmp_double(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;

	return (x > y) - (x < y);
}

static void
report(const char *op, unsigned logn, double *t, int n)
{
	double sum = 0.0, median;
	int i;

	qsort(t, (size_t)n, sizeof *t, cmp_double);
	for (i = 0; i < n; i++) {
		sum += t[i];
	}
	median = (n & 1) ? t[n / 2] : 0.5 * (t[n / 2 - 1] + t[n / 2]);
	printf("%s %u %.6f %.6f %.6f %d\n",
		op, 1u << logn, median, sum / (double)n, t[0], n);
	fflush(stdout);
}

int
main(int argc, char *argv[])
{
	unsigned logn;
	int nkg, nsg, nvf, i;
	shake256_context rng;
	unsigned char seed[48];
	unsigned char *sk, *pk, *sig, *esk, *tmp;
	size_t sklen, pklen, siglen, esklen, tmplen;
	double *t;
	const char *msg = "falcon-jl benchmark message";
	size_t msglen = strlen(msg);

	logn = (argc > 1) ? (unsigned)atoi(argv[1]) : 9;
	nkg  = (argc > 2) ? atoi(argv[2]) : 10;
	nsg  = (argc > 3) ? atoi(argv[3]) : 100;
	nvf  = (argc > 4) ? atoi(argv[4]) : 100;
	if (logn < 1 || logn > 10) {
		fprintf(stderr, "logn must be in 1..10\n");
		return 1;
	}

	printf("# cref_bench fpemu=%d fpnative=%d sizeof_fpr=%zu\n",
		FALCON_FPEMU, FALCON_FPNATIVE, sizeof(fpr));
	printf("# op n median_ms mean_ms min_ms iters\n");

	for (i = 0; i < (int)sizeof seed; i++) {
		seed[i] = (unsigned char)i;
	}
	shake256_init_prng_from_seed(&rng, seed, sizeof seed);

	sklen  = FALCON_PRIVKEY_SIZE(logn);
	pklen  = FALCON_PUBKEY_SIZE(logn);
	esklen = FALCON_EXPANDEDKEY_SIZE(logn);
	sk  = malloc(sklen);
	pk  = malloc(pklen);
	esk = malloc(esklen);
	sig = malloc(FALCON_SIG_PADDED_SIZE(logn));

	/* One buffer big enough for every operation, allocated once: the
	   allocator is not what we are trying to measure. */
	tmplen = FALCON_TMPSIZE_KEYGEN(logn);
	if (FALCON_TMPSIZE_SIGNDYN(logn) > tmplen) {
		tmplen = FALCON_TMPSIZE_SIGNDYN(logn);
	}
	if (FALCON_TMPSIZE_EXPANDPRIV(logn) > tmplen) {
		tmplen = FALCON_TMPSIZE_EXPANDPRIV(logn);
	}
	if (FALCON_TMPSIZE_SIGNTREE(logn) > tmplen) {
		tmplen = FALCON_TMPSIZE_SIGNTREE(logn);
	}
	if (FALCON_TMPSIZE_VERIFY(logn) > tmplen) {
		tmplen = FALCON_TMPSIZE_VERIFY(logn);
	}
	tmp = malloc(tmplen);
	if (sk == NULL || pk == NULL || esk == NULL || sig == NULL
		|| tmp == NULL)
	{
		fprintf(stderr, "out of memory\n");
		return 1;
	}

	t = malloc(sizeof *t * (size_t)(nkg > nsg
		? (nkg > nvf ? nkg : nvf) : (nsg > nvf ? nsg : nvf)));
	if (t == NULL) {
		fprintf(stderr, "out of memory\n");
		return 1;
	}

	/* ---- keygen ---- */
	for (i = 0; i < nkg; i++) {
		double t0 = now_ms();

		if (falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen,
			tmp, tmplen) != 0)
		{
			fprintf(stderr, "keygen failed\n");
			return 1;
		}
		t[i] = now_ms() - t0;
	}
	report("keygen", logn, t, nkg);

	/* The last key generated above is the one everything else uses. */

	/* ---- sign (dynamic: tree rebuilt every call) ---- */
	for (i = 0; i < nsg; i++) {
		double t0;

		siglen = FALCON_SIG_PADDED_SIZE(logn);
		t0 = now_ms();
		if (falcon_sign_dyn(&rng, sig, &siglen, FALCON_SIG_PADDED,
			sk, sklen, msg, msglen, tmp, tmplen) != 0)
		{
			fprintf(stderr, "sign_dyn failed\n");
			return 1;
		}
		t[i] = now_ms() - t0;
	}
	report("sign_dyn", logn, t, nsg);

	/* ---- expand private key ---- */
	for (i = 0; i < (nsg < 20 ? nsg : 20); i++) {
		double t0 = now_ms();

		if (falcon_expand_privkey(esk, esklen, sk, sklen,
			tmp, tmplen) != 0)
		{
			fprintf(stderr, "expand_privkey failed\n");
			return 1;
		}
		t[i] = now_ms() - t0;
	}
	report("expand_privkey", logn, t, nsg < 20 ? nsg : 20);

	/* ---- sign with the expanded key ---- */
	for (i = 0; i < nsg; i++) {
		double t0;

		siglen = FALCON_SIG_PADDED_SIZE(logn);
		t0 = now_ms();
		if (falcon_sign_tree(&rng, sig, &siglen, FALCON_SIG_PADDED,
			esk, msg, msglen, tmp, tmplen) != 0)
		{
			fprintf(stderr, "sign_tree failed\n");
			return 1;
		}
		t[i] = now_ms() - t0;
	}
	report("sign_tree", logn, t, nsg);

	/* ---- verify (the signature left over from the loop above) ---- */
	for (i = 0; i < nvf; i++) {
		double t0 = now_ms();

		if (falcon_verify(sig, siglen, FALCON_SIG_PADDED,
			pk, pklen, msg, msglen, tmp, tmplen) != 0)
		{
			fprintf(stderr, "verify failed\n");
			return 1;
		}
		t[i] = now_ms() - t0;
	}
	report("verify", logn, t, nvf);

	return 0;
}
