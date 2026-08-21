/*
 * cref_fft_dump.c -- dump FFT vectors from the C reference implementation.
 *
 * Used to generate test/vectors/fft_c_kat.jl, so that the Julia test suite can
 * be checked against the *C* reference and not only against the Python one.
 * Committed rather than run automatically: it needs the C sources, which are
 * not vendored here (see docs/build_cref_macos.md for how to get them).
 *
 * Build (from inside a checkout of https://github.com/algorand/falcon, with
 * this file copied in):
 *
 *     clang -O2 -o cref_fft_dump cref_fft_dump.c \
 *         codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c \
 *         sign.c vrfy.c
 *     ./cref_fft_dump 3        # logn = 3, i.e. n = 8
 *
 * Output format, one section per invocation:
 *
 *     logn <logn>
 *     fpemu <0|1>
 *     in <n space-separated integers>
 *     out <n space-separated doubles, %.17g>
 *
 * The `out` line is the C FFT representation: n/2 real parts followed by n/2
 * imaginary parts.  See docs/math/05_fft.md for the mapping onto our (and the
 * Python reference's) n-element complex representation -- it is *not* the
 * identity: C index k corresponds to our index gray(k) = k XOR (k >> 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inner.h"

/*
 * In FPEMU mode `fpr` is a uint64_t holding IEEE-754 binary64 bit patterns.
 * In FPNATIVE mode it is a struct wrapping a double.  Note that the reference
 * function named fpr_double() means "multiply by two", NOT "convert to
 * double" -- an easy and costly misreading.
 */
static double
as_double(fpr x)
{
	double d;
#if FALCON_FPEMU
	uint64_t u = x;
	memcpy(&d, &u, sizeof d);
#else
	d = x.v;
#endif
	return d;
}

int
main(int argc, char **argv)
{
	unsigned logn;
	size_t n, i;
	fpr *f;
	int *in;

	if (argc != 2) {
		fprintf(stderr, "usage: %s <logn>\n", argv[0]);
		return 1;
	}
	logn = (unsigned)atoi(argv[1]);
	n = (size_t)1 << logn;
	f = malloc(n * sizeof *f);
	in = malloc(n * sizeof *in);

	/*
	 * Deterministic input, small enough to stay exact in a double and
	 * signed so that the negacyclic sign errors show up.
	 */
	for (i = 0; i < n; i++) {
		in[i] = (int)((i * 37 + 11) % 23) - 11;
		f[i] = fpr_of((int64_t)in[i]);
	}

	printf("logn %u\n", logn);
	printf("fpemu %d\n", FALCON_FPEMU);
	printf("in");
	for (i = 0; i < n; i++) {
		printf(" %d", in[i]);
	}
	printf("\n");

	Zf(FFT)(f, logn);

	printf("out");
	for (i = 0; i < n; i++) {
		printf(" %.17g", as_double(f[i]));
	}
	printf("\n");
	return 0;
}
