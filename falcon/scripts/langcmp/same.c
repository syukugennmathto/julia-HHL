/* The same algorithm as Falcon.polymulq, written in C, to separate "Julia is
 * slower than C" from "our algorithm is slower than C's algorithm".
 *
 * Deliberately a transliteration, not an optimisation: same loop order, same
 * int64 accumulator, same negacyclic sign, same final mod.  No restrict, no
 * intrinsics, no unrolling by hand -- whatever -O2 does to it, -O2 would do to
 * the equivalent Julia too. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define Q 12289

static void polymulq(const long *f, const long *g, long *out, int n)
{
	long *acc = calloc((size_t)n, sizeof *acc);
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			int k = i + j;
			long p = f[i] * g[j];
			if (k < n) {
				acc[k] += p;
			} else {
				acc[k - n] -= p;
			}
		}
	}
	for (i = 0; i < n; i++) {
		long c = acc[i] % Q;
		out[i] = c < 0 ? c + Q : c;      /* Julia's mod, not C's % */
	}
	free(acc);
}

static double now_ms(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

static int cmpd(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

int main(int argc, char **argv)
{
	int n = argc > 1 ? atoi(argv[1]) : 512;
	int iters = argc > 2 ? atoi(argv[2]) : 300;
	long *f = malloc((size_t)n * sizeof *f);
	long *g = malloc((size_t)n * sizeof *g);
	long *o = malloc((size_t)n * sizeof *o);
	double *t = malloc((size_t)iters * sizeof *t);
	int i;
	unsigned s = 12345;

	for (i = 0; i < n; i++) {                    /* the same LCG as the Julia side */
		s = s * 1103515245u + 12345u;
		f[i] = (long)((s >> 8) % Q);
		s = s * 1103515245u + 12345u;
		g[i] = (long)((s >> 8) % Q);
	}
	polymulq(f, g, o, n);
	for (i = 0; i < iters; i++) {
		double t0 = now_ms();
		polymulq(f, g, o, n);
		t[i] = now_ms() - t0;
	}
	qsort(t, (size_t)iters, sizeof *t, cmpd);
	printf("C   polymulq n=%d: median %.4f ms   (checksum %ld)\n",
		n, t[iters / 2], o[0] + o[n - 1]);
	return 0;
}
