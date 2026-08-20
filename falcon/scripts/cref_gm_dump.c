/*
 * Dump the C reference implementation's FFT root table, fpr_gm_tab.
 *
 * Two things are being established at once:
 *
 *   1. that the emulated-floating-point build and the native-floating-point
 *      build carry the *same* table -- FPEMU stores raw uint64 bit patterns,
 *      FPNATIVE stores 27-digit decimal literals for the compiler to round,
 *      and it is not obvious a priori that those agree in the last bit;
 *
 *   2. what that table actually is, so falcon/src/fft.jl can be compared
 *      against it entry by entry rather than "to within a tolerance".
 *
 * Build both ways and diff the outputs:
 *
 *     cc -O2                                   -o gmdump_emu cref_gm_dump.c fpr.c
 *     cc -O2 -DFALCON_FPEMU=0 -DFALCON_FPNATIVE=1 -o gmdump_nat cref_gm_dump.c fpr.c
 *     ./gmdump_emu > gm_emu.txt && ./gmdump_nat > gm_nat.txt && diff gm_emu.txt gm_nat.txt
 *
 * (config.h in the reference guards both macros with #ifndef, so -D on the
 * command line wins; see docs/build_cref_macos.md.)
 *
 * The table holds 1024 complex numbers as 2048 consecutive fpr, real part
 * first: entry j is (fpr_gm_tab[2*j], fpr_gm_tab[2*j+1]).  Slot 0 is unused
 * and reads as zero.  Reading past j = 1023 walks off the end -- the table is
 * sized for logn <= 10, and there is no sentinel.
 *
 * `fpr` is uint64_t under FPEMU and struct { double v; } under FPNATIVE, and
 * the FPEMU encoding *is* IEEE-754 binary64.  So in both builds the eight
 * bytes of an fpr are the bit pattern of the double it denotes, and memcpy
 * out of it is meaningful without a per-build conversion.
 */
#include <stdio.h>
#include <string.h>
#include "inner.h"

int
main(void)
{
	int j;

	for (j = 0; j < 1024; j++) {
		unsigned long long re, im;

		memcpy(&re, &fpr_gm_tab[2 * j], sizeof re);
		memcpy(&im, &fpr_gm_tab[2 * j + 1], sizeof im);
		printf("%d %016llx %016llx\n", j, re, im);
	}
	return 0;
}
