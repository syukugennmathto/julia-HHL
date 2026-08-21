/* Instrumentation for verifying three claims about the C reference's Babai
 * reduction (docs/debug_log.md #041):
 *   H1 the reduction works on binary base-2^31 limbs, not RNS
 *   H2 it strips 25 bits per iteration
 *   H3 the top-53-bit extraction is a pointer offset, not a shift
 */
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

static int iters[12];
static int dumped[12];

void falcon_instr_enter(unsigned logn, unsigned depth, int maxbl_FG,
	int minbl_fg, int maxbl_fg, size_t slen)
{
	printf("# enter depth=%u n=%u  maxbl_FG=%d minbl_fg=%d maxbl_fg=%d slen=%zu\n",
		depth, 1u << logn, maxbl_FG, minbl_fg, maxbl_fg, slen);
	iters[depth] = 0;
}

void falcon_instr_iter(unsigned depth, int scale_k)
{
	iters[depth]++;
	if (iters[depth] <= 3 || scale_k <= 0) {
		printf("#   depth=%u iter=%d scale_k=%d\n", depth, iters[depth], scale_k);
	}
	if (scale_k <= 0) {
		printf("# depth=%u TOTAL ITERATIONS = %d\n", depth, iters[depth]);
	}
}

/* Print the raw words of coefficient 0 of Ft, once per depth.  If the
 * representation is positional base-2^31 the Julia side can reconstruct the
 * integer from them; if it were RNS it could not. */
void falcon_instr_words(unsigned depth, const uint32_t *Ft, size_t FGlen)
{
	size_t v;

	if (dumped[depth]) {
		return;
	}
	dumped[depth] = 1;
	printf("WORDS depth=%u len=%zu", depth, FGlen);
	for (v = 0; v < FGlen; v++) {
		printf(" %u", Ft[v]);
	}
	printf("\n");
}
