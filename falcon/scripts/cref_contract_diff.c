/* cref_contract_diff.c -- does the FALCON reference implementation, compiled
 * two conforming ways, produce the same signatures?
 *
 * C99 6.5p8 lets an implementation contract `a*b + c` into a single fused
 * multiply-add with one rounding instead of two; FP_CONTRACT governs it.  The
 * reference's FPC_MUL expands (after inlining the native-double `fpr_mul` /
 * `fpr_sub`) to exactly that shape:
 *
 *     d_re = a_re*b_re - a_im*b_im;
 *     d_im = a_re*b_im + a_im*b_re;
 *
 * GCC's default is -ffp-contract=fast, and it contracts across statements.
 * clang's default for standard C is off.  The reference sets neither, and the
 * FALCON specification does not mention FP_CONTRACT at all.
 *
 * This program signs N messages with one fixed key, on a per-message
 * deterministic PRNG tape, and writes one 8-byte FNV-1a hash per signature to
 * stdout.  Build it twice with different contraction settings and `cmp` the two
 * files: the first differing offset names the message index.  Everything else
 * -- key, salt, sampler tape -- is identical by construction, so a difference
 * can only come from the floating point.
 *
 * Note that contraction can only happen if the target ISA has an FMA
 * instruction.  On baseline x86-64 it does not, so a distro build does not
 * contract however it is compiled; with -march=native on any post-Haswell x86,
 * and on *every* aarch64 target (FMA is in the base ISA), it does.
 *
 *     gcc -O2 -march=native -ffp-contract=fast -o cd_fast ...
 *     gcc -O2 -march=native -ffp-contract=off  -o cd_off  ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "falcon.h"

static unsigned long long fnv(const unsigned char *p, size_t n) {
    unsigned long long h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

/* -d i1,i2,... : instead of the hash stream, dump the key, and the message and
 * signature of the listed indices, in hex, so the two builds' divergent pairs
 * can be handed to the key-recovery of ePrint 2024/1709 section 5.1. */
static int wanted(const char *list, long j) {
    if (!list) return 0;
    const char *p = list;
    while (*p) {
        long v = strtol(p, (char **)&p, 10);
        if (v == j) return 1;
        if (*p == ',') p++; else if (*p) return 0;
    }
    return 0;
}

static void hexline(const char *tag, long j, const unsigned char *p, size_t n) {
    printf("%s %ld ", tag, j);
    for (size_t i = 0; i < n; i++) printf("%02x", p[i]);
    printf("\n");
}

int main(int argc, char **argv) {
    unsigned logn = 9;
    long n = argc > 1 ? atol(argv[1]) : 100000;
    const char *dump = (argc > 3 && strcmp(argv[2], "-d") == 0) ? argv[3] : 0;

    shake256_context rng;
    unsigned char seed[32];
    for (int i = 0; i < 32; i++) seed[i] = (unsigned char)i;
    shake256_init_prng_from_seed(&rng, seed, sizeof seed);

    size_t sklen = FALCON_PRIVKEY_SIZE(logn), pklen = FALCON_PUBKEY_SIZE(logn);
    unsigned char *sk = malloc(sklen), *pk = malloc(pklen);
    size_t tklen = FALCON_TMPSIZE_KEYGEN(logn);
    unsigned char *tk = malloc(tklen);
    if (falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen, tk, tklen)) {
        fprintf(stderr, "keygen failed\n"); return 1;
    }
    /* The key itself is a floating-point-dependent object: keygen runs the FFT
     * for its rejection test.  Report its hash on stderr so the two builds can
     * be checked to agree on the key before their signatures are compared. */
    fprintf(stderr, "privkey fnv %016llx\npubkey  fnv %016llx\n",
            fnv(sk, sklen), fnv(pk, pklen));
    if (dump) { hexline("privkey", 0, sk, sklen); hexline("pubkey", 0, pk, pklen); }

    size_t tslen = FALCON_TMPSIZE_SIGNDYN(logn);
    unsigned char *ts = malloc(tslen);
    size_t cap = FALCON_SIG_PADDED_SIZE(logn);
    unsigned char *sig = malloc(cap);

    for (long j = 0; j < n; j++) {
        char msg[64];
        int mlen = snprintf(msg, sizeof msg, "contract/%ld", j);

        /* a fresh, deterministic tape per message: identical in both builds */
        shake256_context srng;
        unsigned char s2[32];
        memset(s2, 0, sizeof s2);
        for (int i = 0; i < 8; i++) s2[i] = (unsigned char)((j >> (8*i)) & 0xff);
        shake256_init_prng_from_seed(&srng, s2, sizeof s2);

        size_t siglen = cap;
        if (falcon_sign_dyn(&srng, sig, &siglen, FALCON_SIG_PADDED,
                            sk, sklen, msg, (size_t)mlen, ts, tslen)) {
            fprintf(stderr, "sign failed at %ld\n", j); return 1;
        }
        if (dump) {
            if (wanted(dump, j)) {
                hexline("message", j, (const unsigned char *)msg, (size_t)mlen);
                hexline("signature", j, sig, siglen);
            }
            continue;
        }
        unsigned long long h = fnv(sig, siglen);
        unsigned char out[8];
        for (int i = 0; i < 8; i++) out[i] = (unsigned char)((h >> (8*i)) & 0xff);
        fwrite(out, 1, 8, stdout);
    }
    fflush(stdout);
    return 0;
}
