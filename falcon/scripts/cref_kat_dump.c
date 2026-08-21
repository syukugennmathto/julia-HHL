/* Emit a FALCON-512 key + signature from a fixed seed, for cross-checking. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "falcon.h"

static void hex(const char *name, const unsigned char *p, size_t n) {
    printf("%s ", name);
    for (size_t i = 0; i < n; i++) printf("%02x", p[i]);
    printf("\n");
}

int main(void) {
    unsigned logn = 9;
    shake256_context rng;
    unsigned char seed[32];
    for (int i = 0; i < 32; i++) seed[i] = (unsigned char)i;
    shake256_init_prng_from_seed(&rng, seed, sizeof seed);

    size_t sklen = FALCON_PRIVKEY_SIZE(logn), pklen = FALCON_PUBKEY_SIZE(logn);
    unsigned char *sk = malloc(sklen), *pk = malloc(pklen);
    size_t tklen = FALCON_TMPSIZE_KEYGEN(logn);
    unsigned char *tk = malloc(tklen);
    int r = falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen, tk, tklen);
    if (r) { fprintf(stderr, "keygen failed %d\n", r); return 1; }
    printf("logn %u\n", logn);
    hex("privkey", sk, sklen);
    hex("pubkey", pk, pklen);

    const char *msg = "falcon-jl cross-check message";
    size_t msglen = strlen(msg);
    hex("message", (const unsigned char *)msg, msglen);

    size_t tslen = FALCON_TMPSIZE_SIGNDYN(logn);
    unsigned char *ts = malloc(tslen);
    size_t siglen = FALCON_SIG_PADDED_SIZE(logn);
    unsigned char *sig = malloc(siglen);
    r = falcon_sign_dyn(&rng, sig, &siglen, FALCON_SIG_PADDED,
                        sk, sklen, msg, msglen, ts, tslen);
    if (r) { fprintf(stderr, "sign failed %d\n", r); return 1; }
    hex("signature", sig, siglen);

    size_t tvlen = FALCON_TMPSIZE_VERIFY(logn);
    unsigned char *tv = malloc(tvlen);
    r = falcon_verify(sig, siglen, FALCON_SIG_PADDED, pk, pklen, msg, msglen, tv, tvlen);
    printf("selfverify %d (0 means OK)\n", r);
    return 0;
}
