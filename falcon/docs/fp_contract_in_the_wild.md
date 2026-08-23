# Is floating-point contraction enabled in real FALCON / FN-DSA builds?

A Deep Research survey (2026-08) of shipped and widely-used FALCON/FN-DSA builds,
with the reference-source claims re-verified here against the vendored copy in
`scripts/cref/`.

## Answer

**No currently-shipped, non-test build enables contraction on an FMA-capable
target.** Every major distribution route is protected, by one of two means:

1. **Integer emulation** (`FALCON_FPEMU`): PQClean's `clean` implementation
   (`typedef uint64_t fpr;` — no `double` at all), pqm4, and the reference's
   *default*. No floating point, so nothing to contract.
2. **Explicit AVX2 intrinsics** (`FALCON_AVX2` without `FALCON_FMA`): the FFT's
   complex products are written as `_mm256_mul_pd` + `_mm256_add_pd` (two
   roundings), which a compiler does not contract. This is liboqs's and
   PQClean's native path, chosen by runtime CPU detection.

liboqs (and its downstream: oqs-provider, wolfSSL, language bindings, distro
packages), PQClean, Rust `pqcrypto-falcon` (binds PQClean `clean`), Bouncy Castle
(pure-Java port; Java does not fuse `a*b+c` without explicit `Math.fma`), and the
reference itself all fall here. OpenSSL 3.5, AWS-LC/BoringSSL, Botan, mbedTLS,
Cloudflare CIRCL and the Go standard library **do not ship FN-DSA/Falcon at all**
as of 2026 (they ship ML-DSA / ML-KEM / SLH-DSA), so they are N/A.

**So the threat this paper measures is latent, not active** — which is the right
time to fix a standard, before native-FP builds proliferate under FIPS 206.

## Why "latent" is not "harmless" — the reference documents the exact danger

The protection does **not** come from the specification. It comes from three
things the *reference* does, none of them normative:

1. **`config.h` forces `FALCON_FPEMU` on by default** and carries a
   **CRITICAL SECURITY WARNING** that native FP and FMA should be *disabled* for
   determinism, "to prevent a potential catastrophic security failure in the
   deterministic mode."
2. The AVX2 path uses hand-written two-rounding intrinsics unless `FALCON_FMA`
   is set.
3. The `fpr` struct wrapper (`typedef struct { double v; } fpr;`) blocks
   compiler contraction of the scalar native path except under clang
   `-ffp-contract=fast` (§6.7).

The reference even ships an **explicit `FALCON_FMA` build option**, and its own
`config.h` says of it, verbatim (lines 124–134):

> "setting this option will slightly modify the values of expanded private keys,
> but will normally not change the values of non-expanded private keys, public
> keys or signatures, for a given keygen/sign seed (non-expanded private keys and
> signatures might theoretically change, but only with low probability, less than
> 2^(-40); produced signatures are still safe and interoperable)."

This is the paper's thesis stated by the reference itself — and then dismissed.
"Produced signatures are still safe and interoperable" is the claim we contest:
in the deterministic / cross-implementation setting a single changed signature,
observed alongside its unchanged counterpart, is a full key-recovery event
(ePrint 2024/1709 §5.1, and §6.7 here). §6.7.x measures `FALCON_FMA`'s actual
signature-change rate against that `2^(-40)` estimate.

## The residual exposure

A re-implementation that (a) uses the native `double` path, (b) omits the `fpr`
wrapper (any port to a language that emits FP directly, or a plain-`double` C
port), and (c) is built with contraction — including plain GCC on aarch64, where
FMA is base-ISA — loses all three protections at once. No such build was found
shipping today, but nothing in the specification prevents one, and FIPS 206's
move toward native-FP, bit-exact implementations is exactly the pressure that
would produce one. `rust-fn-dsa` (Pornin's pure-Rust FN-DSA reimplementation) is
the closest current example of a wrapper-free port; its FP handling was not
verified, and it showed no evidence of contraction, but it is the shape to watch.

## Housekeeping note

PQClean was archived (read-only) on 2026-08-04, and RustSec flagged
`pqcrypto-falcon` as unmaintained on 2026-06-04. Future Falcon C maintenance
moves to successors (e.g. the PQ Code Package); whether they keep the emulated
default and the intrinsic AVX2 path is the thing to track.

## Per-project verdicts

| Project | Ships Falcon? | FP path | Contraction? | Verdict |
|:--|:--|:--|:--|:--|
| reference (algorand/falcon) | yes | FPEMU default; wrapper on native | no (`-O3`, no fast-math) | protected |
| PQClean `clean` | yes | integer (`uint64_t fpr`) | n/a | protected |
| PQClean `avx2` | yes | native double + explicit intrinsics | no (`FALCON_FMA` off) | protected |
| liboqs (+ oqs-provider, wolfSSL, distros) | yes | clean / avx2 by CPU detect | no | protected |
| Rust `pqcrypto-falcon` | yes | binds PQClean clean | no | protected |
| Bouncy Castle (Java/C#) | yes | Java double (no implicit fuse) | no | protected |
| pqm4 | yes | FPEMU + M4 asm | no | protected |
| OpenSSL 3.5 / AWS-LC / Botan / mbedTLS / CIRCL / Go | **no** | — | — | N/A |
| a wrapper-free native-double port + contraction | (hypothetical) | native double | yes | **EXPOSED** |

Source: user-run Deep Research, 2026-08; reference-source quotes re-verified
against `scripts/cref/config.h` and `scripts/cref/inner.h`.
