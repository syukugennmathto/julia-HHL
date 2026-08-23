# Deep Research prompt — Is any real FALCON / FN-DSA deployment built with
# floating-point contraction (`-ffast-math` / `-Ofast` / `-ffp-contract=fast`)?

Paste everything below the line into a Deep Research tool. It is self-contained.

---

## The question

I need to know whether any **real, shipped or widely-used build of FALCON
(a.k.a. FN-DSA, the NIST post-quantum signature being standardized as FIPS 206)**
compiles its signing code with floating-point *contraction* enabled — i.e. with
any of `-ffast-math`, `-Ofast`, `-ffp-contract=fast`, or a global flag that
implies them — on a target where the FMA (fused multiply-add) instruction is
available.

This is a security question, not a performance one. FALCON's signature is a
deterministic function of double-precision floating-point arithmetic in its
Gaussian sampler. If a compiler fuses `a*b + c` into a single rounding (an FMA)
in that arithmetic, the produced signatures can differ, at a rate around
1 in 10,000, from a build that does not — and a pair of differing signatures on
the same message leaks the private key. So I need concrete evidence of *which
builds, if any, enable contraction*, with sources.

## Essential background so you can judge relevance (verify, don't assume)

1. **The reference implementation wraps `double` in a one-field struct**
   (`typedef struct { double v; } fpr;` in `fpr.h`), and does complex
   multiplication through separate inline helpers (`FPC_MUL`). Empirically, this
   wrapper *blocks* contraction for GCC at every optimization level, and for
   clang except at `-ffp-contract=fast` / `-Ofast` / `-ffast-math`. So a build of
   the **unmodified C reference** is exposed **only** when compiled with clang
   (or a clang-derived toolchain) at one of those settings. Confirm whether this
   wrapper is still present in current sources.
2. **The reference’s `config.h` forces the integer-only emulated FP path
   (`FALCON_FPEMU`) on by default**, "because native FPUs may yield slight
   discrepancies that could affect determinism." The emulated path cannot
   contract (no floating-point at all). So the exposure exists only where an
   integrator **turns on the native-`double` path (`FALCON_FPNATIVE`)** *and*
   uses a contracting compiler. Check whether each project below uses the
   emulated or the native path, and whether that default has changed.
3. **A re-implementation that does NOT use the struct wrapper** (a plain-`double`
   C port, or a port to another language that emits FP directly) loses the
   accidental protection and can contract under ordinary flags — including GCC’s,
   and including on aarch64 where FMA is base-ISA. So ports matter as much as the
   reference.
4. Therefore the risk checklist for any project is: **(native FP path?) AND
   (contracting compiler/flags?) AND (wrapper absent or defeated?)**.

## What to investigate — concrete targets

For each, find the *canonical build configuration actually used* (CMake/Make/
autotools files, packaging recipes, CI YAML, release build scripts — not test or
micro-benchmark targets, which often set `-Ofast` harmlessly on throwaway code):

- **liboqs** (Open Quantum Safe) — how it builds Falcon: default `CMAKE_C_FLAGS`,
  `-O` level, whether it passes `-march=native`/`-mfma`, whether Falcon uses the
  AVX2 or the reference path, and whether any build preset enables `-ffast-math`
  or `-Ofast`. Also its downstream: **oqs-provider** (OpenSSL), **oqs-boringssl**,
  and language wrappers.
- **PQClean** — the Falcon `clean` and `avx2` implementations and their Makefiles;
  what flags a consumer inherits.
- **The NIST/PQC-forum reference submission package** and the maintained mirror(s)
  (e.g. `algorand/falcon`) — current default FP path and flags.
- **OpenSSL 3.5+**, **AWS-LC / BoringSSL**, **wolfSSL/wolfCrypt**, **Botan**,
  **mbedTLS**, **Cloudflare CIRCL**, **Bouncy Castle** — for each: do they ship
  FN-DSA/Falcon at all yet (many ship only ML-DSA/ML-KEM)? If yes, native or
  emulated FP, and what build flags. If FN-DSA is absent, say so — that project
  is N/A for this question.
- **Go standard library `crypto`** and the Go ecosystem — confirm whether FN-DSA
  is present (as of 2026 the stdlib shipped ML-DSA/ML-KEM; verify FN-DSA status).
- **Rust `pqcrypto` / `pqclean` crates**, and any pure-Rust Falcon — build
  profiles, `codegen-units`, `fast-math`-like settings (`-C llvm-args`,
  `fadd_fast`), and whether they bind the C reference or reimplement.
- **Distribution packaging**: Debian/Ubuntu, Fedora/RHEL, Arch, Alpine,
  Homebrew, vcpkg, Nixpkgs — for any package that contains Falcon (often via
  liboqs). Do the distro-wide hardening/optimization flags include `-ffast-math`
  or `-Ofast`? (Standard distro flags are usually `-O2` without fast-math, but
  confirm; and check for any package that overrides.)
- **Embedded / constant-time forks**: `pqm4` (Cortex-M4), `pqriscv`, and any IoT
  SDK shipping Falcon — these typically use the emulated/assembly path, but
  confirm, and check whether any uses `-Ofast`.
- **Large applications that vendor a PQC library into one big build** with a
  global `-Ofast` / `-ffast-math` (this is the dangerous pattern: a whole
  application built with `-ffast-math` would sweep Falcon’s translation units
  into contraction). Look for messaging/VPN/TLS products advertising Falcon/
  FN-DSA and any public build flags.

## What counts as evidence

- A **link to the actual build file** (CMakeLists.txt, Makefile, `.mk`, spec/
  rules file, CI config) showing the flags, with the specific lines quoted.
- The **compiler and target** the canonical build uses (GCC vs clang; `-march`/
  `-mfma`; x86-64 vs aarch64).
- Whether **`FALCON_FPNATIVE` or `FALCON_FPEMU`** is selected, or whether the
  port uses native `double` unconditionally.
- Whether the build sets **`-ffp-contract=off`** or **`#pragma STDC FP_CONTRACT
  OFF`** anywhere near the Falcon sources (that would be a deliberate
  mitigation — note it).
- Distinguish **library build flags** (what matters) from **test/benchmark
  flags** (usually irrelevant) and from **flags a downstream consumer could add**.

## Deliverable

A table with one row per project:

| Project | FN-DSA/Falcon shipped? | FP path (native / emulated / port) | Canonical compiler + flags | Contraction on? (`-ffast-math`/`-Ofast`/`-ffp-contract=fast`, FMA target) | `FP_CONTRACT OFF` set? | Verdict (EXPOSED / protected / N/A) | Source links |

Then a short synthesis answering directly:
1. Does **any** real, non-test build of Falcon/FN-DSA enable contraction on an
   FMA-capable target? Name it, or state that none was found.
2. Among projects that ship it, how many use the **native FP path** (the
   necessary precondition), and how many rely on the **emulated path** (immune)?
3. Are there projects using a **wrapper-free port** that would contract under
   ordinary GCC/clang flags?
4. Note any project that already sets `-ffp-contract=off` / `FP_CONTRACT OFF`
   deliberately.

Be skeptical and cite everything. If evidence is ambiguous (e.g. flags depend on
a preset), say which preset and how commonly it’s used. Prefer primary sources
(repository files, package recipes) over blog posts.
