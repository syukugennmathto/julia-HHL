# Two faithful implementations, two signatures: specification underdetermination in FALCON / FN-DSA

*Working manuscript. This document is generated from the artifact in this
repository; every number in it is produced by a script named in the text, and
the scripts are runnable without a C compiler except where noted.*

*Status: all measurements reported here are final. Limitations that remain
open are stated as such in §12.*

---

## Abstract

FALCON — standardized by NIST as FN-DSA — is the only selected signature whose
output is a function of floating-point arithmetic, and its specification's
precision analysis is about the sampled *distribution*, not about whether two
conforming implementations produce the same bytes. We write an independent
implementation in Julia, from the specification, and reconcile it with the C
reference until it reproduces the reference's signatures **byte for byte**. That
reconciliation exposes three classes of place where the specification does not
determine the intermediate values: how an operation is spelled, a constant that
cannot be derived from the published parameters, and a smoothing parameter that
is used but not printed. We then measure, over hundreds of thousands of
signatures with matched controls, which of these reach the signature. Two
spellings — complex division and the `D11` entry of LDL\* — change the signature
at ≈ 7.5 × 10⁻⁶ per signature, always at the last two sampler calls, the
positions where ePrint 2024/1709 (Lin–Tibouchi–Yu–Zhang, EUROCRYPT 2025) turns a
single discrepant pair into full key recovery — and we run that recovery,
reconstructing the private key from one A2 discrepant pair and the public key
alone. Unlike the perturbations that paper studies, these are differences
between the *specification* and the *reference*, present in both of the
reference's signing modes, and not removed by its `sign_dyn`/`sign_tree`
countermeasure. We reproduce that paper's mechanism independently, check its
Heuristic 1 directly (which it does not), confirm the effect at n = 1024, and
offer a **mechanism for the asymmetry it left open** — why the last two sampler
calls diverge far more readily than the first two — as a perturbation that
accumulates along the tree descent and is 32× larger at the last two calls
(§12.3 states what that mechanism does and does not settle). We also
implement and evaluate that paper's own proposed countermeasure. It is a
two-part fix — round instead of floor, *and* force `‖(g,−f)‖²` odd — and we find
that **part 1 alone, the change an implementer makes if they touch only the
sampler, gives no protection where it matters**: `q = 12289` is odd, so rounding
immunises the first two sampler calls completely, while `‖(g,−f)‖²` is even on
every reference key, so it immunises the last two — the full-key-recovery
positions — not at all. With part 1 deployed we still recover keys, 3 of 3
divergent pairs in 160000 signatures. Part 2 *is* deployable (we show a
one-line key-generation change yields odd norms, 40 of 40), so the two parts
must be mandated together and a standard must say so. We also give a simpler
alternative that needs no key-generation change: snapping the centre to its
exact rational value, which the signer can do because it knows the denominator. Performance: this
implementation is 20–300× faster than the Python reference, beats the widely
deployed emulated-floating-point C build on key generation and signing at
n = 512, and is the fastest of every build we measured at verification, while
remaining ≈ 6–7× slower than native-floating-point C at signing. Constant-time
behaviour is explicitly out of scope, for a reason we argue is methodological
rather than a concession.

## 要旨 (Japanese abstract)

FALCON（NIST が FN-DSA として標準化する格子署名）の署名は、仕様書が
**決めていない**書き方の選択に依存する。本稿は、仕様書だけを読んで独立に
書いた Julia 実装が参照実装の署名を**バイト単位で再現**することを示したうえで、
仕様書が中間値を固定していない箇所を 3 クラス同定する：**演算の書き方**、
**導出できない定数**、**印字されていない ε**。そのうち複素除算と `LDL*` の
`D11` の書き方は、参照実装の 2 つの signing mode のどちらにも共通に存在し、
40 万署名あたり 3 本の割合で署名を食い違わせ、その食い違いは ePrint 2024/1709
（Lin–Tibouchi–Yu–Zhang, EUROCRYPT 2025）が秘密鍵の全復元に使う位置
（呼び出し 2n−1, 2n）に集中する。この差は同論文 §7.2 の対策では塞がらない。
機構は同論文 Lemma 1（`floor(mu)` の不連続）そのものであり、本稿はそれを
**再発見ではなく別経路からの追認**として位置づけ、同論文の Heuristic 1 を
（同論文が測っていない側から）直接検証する。性能面では、この実装は Python
参照実装より 20〜270 倍速く、広く配備されている整数エミュレーション版 C を
鍵生成・署名で上回り、**検証ではすべての C ビルドより速い**が、ネイティブ
浮動小数点版 C には署名で約 6 倍、鍵生成で約 2 倍遅い。定数時間性は明示的に
スコープ外とし、その選択自体が「処理系由来の変時間性とアルゴリズム由来の
それを分離する」という本研究の方法論の一部であることを述べる。

---

## 1. Introduction

FALCON is the compact lattice signature NIST selected for standardization as
FN-DSA. Alone among the selected signatures it computes its output through
**floating-point arithmetic**: its Gaussian sampler runs a fast-Fourier
nearest-plane recursion over `Complex{Float64}`, and the sampled integers, and
therefore the signature, are a function of the exact floating-point values that
recursion produces. ML-DSA (FIPS 204) and SLH-DSA (FIPS 205) are integer- and
hash-based; given the same key, message and randomness they return the same
bytes on any conforming implementation, by construction. FALCON carries no such
guarantee, and its own specification is explicit that it does not intend to:
§2.5.2 analyzes the required floating-point precision through a Rényi-divergence
argument about the *output distribution* (concluding that 53 bits lose no
security because the sampled distribution is close enough), and §4.1 notes that
C's `double` is only required to "match at least" IEEE-754 and that "FALCON
works properly with such limited floating-point types." Both statements are
about the distribution the sampler induces, not about whether two conforming
implementations produce the same signature on the same input.

That distinction is now load-bearing. The FIPS 206 status update (Perlner,
NIST, September 2025) indicates the standard will fix the order of
floating-point operations, forbid fused multiply-add, and require
implementations to match known-answer tests bit for bit. **If the standard
requires bit-exact agreement, then the question "does the specification
determine the intermediate values?" becomes a normative one.** This paper
answers it empirically, by building a second implementation from the
specification and measuring exactly where it disagrees with the reference.

### Why Julia, honestly

We chose Julia because we expected it to be fast — specifically because its
native arrays and its reputation for numerical and linear-algebra performance
suggested the matrix-heavy parts of FALCON (the FFT, the NTT, the Gram matrix)
would be quick to write and quick to run. That expectation was only partly
borne out (§9): the implementation beats the Python reference everywhere and
the emulated C build on two of three operations and every C build on
verification, but it loses to native-floating-point C on signing and key
generation.

There is, however, a second reason that is genuine but was recognized only
after the fact, and we state it as such rather than dressing it up as the
original motivation. FALCON's security-relevant timing variability has two
sources: **algorithmic** (the rejection sampler's data-dependent loop count)
and **implementational** (whatever the language and its runtime do — bounds
checks, garbage collection, JIT compilation, dynamic dispatch). A constant-time
C implementation works hard to remove *both*. Writing FALCON in a
garbage-collected, JIT-compiled language makes the second source impossible to
close — and that makes it a clean instrument for studying the *first* in
isolation, and for studying a third thing that neither a hardened C
implementation nor the specification's distributional analysis makes visible:
the places where two faithful implementations, in different languages, disagree
because the specification never said which of two equal formulas to use. A
language that cannot hide implementation-level variance is exactly the language
in which specification-level variance shows up undisguised. We did not pick
Julia for that reason, but it is why the choice turned out to be the right one
for this particular question.

### Contributions

1. **A byte-exact independent reproduction** of the FALCON-512 reference
   signatures, written from the specification and the Python reference and then
   reconciled with the C reference down to the last bit of `s2` and of the
   compressed signature (§4). This is the instrument the rest of the paper
   uses.
2. **Three classes of specification underdetermination** (§5): how an operation
   is spelled, a constant that cannot be derived from the published parameters,
   and a smoothing parameter that is used but not printed. For the first class
   we exhibit both implementation patterns — the reference's spelling and the
   specification-reader's — as running code.
3. **A measurement of which of these reach the signature, and how often** (§6),
   with matched positive and negative controls, and a direct check of the
   arithmetic mechanism ePrint 2024/1709 identifies.
4. **One instance that ePrint 2024/1709 does not cover** (§6.1, §7): complex
   division and the `D11` entry of LDL\* differ between the specification and the
   reference, are present in **both** of the reference's signing modes, reach
   the signature at ≈ 7.5 × 10⁻⁶ per signature at the key-recovery position, are
   not removed by that paper's `sign_dyn`/`sign_tree` countermeasure, and
   **from one such pair we recover the private key** (§6.1). The effect
   replicates at n = 1024.
5. **An honest three-way performance comparison** (§9) against a range of C
   builds (two compilers, three optimization levels, emulated and native
   floating point) and the Python reference, reported as distributions rather
   than single numbers.
6. **A cross-scheme scoping** (§11.1): from the reference code of Mitaka,
   Antrag and HAWK, the sensitivity requires both floating point in the signing
   sampler and a small-denominator rational centre — Falcon alone has both, so
   the hazard is specific to its design, not generic to lattice hash-and-sign.
7. **An evaluation of the proposed countermeasure** (§7.1): we implement ePrint
   2024/1709's Algorithm 4, verify the rational structure of the centres
   directly, and show that its *first part alone* immunises only the harmless
   sampler positions — recovering keys, 3 of 3, with that part deployed. We also
   show its second part is deployable (contrary to what the reference key
   generator suggests), so the normative point is that the two must be mandated
   **together**. §7.2 gives a simpler alternative requiring no key-generation
   change.

We are equally explicit about what is **not** a contribution. The arithmetic
mechanism (an integer centre passing through `floor`), the fact that the
reference's two signing modes diverge, and the insensitivity of the sampler to
its standard deviation are all results of ePrint 2024/1709; §7 states precisely
which of our findings reproduce theirs and which extend them. Constant-time
behaviour is out of scope (§10).

---

## 2. Background

### 2.1 FALCON signing

A FALCON private key is an NTRU trapdoor basis `B = [[g, -f], [G, -F]]` with
`fG - gF = q` over `Z[x]/(x^n + 1)`, `q = 12289`, `n ∈ {512, 1024}`. The public
key is `h = g/f mod q`. To sign a message, FALCON hashes it (with a 320-bit
salt) to a point `c`, computes the target `t = (FFT(c), 0) · FFT(B)^{-1}`, and
samples a lattice point near `t` with a fast-Fourier Gaussian sampler
`ffSampling`, obtaining a short `s = (s1, s2)` with `s1 + s2 h ≡ c (mod q)` and
`‖s‖² ≤ ⌊β²⌋`. The signature transmits the salt and a compressed `s2`.

`ffSampling` walks a binary tree (the LDL\* decomposition of the Gram matrix
`BB*`, in the FFT domain), and at each of its `2n` leaves draws one integer from
a discrete Gaussian `D_{Z, μ, σ}` with a leaf-dependent centre `μ` and standard
deviation `σ`. The one-dimensional draw is `SamplerZ` (Algorithm 15 of the
specification): it takes `r = μ − ⌊μ⌋`, draws a base sample and a sign, forms a
candidate `z`, and accepts it with a probability computed by `BerExp`. Its very
first line is
```
    s = ⌊μ⌋
```
and this `floor` is the mechanism the whole paper turns on: it is discontinuous,
so a perturbation of `μ` by 10⁻¹³ changes `s` by 1 whenever the two values
straddle an integer.

### 2.2 The specification's precision argument is distributional

The specification (§2.5.2) bounds the effect of finite precision on the
*distribution* the sampler outputs: with centres known to absolute error `δc`
and deviations to relative error `δσ`, the Rényi-divergence argument gives no
security loss when `δc + δσ ≤ 2⁻⁴⁶`, and the authors measure `≤ 2⁻⁴⁰` at 53
bits, calling residual leakage "a purely theoretic threat." This is correct and
is about the distribution. It says nothing about whether two implementations
that both stay within that precision produce the *same* sample on the *same*
randomness — and, as §6 shows, they do not always, because the quantity that
decides the sample passes through a `floor` the divergence argument never sees.

### 2.3 ePrint 2024/1709

Lin, Tibouchi, Yu and Zhang ("Do Not Disturb a Sleeping Falcon", EUROCRYPT
2025) identify precisely this sensitivity. Their Lemma 1 shows `SamplerZ` is
sensitive to a floating-point error on the centre only when the centre is
near-integer (and then a `floor` straddle flips the output by 1); their Lemma 2
shows it is essentially insensitive to an error on the standard deviation
(inconsistency probability `≤ 160ε`, so ~10⁻¹¹ at the ε ≈ 2⁻⁵² of real
discrepancies). Their Heuristic 1 locates the near-integer centres: at the first
two and last two of the `2n` calls, with probability `1/q` and
`1/‖(g,−f)‖²` respectively, and negligibly elsewhere. A divergence in the **last
two** calls differs in two coordinates over the trapdoor basis and yields full
key recovery (their §5); in the first two it yields only a short lattice vector.
They exhibit discrepancies from two sources: fused multiply-add, and the
reference's two signing modes `sign_dyn` / `sign_tree`, which order the bottom
of the tree traversal differently (their §6.1); their §7.2 gives a
countermeasure that aligns the two modes. Their attack matters for
**derandomized** FALCON (IBE key extraction, SNARK-friendly signatures,
aggregation), where the same input is signed twice.

The present paper's relationship to this one is stated in full in §7. In brief:
we reproduce their mechanism from an independent implementation, we add one
perturbation source they do not consider (which their countermeasure does not
fix), and we check their Heuristic 1 directly, which they do not.

---

## 3. Methodology

The measurement discipline here is not incidental; the project's log records six
occasions on which the measuring apparatus, not the implementation, produced the
wrong number, each caught by a control. We promote the discipline to a method.

- **Compare raw 64-bit patterns, not rounded values.** "Agrees to 10⁻¹²" is
  meaningless once `SamplerZ` rounds its input through `floor`: a difference far
  below 10⁻¹² can still cross an integer. All floating-point comparisons in
  §4–§6 are on `reinterpret(UInt64, ·)`.
- **Measure the effect of a change by reverting it.** The claim "respelling the
  operations is what reproduces C" was tested by reverting the respelling; it
  survived, so the claim was false (§4). The claim "the respellings do not reach
  the signature" was tested by measuring at a sample size that could see 10⁻⁵;
  it failed (§6).
- **Keep a positive and a negative control on the instrument, and confirm each
  arm actually perturbs what it claims to.** The divergence harness (§6) prints,
  before measuring, how many of the tree's 9728 doubles each arm moves; an arm
  that moves none is flagged as a no-op rather than reported as a null result.
- **Report C as a set of builds, not one number.** "Faster than C" with one
  `gcc -O2` build is not a comparison; the spread across compilers and
  optimization levels is itself data (§9).
- **Warm up a JIT before timing it, and do not disable the garbage collector.**
  The first call to any Julia function is compilation; timing it inflates the
  result by orders of magnitude. Allocation is a real cost of this
  implementation and is left in.

### 3.1 The instrument that made byte-exact comparison possible

The obstacle to reproducing the reference's signatures was not the arithmetic
but the randomness. The reference's public API derives the sampler's ChaCha20
stream internally from SHAKE, with no way to fix it from outside. We therefore
compiled the reference as a shared library with a shim
(`scripts/cref_shim.c`) that seeds the sampler's 56-byte ChaCha20 state
directly, and recorded, for eight keys, the reference's expanded key, its
sampled `s2`, and its compressed bytes (`scripts/gen_cref_sign_kat.jl`,
`test/vectors/cref_sign_kat.jl`). With the randomness pinned, "does our
implementation compute the same signature?" becomes a bit comparison that needs
no C compiler to re-check.

---

## 4. Byte-exact reproduction of the reference

Starting from an implementation written to the specification and the Python
reference — which produces valid FALCON signatures but different ones from the C
reference on the same input — we reconciled the two by a stage-by-stage bit
comparison (`scripts/cmp_cref_fp.jl`). The result splits into three groups:

- **Agree with no change:** `fft`, `ifft`, `split_fft`, `merge_fft`, `add_fft`,
  `sub_fft`, `mul_fft`, `adj_fft`, `mul·adj`, `mul·selfadj`, and the entire `B0`
  matrix of the expanded key — 5120 of 5120 doubles at `logn = 9`.
- **Agree once spelled the reference's way:** complex division, the `D11` entry
  of LDL\*, and the hand-unrolled bottom two levels of `ffSampling`. These are
  §5's first class.
- **A different storage convention:** the tree leaves hold the reciprocal
  width, and are multiplied by the `fpr_inv_sigma` table (§8).

With these reconciled, signing reproduces the reference **byte for byte**:
`test/test_falcon.jl` ("signing reproduces the C reference byte for byte")
replays the eight recorded vectors and asserts our `s2` equals C's `s2` and our
compressed bytes equal C's, on every vector, and that each is a signature that
verifies. Reproducing the *bytes* — not just the sampled integers — requires the
four `samplerz` calls of the hand-unrolled block to consume randomness in the
reference's order, so that order is verified against the real sampler, not
merely transcribed (`src/ffsampling.jl`, `_ffsampling_c4`).

**A caution that became a result.** It is tempting to conclude that respelling
those three operations is *why* signing reproduces C. It is not: with all three
reverted to the specification's spelling, signing still reproduces C on all
eight vectors, and over 480 signatures no coefficient differed. We reported that
as "the respellings do not reach the signature." That conclusion was wrong, and
§6 is the correction; the 480-signature sample was simply too small to see a
10⁻⁵ event.

---

## 5. Three classes of specification underdetermination

### 5.1 Class I — how an operation is spelled

Three computations are algebraically identical between the specification (with
the Python reference) and the C reference, and round differently. We show the
two implementation patterns for the first two as running code; the third is
described.

**(a) Complex division.** The specification gives no formula; a reader uses the
language's complex division. Julia's and Python's are both Smith's algorithm,
which scales by the larger of the real and imaginary parts to avoid overflow.
The C reference forms the reciprocal explicitly. Both patterns are in the
codebase, selected by a flag:

```julia
# Pattern S — the specification-reader's spelling (Julia's built-in /,
# i.e. Smith's algorithm).  Numerically the better choice.
ComplexF64[ f[i] / g[i] for i in eachindex(f) ]

# Pattern C — the reference's spelling (scripts/cref/fft.c, FPC_DIV):
# form the reciprocal, then one complex multiply.  This is the one the KATs
# encode.
@inline function _cdiv_cref(a::ComplexF64, b::ComplexF64)
    br = real(b); bi = imag(b)
    m  = 1.0 / (br * br + bi * bi)
    br *= m; bi *= -m
    ar = real(a); ai = imag(a)
    return ComplexF64(ar * br - ai * bi, ar * bi + ai * br)
end
```

These are equal in exact arithmetic and differ in the last bit; on the toy
inputs `a = 1+3i`, `b = 7+11i` they already differ, which the test at
`test/test_falcon.jl` asserts so that the downstream comparison is not vacuous.

**(b) The `D11` entry of LDL\*.** Algorithm 8 of the specification and the
Python reference compute `D11` with three multiplications; the C reference keeps
the quotient the other way round and folds, doing one. Equal because `G00` is
Hermitian and therefore real:

```julia
# Pattern S — Algorithm 8 / ffsampling.py.  Three multiplications.
L10 = div_fft(G[2, 1], G[1, 1])
D11 = sub_fft(G[2, 2], mul_fft(mul_fft(L10, adj_fft(L10)), G[1, 1]))

# Pattern C — Zf(poly_LDL_fft).  One multiplication; keeps mu = G01/G00.
G01 = adj_fft(G[2, 1])
mu  = div_fft(G01, G[1, 1])
L10 = adj_fft(mu)
D11 = sub_fft(G[2, 2], mul_fft(mu, G[2, 1]))
```

**(c) The bottom two levels of `ffSampling`.** The reference special-cases
`logn == 2`, inlining two `split`/`merge` pairs and folding their twiddle
factors into the constants `1/√2` and `1/√8`, one multiplication where the
generic recursion does two (`src/ffsampling.jl`, `_ffsampling_c4`; the generic
path is the `else` branch of `ffsampling_fft`).

All three are switched independently by
`with_spec_spelling(f; cdiv, ldl, ffsampling)`, defaulting to the reference's
spelling. This independence is the point: it lets §6 attribute a divergence to a
specific formula, which an earlier version of this work could not do because all
three sat behind a single flag.

### 5.2 Class II — a constant that cannot be derived

The reference multiplies every tree leaf by `fpr_inv_sigma[logn]`, a stored
`1/σ` per degree. This constant is not the correctly rounded reciprocal of the σ
that Table 3.3 of the specification prints (§8): +941 ulp at `logn = 9`, +18545
at `logn = 10`. An implementation that computes `1/σ` from the published σ
builds a different expanded key.

### 5.3 Class III — a smoothing parameter that is used but not printed

`σ` is Table 3.3's value, reproducible from equation (2.13) at
`ε ≤ 1/√(Q_s·λ)`. `σ_min` is printed as a value with no derivation; it is the
smoothing parameter of `Z` at `ε = 1/√(2^64·n³)`, a *different* ε that the text
never states, and assuming the two share an ε gives a value 10 % wrong (§8).

---

## 6. Which underdetermined choices reach the signature

We sign each of 100 keys many times with two spellings on identical key, message
and PRNG state, and count how often the resulting signatures differ
(`scripts/divergence_rate.jl`). Five arms:

| arm | perturbs | quantity | events | rate | Poisson 95 % |
|:----|:---------|:---------|-------:|-----:|:-------------|
| **A2** | complex division + `D11` (Class I a,b) | **centre** | **3 / 400000** | **7.5 × 10⁻⁶** | 1.6 × 10⁻⁶ … 2.2 × 10⁻⁵ |
| A2 · n = 1024 | complex division + `D11` | centre | **1 / 350000** | 2.9 × 10⁻⁶ | 7 × 10⁻⁸ … 1.6 × 10⁻⁵ |
| A1 | hand-unrolled bottom levels (Class I c) | centre | 10 / 525000 | 1.9 × 10⁻⁵ | 9.1 × 10⁻⁶ … 3.5 × 10⁻⁵ |
| A2a | complex division alone | centre | 0 / 100000 | 0 | … 3.7 × 10⁻⁵ |
| A2b | `D11` alone | centre | 0 / 100000 | 0 | … 3.7 × 10⁻⁵ |
| B | `fpr_inv_sigma` vs `1/σ` (Class II) | **width** | 0 / 100000 | 0 | … 3.7 × 10⁻⁵ |

The rates are accumulated from independent foreground chunks
(`scripts/divergence_accum.jl`); background jobs do not survive this
environment's between-turn suspension (docs/debug_log.md #057), and summing
independent Poisson chunks is exact.

Before measuring, the harness reports how many of the tree's 9728 doubles each
arm moves — A2 moves 6800, A2a 6345, A2b 5400, B 452 — so none of the zero
results is a silent no-op. When a divergence occurs, ≈ 470 of the 512
coefficients differ, because the byte stream desynchronizes.

**A1 is a positive control**: it is exactly the `sign_dyn` / `sign_tree`
difference of ePrint 2024/1709 §6.1, and its rate (1.9 × 10⁻⁵ over 525000
signatures) sits on their Table 2 (≈ 3 × 10⁻⁵). **B is a negative control**, and
its zero is *expected from theory*: it perturbs the width, which Lemma 2 of that
paper bounds at ~10⁻¹¹ per signature — well below our resolution, so the zero
carries no information and we do not present it as a finding.

**A2 and A1 are not distinguishable in rate.** A2 at 7.5 × 10⁻⁶ and A1 at
1.9 × 10⁻⁵ differ by a point-estimate factor of 2.5, but a two-sample Poisson
test (conditional binomial on the 13 pooled events) gives p ≈ 0.23: not
significant. We therefore do **not** claim A2 is rarer than A1; both are of
order 10⁻⁵, which is consistent with their sharing the mechanism. An earlier
version of this work quoted A1 at 5 × 10⁻⁵ from a single 5-in-100000 point —
a high fluctuation; the 525000-signature estimate is 1.9 × 10⁻⁵.

**The finding replicates at n = 1024.** One divergence in 350000 signatures,
and reproducing it (`scripts/divergence_accum.jl` chunk 22, key 34, sig 454)
puts it at call 2047 of 2048 — the last two calls again — with `floor` straddling
427/426 and the leaf width bit-identical between runs. The mechanism is
degree-independent, as expected; the rate is comparable because
`1/‖(g,−f)‖²` is comparable at the two degrees (both ≈ 6 × 10⁻⁵).

**A2 is the arm that is not in that paper**, and it is not zero. All three of
its divergences carry the mechanism exactly
(`scripts/first_divergence.jl`):

```
key 70 sig  534   call 1023 of 1024   mu -267.00000000000006  vs -267
key 65 sig 1239   call 1023 of 1024   mu  337                 vs  336.99999999999994
key 68 sig 1885   call 1024 of 1024   mu  285.99999999999989  vs  286
```

Three of three at call 2n−1 or 2n; three of three straddling an integer with
`floor` on opposite sides; three of three with the leaf standard deviation
bit-identical between the two runs — so the divergence is carried entirely by
the centre. This is Lemma 1, not Lemma 2, at the position where §5 of that paper
recovers the whole private key from a single pair.

### 6.1 The divergence recovers the private key

"At the key-recovery position" is a claim about consequences, and we discharge
it rather than cite it. Taking the A2 divergent pair `(s, s')` at key 70,
signature 534 — two signatures on the same syndrome, differing only in the last
two sampler outputs — we run ePrint 2024/1709 §5.1's recovery
(`scripts/key_recovery.jl`). The signature difference is `Δs0 = Δz0·g`,
`Δs1 = −Δz0·f` with `Δz0 = a + b·x^{n/2}`, and since
`(a + b x^{n/2})^{-1} = (a − b x^{n/2})/(a² + b²)` the recovery is a search over
`(a, b) ∈ {−19,…,19}²` — no ring inversion. From `(s, s')` and the public key
`h` **alone** (not the secret, not `z`), the search recovers a short `(f, g)`
with `max|g| = 14`, `max|f| = 13` — matching the true key's coefficient sizes —
that reproduces `h` mod `q`. It equals the stored secret up to the negation
symmetry `(f,g) ↦ (−f,−g)`, which is the same signing key. **This is full key
recovery from a single A2 discrepant pair.** Complex division and `D11` — which
ePrint 2024/1709 does not consider — therefore suffice not merely to change the
signature but to expose the key, in the derandomized settings where the same
syndrome is signed twice.

### 6.2 The near-integer centres, measured directly

ePrint 2024/1709 measures the *consequence* of Heuristic 1 (its Tables 2, 4) and
notes (Remark 1) that the heuristic cannot be made a theorem. We measure the
*cause* — how often each call's centre is within 10⁻⁹ of an integer — in an
implementation written from the specification, not derived from theirs
(`scripts/heuristic1_check.jl`):

| position | measured | Heuristic 1 |
|:---------|---------:|------------:|
| calls 1, 2 (vary message) | 17 / 200000 = 8.5 × 10⁻⁵ | 8.14 × 10⁻⁵ |
| calls 2n−1, 2n (vary PRNG) | 14 / 200000 = 7.0 × 10⁻⁵ | 6.08 × 10⁻⁵ |
| calls 3 … 2n−2 | 0 / 102000000 | ≈ 0 |

Both ends agree with the heuristic inside Poisson error; the interior is empty
across 1.02 × 10⁸ draws, as the `q²` denominators require. The two ends must be
measured with *different* experiments — varying the message for the first two,
the PRNG for the last two — because the centres of calls 1 and 2 do not depend
on the signature randomness at all (one distinct value over 8 signatures, where
call 3 takes seven). Overlooking that determinism briefly read as a
contradiction with the heuristic, and was a sampling error in our experiment,
not a fact about FALCON.

A separate confirmation falls out: all 100 keys have **even** `‖(g,−f)‖²`,
reproducing the idiosyncrasy that §7.1 of that paper attributes to the C
reference's key generation (and which blocks its own countermeasure), from an
implementation that followed the C key generation without knowing it mattered.

### 6.3 Why the last two calls dominate — a mechanism for an open question of 2024/1709

Section 6.1 of ePrint 2024/1709 reports, and states it does "not fully
understand," that conditional on an integer centre the **last** two calls
diverge more readily than the first two, specifically for the `sign_dyn` /
`sign_tree` (our A1) difference and not for FMA. We give the mechanism, by
measuring the perturbation `|Δμ| = |μ_C − μ_spec|` the respelling puts on the
centre at **every** position, over 10⁴ signatures per position
(`scripts/position_profile.jl`).

For A1, the perturbation is not uniform across the traversal — it accumulates,
with a sharp step at the midpoint:

```
positions   0.. 511  (z1, first descent) : mean |Δμ|  2.7e-15 → 1.0e-14
positions 512..1023  (z0, second descent): mean |Δμ|  7.4e-14 → 8.7e-14
first two  (0,1):   mean |Δμ| = 2.7e-15
last two (1022,1023): mean |Δμ| = 8.6e-14      ratio ≈ 32×
```

The step is structural. `ffSampling` samples `z1` first, directly from the
split target; then forms `t0 ← t0 + (t1 − z1)·l10` and samples `z0`. Every
centre in `z0` (the second half, including the last two calls) therefore passes
through the top-level `l10` correction, which carries the bottom-level rounding
difference; the centres in `z1` (the first half, including the first two calls)
do not. A centre diverges only when `|Δμ|` exceeds its distance to the nearest
integer, so a 32× larger `|Δμ|` yields a ~32× higher conditional divergence
probability at the last two calls — exactly the asymmetry the paper observed.
It is A1-specific because FMA perturbs every floating-point operation
throughout, giving a flat `|Δμ|` profile with no such step, consistent with the
paper's report that FMA shows no asymmetry.

The same instrument sharpens the A2 result. For A2 (the tree-only respelling),
`|Δμ|` at the first two calls is **exactly zero** over all 10⁴ draws:

```
A2   first two: mean |Δμ| = 0.000 (exactly)   last two: 1.1e-14
```

The first leaf's centre is a pure split of the target `t` (built from `B0`, not
the tree), reached before any `l10` correction applies, so a respelling that
changes only the tree cannot move it. A2 divergences are therefore confined to
the last two calls **structurally**, not merely probabilistically — which is
why all four A2 divergences we observed (three at n = 512, one at n = 1024) are
at position 2n−2 or 2n−1 and none at 0 or 1. The sensitive set is not a fixed
property of FALCON; it depends on where the perturbation source enters the
computation.


A follow-up experiment sharpens this into a quantitative rule
(`scripts/first_two_probe.jl`). Because the first-two centre is `n/q`, an
attacker who chooses messages can *manufacture* an exact integer centre there,
one message in about `q/2`. We did: 10 integer centres in 48000 (key, message)
pairs, matching Heuristic 1's `2/q`. **None of them diverged.** Printing both
runs' centres shows why: the computed centre is never exactly the integer, it
sits 1–3 ulp away, and — crucially — *both spellings sit on the same side of it*,
because their difference is smaller than that shared offset:

```
call 1   C 22.000000000000007   spec 22.000000000000007   |Δμ| = 0
call 2   C 42.000000000000021   spec 42.000000000000014   |Δμ| = 7.1e-15
call 2   C -47.000000000000021  spec -47.000000000000014  |Δμ| = 7.1e-15
```

So an integer centre is necessary but not sufficient: the straddle also needs
`|Δμ|` to exceed the spread of the shared rounding offset (~2 ulp ≈ 2.8×10⁻¹⁴ at
these magnitudes). That gives

    P(straddle | integer centre) ≈ |Δμ| / spread

which is ≈ 0.1 at the first two calls (`|Δμ| ≈ 2.7×10⁻¹⁵`) and saturates at the
last two (`|Δμ| ≈ 8.6×10⁻¹⁴`). This is how §6.3's 32× perturbation ratio becomes
a divergence-rate ratio, and it is consistent both with our 0-of-10 and with
2024/1709's Table 2, where about 30 % of the dyn/tree discrepancies fall at the
first two calls. It also means a chosen-message adversary gains nothing at the
first two positions — the only exploitable end is the last two.

### 6.4 The first two centres are a public linear form in the secret

One structural fact underlies both Heuristic 1's first-two case and a different
attack shape. In `sample_preimage` the second target component is
`t1 = (−ĉ·b)/q` with `b = B0[1,2] = −f`, so as a ring element `t1 = (c·f)/q`
exactly, where `c` is the *public* hashed message. The first descent samples
`z1` from `t1`, so the first two sampler centres are coefficients of `(c·f)/q`.
Checked against exact `Int128` negacyclic arithmetic (`scripts/linear_form.jl`),
they are, at fixed indices and to floating-point error alone:

```
call 1 : mu = 41.9012124664   (c*f)[256]/q = 41.9012124664   |diff| = 0
call 2 : mu =  8.6159166734   (c*f)[512]/q =  8.6159166734   |diff| = 2.5e-14
```

Two consequences. First, the *event* is exactly characterised at these two positions: "integer
centre" is precisely `(c·f)_{n/2−1} ≡ 0 (mod q)`, with no heuristic involved,
because the denominator is exactly `q` and the numerator is an explicit integer.
Turning that into Heuristic 1's *probability* `1/q` needs, in addition, that
`(c·f) mod q` be uniform over the random oracle's output and that `f` be
invertible mod `q`; we do not prove either, so we claim the characterisation of
the event, not a theorem replacing the heuristic.

Second, that condition is **linear in the secret, with public coefficients**.
Each *detected* integer-centre event yields one linear equation on `f mod q`:
writing the coefficient of `x^k` in `c·f` over `Z[x]/(x^n+1)` as `⟨c, x^k f*⟩`
with `f*` the adjoint (`f*(x) = f(x^{-1})`), the condition at call 1 is
`⟨c, x^{n/2−1} f*⟩ ≡ 0 (mod q)`. `n` independent events determine `f mod q`, and
since `‖f‖ ≪ q` that determines `f`. This is a different shape from §5.1, which reads the key out of
the difference *vector* of two signatures — here the leak is the *event*, and no
difference vector is needed. It applies at the first two calls, the positions
2024/1709 dismisses as yielding only a short lattice vector.

We are explicit about what is missing: an oracle for the event. The respelling
perturbation does not provide one — §6.3 shows `|Δμ|` at the first two calls
sits below the shared rounding offset, and `scripts/first_two_probe.jl`
manufactured 10 integer centres there with zero straddles. A larger perturbation
source (FMA), a conformance-mismatch report, or a side channel would. We state
the shape and its requirement; we do not claim a working attack on plain Falcon.

### 6.5 Two directions that did not pan out (recorded honestly)

Two hypotheses we tested and rejected, since the boundary they probe is part of
the result. **A weak-key class by `‖(g,−f)‖²`**: the last-two rate is
`1/‖(g,−f)‖²`, but over 3000 keys `‖(g,−f)‖²` is tightly bounded (15078 to
16822, a 1.09× spread), because key generation rejects anything above `1.17²q`.
The rate is essentially key-independent; there is no weak-key tail on this
invariant. **A fifth sensitive position in the interior, for random messages**:
over ~10⁷ interior draws the near-integer (< 10⁻⁶) count is 21, exactly the ~20
expected from uniform fractional parts alone, i.e. no integer centres occur in
the interior — confirming 2024/1709's restriction to the first and last two
positions holds for random messages. (Whether a chosen-message adversary can
reach the interior positions, whose denominators are within double precision up
to the first and last *six* calls, is left open.)

### 6.6 The leak is structural, not generic rounding

Injecting a uniform perturbation `ε` at every sampler centre and measuring the
divergence rate (`scripts/precision_law.jl`) confirms Lemma 1's linear law
empirically — `rate ≈ 2n·ε`, saturating to 1 for `ε ≳ 1/2n`:

| ε          | 2.4e-7 | 9.5e-7 | 3.8e-6 | 1.5e-5 | 6.1e-5 | 9.8e-4 |
|:-----------|-------:|-------:|-------:|-------:|-------:|-------:|
| rate       | 4.2e-4 | 3.1e-3 | 1.4e-2 | 5.4e-2 | 1.6e-1 | 9.2e-1 |

The instructive comparison is with the real respellings. A uniform `ε = 10⁻¹⁴`
would give `rate ≈ 2n·ε ≈ 10⁻¹¹`; the actual respellings, whose perturbation is
the same size (`|Δμ| ≈ 10⁻¹⁴`, §6.3), produce `≈ 10⁻⁵` — six orders of
magnitude more. The difference is not the size of the perturbation but where it
lands: a uniform injection mostly hits generic centres (which straddle with
probability `ε`), whereas the respelling's perturbation reaches the same
structured **exact-integer** centres that make the mechanism work. The leak is
driven by that structure, not by generic floating-point noise — which is why
improving precision does not remove it (an exact-integer centre is sensitive to
an arbitrarily small perturbation) and only the structural fix of §7.1 does.

---

## 7. Relationship to ePrint 2024/1709

We state this in full because the honest positioning is most of the paper's
contribution.

| | ePrint 2024/1709 | this work |
|:--|:--|:--|
| perturbation source | FMA; and the reference's two signing modes (`sign_dyn`/`sign_tree`) | the **specification vs the reference** |
| mechanism | Lemma 1 (centre) / Lemma 2 (width) | **the same — reproduced independently** |
| hand-unrolled difference (Class I c) | identified §6.1, countermeasure §7.2 | reproduced as positive control (A1) |
| complex division, `D11` (Class I a,b) | **not present** | **7.5 × 10⁻⁶ at the key-recovery position; not fixed by §7.2** |
| Heuristic 1 | consequence measured | **cause measured directly (§6.2)** |
| last-two-vs-first-two asymmetry | observed, "not fully understood" (§6.1) | **mechanism given: depth-accumulated `l10` perturbation, 32× (§6.3); not calibrated (§12.3)** |
| key recovery | described (§5) | **run on our own A2 pair (§6.1)**, and on a pair produced with part 1 of the §7.1 countermeasure deployed (§7.1 iv, arm A1) |
| even `‖(g,−f)‖²` (§7.1) | noted as a C idiosyncrasy | reproduced independently |

The one genuinely new empirical fact is the A2 arm. Its significance is not that
it is a new mechanism — it is Lemma 1 — but that its *source* is different.
ePrint 2024/1709's two sources both live inside the reference implementation:
FMA is a build flag, and `sign_dyn`/`sign_tree` are two entry points of the same
code, which their §7.2 countermeasure reconciles with each other. Complex
division and `D11` differ between the **specification** and the reference, so
they are present in *both* signing modes and survive that countermeasure. For a
standard that intends bit-exact KATs, this is the difference between "the two
implementations should be made to agree" and "the text must say which of two
equal formulas to use."

We were, for three revisions, wrong about this in the direction of claiming too
much: an earlier framing ("specification underdetermination alone suffices — no
different build needed") was a rediscovery of §6.1, because a single flag we had
not read carefully was moving only Class I(c). Reading the paper, and then
separating the three spellings, is what produced the correct and narrower claim.

---

### 7.1 The proposed countermeasure, implemented and evaluated

Section 7.1 of ePrint 2024/1709 proposes the structural fix: replace `SamplerZ`
by `NewSamplerZ` (Algorithm 4), which splits the centre with **round** instead
of **floor**, moving the instability from the integers to the half-integers; and
restrict `‖(g,−f)‖²` to be **odd** in key generation, which makes half-integer
centres impossible. The paper notes, in one sentence, that the C reference only
generates keys with `‖(g,−f)‖²` even, calling this "an idiosyncrasy of the C
implementation itself, that is easily fixable." We implemented the countermeasure
and evaluated it (`scripts/countermeasure_eval.jl`,
`scripts/denominator_check.jl`). Three findings.

**(i) The countermeasure is distributionally sound.** We implemented Algorithm 4
faithfully — `r ← c − round(c)`, `y ← (2b−1)y₊` (not `b + (2b−1)y₊`),
`x ← (y−r)²/2σ² − (y₊²−y₊)/2σ_max²` — building `NewBaseSampler`'s reverse CDT at
256 bits from the weights `w(0) = ½`, `w(i) = exp(−(i²−i)/2σ_max²)`. A χ² test
over 200000 samples at each of five `(μ,σ)` pairs gives `χ²/df` of 0.72, 0.96,
0.56, 0.62 and 1.02 — the rounding sampler produces the correct discrete
Gaussian. This is an independent confirmation of that paper's Algorithm 4.

**(ii) The rational structure of the centres, verified directly.** Theorem 1 and
Heuristic 1 assert the exact centre is `n_k/g_k` with `g₀ = q` at the first two
calls and `g_{n−1} = ‖(g,−f)‖²` at the last two, but check this only indirectly.
Measuring `|g·c − round(g·c)|` over 100000 samples at each end gives a maximum
residual of 1.4×10⁻⁹ and 2.4×10⁻⁸ against values of magnitude ~5×10⁶ — a relative
2×10⁻¹⁵, i.e. exactly integral up to floating-point error. The structure holds.

**(iii) Deployed on the reference key generator, part 1 immunises only the
harmless positions.** This is the consequence the paper does not draw. A centre
`n/g` is an integer iff `g | n`, which is possible for any `g`; but it is a
half-integer iff `2n = g(2j+1)`, which is possible **only when `g` is even** (for
odd `g`, `g | 2n` and `gcd(g,2)=1` force `g | n`, an integer). Now:

- at the **first two** calls `g₀ = q = 12289`, which is **odd**, so half-integer
  centres cannot occur at all — rounding immunises these positions completely;
- at the **last two** calls `g_{n−1} = ‖(g,−f)‖²`, which the reference key
  generator always makes **even**, so half-integer centres occur at the same
  density `1/g` that integer centres did — rounding immunises them not at all.

Measured, over 100000 samples at each end:

| position | `g` | integer centres | half-integer centres |
|:---------|:----|----------------:|---------------------:|
| calls 1, 2 | `q = 12289` (odd) | 5 | **0** |
| calls 2n−1, 2n | `‖(g,−f)‖²` (even) | 7 | **11** |

(Predicted densities `1/q = 8.1×10⁻⁵` and `1/t = 6.1×10⁻⁵` give expectations
8.1, 6.1 and 6.1; the observed 11 is high at `p ≈ 0.04` but of the predicted
order.) And the first two calls are precisely the *harmless* ones — §5 of that
paper obtains from them only "a short lattice vector… not expected to be short
enough to make key recovery feasible" — while the last two are the full
key-recovery positions. So under `floor`, 7 of 12 exposures sit at the dangerous
end; under `round` without the key-generation change, **all 11 do**. The
countermeasure does not reduce key-recovery exposure; it concentrates what
remains onto the positions that matter.

**(iv) With part 1 alone deployed, key recovery still succeeds.** We ran the
whole attack against it (`scripts/countermeasure_break.jl`): sign with
`NewSamplerZ` on reference keys, look for a half-integer straddle at the last
two calls, and when one occurs run §5.1's recovery. Over 160000 signatures, 29
centres came within 10⁻¹² of a half-integer, **3 produced divergent signature
pairs, and all 3 yielded the private key** — each verified to be exactly the
stored key up to the NTRU lattice's symmetries (`x²⁵⁶(f,g)`, `x⁰(f,g)`,
`−x²⁵⁶(f,g)`).

**Disclosure, because it materially qualifies this result:** the perturbation
source used here is arm **A1** (the hand-unrolled bottom levels), which is the
`sign_dyn`/`sign_tree` difference that the *other* countermeasure of that paper,
§7.2, removes. So the configuration demonstrated is "deploy half of §7.1, do not
deploy §7.2, and attack with the source §7.2 kills" — not a configuration an
implementer following that paper would be in. Our own arm A2 survives §7.2, but
has never been observed under `NewSamplerZ`: its expected count in 160000
signatures is ≈1.2, so demonstrating it would take of order 10⁶ signatures,
which we have not run. What the experiment does establish is the mechanism —
that rounding leaves the last-two exposure intact — not a break of the
countermeasure as a whole.

```
DIVERGENCE key 20 sig  562 : 909 of 1024 coefficients differ
   *** KEY RECOVERED *** (a,b)=(-1,0)  = x^256 (f,g)
DIVERGENCE key 41 sig 1346 : 900 of 1024 coefficients differ
   *** KEY RECOVERED *** (a,b)=(-1,0)  = x^0 (f,g)
DIVERGENCE key 82 sig 1054 : 936 of 1024 coefficients differ
   *** KEY RECOVERED *** (a,b)=(-2,0)  = -x^256 (f,g)
```

The end-to-end rate is unchanged: 3 in 160000 (1.9×10⁻⁵) with the countermeasure
deployed, against the `floor` baseline's 10 in 525000 (1.9×10⁻⁵).

**(v) Part 2 is deployable, so the point is that both parts must be mandated
together.** An earlier version of this work claimed part 2 was unavailable,
on the strength of observing 0 odd-norm keys in 3200. That was an argument from
absence, and it was wrong. The reason the reference never produces one is exact:
`gen_poly_cdt` forces the coefficient sum of **both** `f` and `g` odd (so that
`Res(f, x^n+1)` is odd and the binary GCD at the bottom of the Pornin–Prest
descent does not fail on a factor of 2), whence
`‖(g,−f)‖² = Σg_i² + Σf_i² ≡ g(1) + f(1) ≡ 0 (mod 2)` in one line. But
2024/1709 says only *one* of the two needs to be odd, and that is correct: with
`f` odd and `g` even (`scripts/odd_norm_keygen.jl`) the solver succeeds **40 of
40** and every resulting key has an odd norm, against 0 of 40 in the control.
Their "easily fixable" stands.

The normative consequence is therefore not that the countermeasure fails, but
that **its two parts are not independently useful**. The sampler change is the
obvious one and the key-generation parity change is easy to overlook; an
implementer who makes only the first gets, by the parity argument above, no
protection at exactly the positions that leak the key. A standard adopting
Algorithm 4 must mandate the key-generation change in the same clause and say
explicitly that neither part alone suffices.

With `t = ‖(g,−f)‖²` odd,
`m₂ = t² − 2u²` and `m₃ = t³ − 2t(u²+v²+w²) + 2u(v−w)²` are odd as well
(odd − even = odd), so all six in-precision denominators are odd and half-integer
centres cannot occur anywhere; that argument is sound. The gap is one of
deployment, not of mathematics: part 2 is unavailable on the reference key
generator (0 odd keys in 3200 generated), so what an implementer can actually
deploy today is part 1 alone, and part 1 alone leaves the key-recovery exposure
untouched. Any FIPS 206 text that adopts the rounding sampler **must** adopt the
key-generation parity change in the same breath, and must say so explicitly.

### 7.2 A countermeasure that needs no key-generation change

The break in §7.1 is a deployment failure, not a mathematical one, and it
suggests its own repair. At the sensitive positions the exact centre is a
rational `n/g` whose denominator the **signer already knows**: `g = q` at the
first two calls, `g = ‖(g,−f)‖²` at the last two. The floating-point centre
`μ̂` approximates `n/g` to about 10⁻¹³, so `g·μ̂` sits within `g·10⁻¹³ ≈ 1.6×10⁻⁹`
of the integer `n`. Therefore

    n = round(g · μ̂)        recovers the numerator exactly,
    s = fld(n, g)           is then exact integer arithmetic,

and the split of the centre stops depending on floating point at exactly the
places where that dependence is dangerous. The margin is not marginal: §7.1(ii)
measured `max |g·μ̂ − round(g·μ̂)|` at 1.4×10⁻⁹ (first two calls) and 2.4×10⁻⁸
(last two) over 100000 samples each, against the 0.5 that would be needed to
break the rounding — seven orders of magnitude of slack.

Replaying every divergent pair this project has found, with the exact split off
and then on (`scripts/snapping.jl`):

```
A2 key 70 sig  534   off: 463 of 512 differ  DIVERGE   on: 0  AGREE
A2 key 65 sig 1239   off: 498 of 512 differ  DIVERGE   on: 0  AGREE
A2 key 68 sig 1885   off: 457 of 512 differ  DIVERGE   on: 0  AGREE
A1 key  4 sig  105   off: 455 of 512 differ  DIVERGE   on: 0  AGREE
```

And it changes nothing else: over 1200 ordinary signatures with the same
spelling and the same PRNG state, snapping altered **0** of them. That is the
point — it does not move the split, it only makes the split exact. The
fractional part `r = μ − s` is still computed in floating point; what Lemma 1
establishes is that the sampler is sensitive *only* at the discontinuity, so a
perturbation of `r` that does not cross it changes the execution with
probability of order that perturbation. (An earlier draft cited Lemma 2 here;
Lemma 2 is about the standard deviation, and the relevant statement is the
continuous half of Lemma 1.)

Compared with §7.1 this is a single change rather than two that must be adopted
together — which, given that §7.1's two parts are not independently useful, is
the practical difference:

| | 2024/1709 §7.1 | rational snapping |
|:--|:--|:--|
| sampler change | replace `SamplerZ` with Algorithm 4 | 3 lines at 4 of 1024 calls |
| key-generation change | **required** (odd `‖(g,−f)‖²`); deployable, but a second change an implementer must not omit | none |
| output distribution | new base sampler, new proof obligation | unchanged (0 of 1200 signatures moved) |
| effect on the deployable path | key recovery still succeeds (§7.1 iv) | all known divergences removed |

Four honest caveats, all of which a standards body would have to see closed
before adopting this. First, the 10⁻⁹ margin is **measured, not proved**; a
rigorous forward-error bound on the centre computation is needed before a
standard could rely on it. Second, snapping is only available where `g` is both
known and below the floating-point precision — the first six and last six calls;
in the interior `g ≳ q³` exceeds `2⁵³`, but there an integer centre occurs with
probability below 10⁻¹⁶ and §6.2 measured zero in 1.02×10⁸ draws, so those
positions do not need protecting. Third, the "0 of 1200 signatures moved" check
has **almost no statistical power**: with near-integer densities of ≈8×10⁻⁵ and
≈6×10⁻⁵ at four snapped positions, the expected number of changes in 1200
signatures is ≈0.3, so the observation confirms only that snapping is a no-op
away from a boundary, which holds by construction. A distributional test with
the power of the χ² we ran for Algorithm 4 has not been done. Fourth, our
implementation (`fld(round(Int128, big(g)*mu), big(g))`) allocates and is not
constant-time, and snapping with `floor` admits `r = μ − s` outside `[0,1)` when
`μ̂` falls just below the true integer centre; a deployable version needs a
constant-time, allocation-free formulation and a statement of what `BerExp` does
with a negative `r`. We present snapping as a promising direction, not a
finished countermeasure.

### 7.3 The floating-point environment: a determinism failure the specification never mentions

Every perturbation studied so far rewrites *source*. One does not. The IEEE-754
rounding direction lives in the x87 control word and in MXCSR; it is per-thread
process state; and **nothing in the FALCON specification, the C reference, or
this implementation ever sets it.** Any library in the address space may change
it, and some do.

With the key and the PRNG state held identical and only the rounding direction
changed for the signing computation (`scripts/rounding_mode.jl`):

| mode | signatures differing |
|:--|--:|
| FE_TONEAREST (control) | **0 of 1800** |
| FE_UPWARD | **1800 of 1800** |
| FE_DOWNWARD | **1800 of 1800** |
| FE_TOWARDZERO | **1800 of 1800** |

The rate is 1. Against the 1.9×10⁻⁵ of the strongest source-level difference we
measured, this is five orders of magnitude larger, and it needs no second
implementation at all — one implementation suffices, run twice in the same
binary. The consequence for FIPS 206 is direct: a signature's bytes are not a
function of (key, message, randomness); they are a function of (key, message,
randomness, **rounding mode**). A bit-exact known-answer-test requirement is
unsatisfiable unless the standard also fixes the floating-point environment,
and no draft text does. (ePrint 2024/1709 §6 mentions the "weak determinism" of
floating point as "a first way, which we do not explore further"; this
quantifies it, and it is by far the largest effect in the class.)

**It does not, however, make the attack easier — and that bounds the whole
attack family.** A whole-signature flip desynchronises the sampler, so the
difference is unstructured, which §5 of that paper already identifies as
useless for key recovery. Confining the flip to a brief window at the tail of
the traversal, so that only the last two centres are perturbed
(`scripts/rounding_attack.jl`), gives 0 recoveries in 48000 signature pairs —
consistent with the ~10⁻⁵ rate we measure for source-level perturbations, not
better. The reason is structural: key recovery needs the *exact* centre at one
of the last two calls to be an integer, and that is a property of the key and
the message, not of the perturbation. Its probability is `2/‖(g,−f)‖² ≈
1.2×10⁻⁴`, and no perturbation, however large, exceeds it. **Every attack in
this family is bounded above by that ceiling**, so an adversary needs of order
10⁴ signature pairs on the same syndrome whatever tool they bring — a bound
worth stating in a risk assessment, and one that also explains why increasing
the perturbation strength (§6.6) buys so little.

---

## 8. The constants

### 8.1 `fpr_inv_sigma` is not a reciprocal of the published σ

Computed at 400 bits and compared on raw bit patterns
(`scripts/inv_sigma_audit.jl`):

| derivation of `1/σ` | logn = 9 | logn = 10 |
|:--|--:|--:|
| reciprocal of Table 3.3's printed σ | +941 ulp | +18545 ulp |
| reciprocal of the Python reference's Float64 σ | +1 ulp | +2 ulp |
| Float64 `1.0 / σ` | +1 ulp | +2 ulp |
| reciprocal of (2.13) at 400 bits | +1 ulp | +1 ulp |
| (2.13) in Float64, then `1.0 /` | +2 ulp | +1 ulp |

The table is above every candidate, and its 35-digit decimal literal is not the
35-digit expansion of `1/σ` (it agrees to ~16 figures then diverges), so it
appears to be some Float64 printed at length. An implementation that recomputes
the reciprocal differs from the reference in 452 of 512 tree leaves — but, as
arm B shows and Lemma 2 explains, not in the signature. This constant is a KAT
issue for intermediate values and the expanded-key format, not for signatures.
(A related documentation slip: `inner.h:749` calls `fpr_sigma_min[]` the array
`1/σ_min`; it holds `σ_min`.)

### 8.2 `σ_min`'s ε is not printed

`η_ε(Z) = (1/π)·√((1/2)·ln(2(1 + 1/ε)))` at `ε = 1/√(2^64·n³)` reproduces both
listed `σ_min` to 3 × 10⁻¹³. That ε is not the ε of (2.13); conflating them (the
`4n(1+1/ε)` of (2.13) is the smoothing of `Z^{2n}`, the `2(1+1/ε)` here is the
smoothing of `Z`) gives a value 10–11 % off. `test/test_spec.jl` verifies both
numerically.

### 8.3 The SamplerZ test-vector byte convention

Table 3.2's sixteen `SamplerZ` vectors assemble nine random bytes into a 72-bit
word in an order the text does not give; read left-to-right four of sixteen
reproduce, with each chunk reversed all sixteen do. The convention lives only in
the reference test harness (`bytes.fromhex(oc)[::-1]`).

---

## 9. Performance

One machine (4 cores, x86-64 with FMA; the single-machine limitation is real and
we do not hide it). We vary what we *can*: two compilers, three optimization
levels, emulated vs native floating point for C; three optimization levels for
Julia; and we report the run-to-run distribution. Full data:
`docs/benchmarks.md`; harness: `scripts/bench_all.sh`, `scripts/cref_bench.c`,
`scripts/bench.jl`, `scripts/pyref_bench.py`.

### 9.1 FALCON-512, median milliseconds per operation

| build | keygen | sign | verify |
|:--|--:|--:|--:|
| **this work, Julia -O2** | **9.78** | **0.89** | **0.015** |
| this work, Julia -O3 | 11.50 | 1.03 | 0.016 |
| this work, Julia -O1 | 15.49 | 1.44 | 0.047 |
| C native, clang -O3 -march=native (fastest) | 5.02 | 0.135 | 0.017 |
| C native, gcc -O2 | 5.91 | 0.171 | 0.026 |
| C emulated, clang -O3 (Algorand-style deploy) | 10.71 | 1.58 | 0.020 |
| C emulated, gcc -O2 | 12.33 | 1.79 | 0.024 |
| Python reference | 2657 | 18.39 | 3.52 |

Reading the ratios from the Julia -O2 row:

- **Verification: this implementation is the fastest measured, C included** —
  0.015 ms against the best C build's 0.0165 ms (≈ 1.1×) and gcc -O2's 0.026 ms
  (≈ 1.6×), and **235× faster than the Python reference**. Verification is
  NTT-heavy integer work with no sampler and no big integers, which is where
  Julia's arrays and LLVM are at their best.
- **Signing: ≈ 1.8× faster than the emulated C build** (the one most deployed
  FALCON actually runs, since native FP is discouraged for determinism), **20×
  faster than Python**, and **≈ 6× slower than native-FP C**. The gap to native
  C is the floating-point sampler, where C's SSE2 doubles beat going through
  Julia's dispatch and the `Ref`-guarded spelling flags.
- **Key generation: ≈ 9.78 ms, faster than every emulated C build** (10.7–12.3
  ms) and **272× faster than Python**, but **≈ 1.9× slower than native C** (5.0
  ms). Key generation carries its big integers in `BigInt`, against the
  reference's hand-rolled 31-bit-limb RNS; that is the residual gap.

Note that the best Julia optimization level is not fixed: at n = 512, `-O2`
beats `-O3` on all three operations, while at n = 1024 (§9.2) `-O3` wins on
keygen and sign. `-O3` enables LLVM passes (aggressive vectorization, loop
transforms) whose payoff depends on the working-set size, so it helps the larger
degree and hurts the smaller. Reporting a single optimization level would have
misrepresented the implementation in one direction or the other; we report all
three and take the best per operation.

The key-generation figure is also the endpoint of a long optimization: the first
working version took 9610 ms, and the path to 9.78 ms was almost entirely the
removal of allocation (a callable-struct buffer to kill a closure box, fused
GMP multiply-accumulate with no temporaries, schoolbook multiplication which
beats Karatsuba at every degree from 2 to 2048 here), not high-level array
rewriting (`docs/debug_log.md` #032–#044).

### 9.2 FALCON-1024, median milliseconds per operation

| build | keygen | sign | verify |
|:--|--:|--:|--:|
| **this work, Julia -O3** | **42.04** | **2.05** | **0.031** |
| this work, Julia -O2 | 46.85 | 2.12 | 0.031 |
| this work, Julia -O1 | 65.58 | 3.23 | 0.094 |
| C native, clang -O3 -march=native (fastest) | 14.86 | 0.276 | 0.035 |
| C native, gcc -O2 | 19.74 | 0.339 | 0.055 |
| C emulated, clang -O2 (deploy) | 31.45 | 3.43 | 0.042 |
| C emulated, gcc -O2 | 35.98 | 3.85 | 0.056 |
| Python reference | 12457 | 38.52 | 7.48 |

The pattern holds and sharpens with degree. **Verification is again the fastest
measured** — 0.031 ms against the best C build's 0.035 ms and gcc -O2's 0.055 ms,
and **245× faster than Python**. **Signing** is **≈ 1.7× faster than emulated C**
and **18.8× faster than Python**, and **≈ 7.4× slower than native C**. **Key
generation** is **296× faster than Python** but now **≈ 1.3× slower than even
emulated C** and **≈ 2.8× slower than native C** — the `BigInt` versus RNS gap
widens with degree, as expected, since `ntru_solve`'s big-integer work grows
faster than the parts Julia does well. At this degree `-O3` is the better Julia
level for keygen and sign (the opposite of n = 512), which is why the table
takes the best level per operation rather than fixing one.

### 9.3 Distributions

Key generation has a long right tail (`ntru_solve` retries), so the median is
the honest statistic and the p25/p75 columns of the raw files show the spread;
the quartiles are recorded for all three back ends and both degrees in
`docs/benchmarks.md`.

---

## 10. Constant time: explicitly out of scope

This implementation is **not** constant-time and does not attempt to be, and the
choice is deliberate in the sense argued in §1. Where the code does something a
constant-time implementation could not, a `CONSTANT TIME:` comment marks it —
`SamplerZ`'s data-dependent rejection loop, the secret-dependent `μ` that is a
plain `Float64`, the `Ref`-guarded spelling flags read on the hot path. A
security venue would require these closed; a Julia implementation cannot close
the implementation-level ones (GC, JIT, bounds checks, dispatch) at all, which
is precisely why it is a clean instrument for the specification-level question
this paper asks and a poor vehicle for a deployable signer. We state this rather
than claim a side-channel contribution we did not make.

---

## 11. Related work

**NTRU key generation.** The key generation implemented here follows Pornin and
Prest's field-norm tower solver (ePrint 2019/015), which reduces the time and
space of solving `fG − gF = q` by quasilinear factors and is what makes onboard
key generation feasible; our performance discussion of key generation (§9) is
against that algorithm as the reference realizes it (a 31-bit-limb RNS), and our
residual 1.9× gap to native C is the cost of `BigInt` over that RNS.

**Isochronous sampling.** Howe, Prest, Ricosset and Rossi (PQCrypto 2020) give
the isochronous `SamplerZ`/`BerExp` that the specification adopts; its purpose
is to make the sampler's *timing* independent of its secret inputs. It does not
make the sampler's *output* independent of floating-point spelling — that is a
different axis, and the one this paper and ePrint 2024/1709 measure.

**Side-channel leakage.** Fouque et al. (EUROCRYPT 2020, ePrint 2019/1180)
recover the key from Gram–Schmidt norm leakage; Karabulut and Aysu (DAC 2021)
and Guerreau et al. (TCHES 2022) mount power and timing attacks on the sampler.
These exploit *implementation* leakage under physical observation. This paper's
divergences need no side channel: two conforming implementations, given
identical inputs, simply return different bytes, which in derandomized mode is
already a key-recovery oracle (ePrint 2024/1709 §5).

**Floating-point error sensitivity.** ePrint 2024/1709 is the closest work and
is treated in full in §7. The one methodological neighbour is any study of
"specification vs reference rounding" as a subject in its own right; we are not
aware of one for a NIST PQC standard, and this is where the paper's independent
byte-exact reproduction (§4) does work the attack papers do not need to.

**Bit-exactness in sibling standards.** FIPS 204 (ML-DSA) contains an explicit
"No Floating-Point Arithmetic" provision; FIPS 205 (SLH-DSA) is hash-based. Both
therefore guarantee cross-implementation bit-exactness by construction, which is
the property FIPS 206 must engineer for FALCON rather than inherit, and the
reason §5's underdetermined choices are normative for it. Go 1.27 (August 2026)
shipped `crypto/mldsa` in its standard library, following ML-KEM in 1.24, while
no standard library ships FN-DSA — a concrete illustration that the
integer/bit-exact schemes are the ones deploying first.

### 11.1 The sensitivity is specific to Falcon's design (cross-scheme)

ePrint 2024/1709's footnote 6 states that Mitaka and Antrag, though they reuse
Falcon's `SamplerZ`, are not sensitive, because they call it with continuously
distributed centres. We tested this — and extended it to HAWK — by reading the
schemes' reference implementations directly (rather than porting, which risks
attributing a port bug to the scheme): `espitau/Mitaka-EC22`, `mti/antrag`, and
the KAT-validated `mjosaarinen/lil-hawk-py`.

| scheme | signing sampler | centre into the integer sampler | FP in signing | Falcon-class sensitive |
|:-------|:----------------|:--------------------------------|:-------------:|:----------------------:|
| Falcon | ffSampling tree (Klein–GPV) | rational, denominator `q` or `‖(g,−f)‖²` at four positions | yes | **yes** |
| Mitaka | hybrid (Peikert), continuous Box–Muller perturbation | continuous | yes | no |
| Antrag | hybrid, continuous perturbation (builds on Mitaka) | continuous | yes | no |
| HAWK | integer CDT, parity-bit centre | `{0, ½}` via a bit-selected fixed table | **no** (signing) | no |

Two facts from the code. Mitaka's `sampler` (and Antrag's, which is built on it)
adds a continuous Box–Muller normal (`normaldist`) to the target before the
discrete sampler, so its centre is continuous and integer centres occur with
probability zero — footnote 6 is correct by construction. HAWK's signing is
entirely integer: its `SamplerSign` centre is a single parity bit selecting one
of two fixed cumulative-distribution tables, and its only floating-point FFT is
in the *verifier's* `RebuildS0`, which does not affect the produced signature.

This yields a two-condition characterisation: the sensitivity requires **both**
(a) floating point in the signing sampler **and** (b) a centre that takes exact
rational values with small denominators (hence near-integer with non-negligible
probability). Falcon is the only one of the four deployed NTRU/lattice
hash-and-sign schemes with both; Mitaka and Antrag have (a) but not (b), HAWK
has neither. The sensitivity is therefore not a generic hazard of lattice
hash-and-sign signatures but a specific consequence of Falcon's design — which
is also why the fixes belong in the Falcon/FN-DSA specification rather than in a
general lattice-signature guideline.

---

## 12. Limitations

1. **A2a / A2b are not separated.** Complex division alone and `D11` alone each
   gave 0 / 100000; separating a ≈ 3.75 × 10⁻⁶ single-arm rate from zero would
   take millions of signatures, which we did not reach. So we cannot say whether
   one respelling suffices or both are needed — only that the two together reach
   the signature.
2. **A2 and A1 rates are not distinguished.** They differ by a point factor of
   2.5 but not significantly (p ≈ 0.23 over 925000 pooled signatures); both are
   order 10⁻⁵. We report this as "comparable," not as an ordering.
3. **The asymmetry explanation is a mechanism, not a calibrated model.** §6.3
   measures the perturbation magnitude `|Δμ|` per position (32× larger at the
   last two calls) and §6.3's follow-up converts it to a straddle probability
   `≈ |Δμ|/spread`. Those two give ratios of 32× and 5–10× respectively, and
   2024/1709's Table 2 30/70 split implies ≈3× after correcting for the density
   difference — three numbers spanning a factor of ten, which we call consistent
   in direction but have not reconciled quantitatively. The only direct test of
   the conditional probability at the first two calls is 0 of 10 manufactured
   integer centres, which cannot distinguish `P = 0.1` from `P = 0`. The
   mechanism (perturbation accumulating through the `l10` correction) we regard
   as established; the quantitative law we do not.
4. **Our measured rates are last-two-only.** `divergence_rate.jl` and
   `divergence_accum.jl` hold the message fixed, and §6.2 shows the first-two
   centres do not depend on the signing randomness, so across 525000 A1
   signatures there are only 100 independent first-two centres. The A1 rate we
   report is therefore not directly comparable to 2024/1709's Table 2, which
   includes first-two events; the comparison in §6 should be read with that
   caveat.
5. **A withdrawn claim.** An earlier version of this work stated that part 2 of
   the §7.1 countermeasure was undeployable, from 0 odd-norm keys in 3200. That
   was an argument from absence and it was wrong; §7.1(v) reports the
   measurement that refutes it. The corrected claim is narrower and is the one
   made here.
4. **One machine.** The performance spread is across compilers and optimization
   levels, not hardware.
5. **FALCON-1024** is measured at lower volume (350000 signatures, one event);
   the rate estimate there is loose (95 % CI 7 × 10⁻⁸ … 1.6 × 10⁻⁵), but the one
   event carries the full mechanism, so replication is established even though
   the rate is not pinned.

---

## 13. Implications for FIPS 206

If FIPS 206 requires bit-exact KAT agreement, then:

1. **Pin the spellings of §5.1.** Give complex division, the `D11` computation
   and the bottom-level `ffSampling` block as explicit operation sequences, or
   adopt ePrint 2024/1709 §7.1's structural fix (sample the centre with `round`
   rather than `floor`, and require `‖(g,−f)‖²` odd) — the latter removes the
   sensitivity itself rather than pinning each instance, and is the more robust
   choice. Note that adopting it also requires relaxing the reference's
   key-generation parity condition, which currently forces `‖(g,−f)‖²` even
   (§6.1).
2. **Publish `fpr_inv_sigma` at full precision** rather than leaving it to be
   derived from a rounded σ, and state whether intermediate values are within
   KAT scope.
3. **State σ_min's ε**, or its derivation.
4. **State the SamplerZ byte-assembly convention** alongside the vectors.

The draft NIST comment built from these points is in `docs/nist_comment_draft.md`
(to be posted when the FIPS 206 IPD appears).

---

## 14. Conclusion

An implementation written from the FALCON specification reproduces the reference
signatures byte for byte, and in doing so makes visible three places where the
specification does not determine the intermediate values. Two of them — complex
division and the `D11` spelling — reach the signature at ≈ 7.5 × 10⁻⁶ per
signature, at the position where a known attack turns one discrepant pair into
full key recovery, and are not removed by the countermeasure for the previously
known instance. The mechanism is the one ePrint 2024/1709 identified; the source
is not. The right fix is structural, and the specification should either adopt it
or pin the arithmetic. Along the way, an implementation in a language that cannot
hide its own runtime variance turns out to be an unexpectedly clean instrument
for seeing specification variance — and, incidentally, one that verifies FALCON
signatures faster than any C build we measured, while remaining a factor slower
at the floating-point-bound and big-integer-bound operations where a decade of C
optimization shows.

---

## Appendix A. The two implementation patterns, in full

The complete side-by-side of the two spellings is in the source:
`src/fft.jl` (`div_fft` / `_cdiv_cref`, complex division), `src/ffsampling.jl`
(`ldl_fft`, the `LDL_CREF` branch for `D11`; `_ffsampling_c4` and the generic
recursion for the bottom levels), and `src/Falcon.jl` (the exported
`with_spec_spelling`, `CDIV_CREF`, `LDL_CREF`, `FFSAMPLING_CREF`). The flags
default to the reference's spelling, so the byte-exact reproduction of §4 is what
runs unless a measurement asks otherwise.

## Appendix B. Reproducing the measurements

| result | command |
|:--|:--|
| byte-exact reproduction | `julia --project=falcon -e 'using Pkg; Pkg.test()'` (the "byte for byte" testset) |
| divergence rates (§6) | `julia --project=falcon falcon/scripts/divergence_rate.jl 100 1000` |
| one divergence, mechanism (§6) | `julia --project=falcon falcon/scripts/first_divergence.jl A2 70 534` |
| key recovery from an A2 pair (§6.1) | `julia --project=falcon falcon/scripts/key_recovery.jl 70 534` |
| per-position perturbation profile (§6.3) | `julia --project=falcon falcon/scripts/position_profile.jl 40 250 A1` |
| precision→rate law (§6.6) | `julia --project=falcon falcon/scripts/precision_law.jl 10 1200` |
| countermeasure, distribution (§7.1) | `julia --project=falcon falcon/scripts/countermeasure_eval.jl chisq` |
| countermeasure, rate (§7.1) | `julia --project=falcon falcon/scripts/countermeasure_eval.jl rate A1 100 1500 1 out.txt` |
| centre denominators (§7.1) | `julia --project=falcon falcon/scripts/denominator_check.jl 100 500` |
| **break the countermeasure** (§7.1) | `julia --project=falcon falcon/scripts/countermeasure_break.jl 100 1600` |
| **the replacement countermeasure** (§7.2) | `julia --project=falcon falcon/scripts/snapping.jl` |
| centres as a linear form in f (§6.4) | `julia --project=falcon falcon/scripts/linear_form.jl` |
| rounding mode changes every signature (§7.3) | `julia --project=falcon falcon/scripts/rounding_mode.jl 15 120` |
| the attack-family ceiling (§7.3) | `julia --project=falcon falcon/scripts/rounding_attack.jl 20 2400 1022` |
| Heuristic 1 (§6.1) | `julia --project=falcon falcon/scripts/heuristic1_check.jl 100 1000` |
| `fpr_inv_sigma` audit (§8) | `julia --project=falcon falcon/scripts/inv_sigma_audit.jl` |
| performance (§9) | `sh falcon/scripts/bench_all.sh` |
