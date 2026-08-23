# NIST への意見（下書き）― FIPS 206 の IPD 公開時に投稿する用

**まだ投げない。** FIPS 206 の initial public draft は 2026-08-21 時点で
未発行である（`docs/refs.md`）。IPD が出た時点で、
本文の節番号を確認したうえで pqc-forum に投稿する。

宛先: `pqc-forum@list.nist.gov`（または IPD が指定するコメント窓口）

---

## 投稿前チェックリスト

- [ ] IPD が実際に公開されたか（csrc.nist.gov の FIPS 一覧）
- [ ] 下記の指摘が IPD で**既に解決されていないか**を本文で確認
- [ ] 節番号・アルゴリズム番号を IPD の実際の番号に差し替え
- [ ] `docs/related_work_prompt.md` の調査を済ませ、既知の指摘でないことを確認
- [ ] 再現用リポジトリを公開状態にする（現在は private branch）
- [ ] §1 の縮約の結果を、**投稿環境とは別のマシン／別のコンパイラ版**でもう一度
      走らせる（現在 x86-64・gcc 13.3 / clang 18.1 の 1 点のみ。aarch64 では
      FMA が基本 ISA なので挙動が変わりうる ― 変わるなら §1 はもっと強くなる）

---

## 本文（英語・草案）

Subject: FIPS 206 (FN-DSA) — a build flag that changes the signature and gives
up the private key, two formulations the specification leaves open, and two
under-specified constants

Dear NIST PQC team,

We have written an independent implementation of Falcon (round-3
specification v1.2, 01/10/2020) in Julia, from the specification, and then
made it reproduce the C reference implementation's signatures bit for bit.
Since the FIPS 206 status update (Perlner, September 2025) indicates that the
standard will fix the order of floating-point operations, forbid fused
multiply-add, and require implementations to match KATs exactly, we think the
following observations are relevant to the draft. All of them are reproducible
from the artifact linked at the end.

**Summary.** The reference implementation itself, compiled two conforming ways,
produces different signatures on about one message in ten thousand, and two
thirds of those disagreements give up the private key. The difference is one
optimization flag that C99 permits and the specification does not mention.
Separately, two places where the specification's formulas and the reference's
code are algebraically identical also reach the signature, at about 7.5e-6, and
are not removed by the countermeasure of ePrint 2024/1709 section 7.2. Two
further items concern constants the specification does not let a reader derive.
If the standard requires bit-exact agreement, all of these need to be pinned in
the text — the first one normatively and first.

### 1. `FP_CONTRACT` is not set, and it changes the signature

C99 6.5p8 permits an implementation to contract `a*b + c` into a single fused
multiply-add, with one rounding instead of two; `FP_CONTRACT` governs it. The
reference sets neither the pragma nor a compiler flag, and neither the round-3
specification nor the reference's `README.txt` mentions contraction. The
September 2025 status update indicates FIPS 206 will forbid fused multiply-add;
we are writing to say that the prohibition has to be normative and explicit,
because of what we measured when it is absent.

We built the reference three ways (`-DFALCON_FPNATIVE=1 -march=native`,
everything else equal), generated one key from one seed, and signed 100000
messages under a per-message deterministic PRNG tape:

    clang -O2 -march=native -ffp-contract=fast : 143 fma instructions
    clang -O2 -march=native                    :   0
    gcc   -O2 -march=native                    :   0

    gcc default == clang default : identical over 100000 signatures
    clang -ffp-contract=fast vs clang default  : 12 of 100000 differ (1.2e-4)

All three builds produce the identical private key from the identical seed, so
the comparison is about signing alone. We then handed each of the twelve
divergent pairs to the key recovery of ePrint 2024/1709 section 5.1, which sees
only the two signatures, the public key and the message:

    8 of 12 : private key recovered, each an exact NTRU symmetry +-x^k (f,g)
              of the generated key, with ||(g,-f)||^2 = 16676 <= 1.17^2 q
    4 of 12 : not recoverable by that route -- these are divergences at the
              FIRST two sampler calls (see section 2 below)

We then checked that the recovered keys are usable rather than merely valid: for
each, we completed the trapdoor basis by solving fG - gF = q (the same descent
key generation runs), signed a message the original signer had never signed, and
offered it to the original public key. All eight were accepted. This is
universal forgery, and its cost follows from the rate: a key-yielding divergence
occurs on 8 of 100000 messages, so an adversary able to obtain both builds'
signatures on the same messages needs on the order of 1.3e4 messages.

The result replicates at FALCON-1024 (6 of 60000).

**The negative control matters as much.** Ordinary build choices do not
diverge. Over 30000 signatures each, these nine are byte-identical to one
another: `gcc -O2` native, `gcc -O2 -march=native`, `gcc` with the reference's
hand-vectorised AVX2 code path, `gcc -O2` emulated floating point,
`gcc -march=native` emulated, `clang -O2` native, `clang` AVX2, `clang -O2`
emulated, and `gcc` AVX2 with `-ffp-contract=fast` (GCC does not contract the
intrinsics either, emitting zero fma instructions there). Only clang with
`-ffp-contract=fast` differs — and contracting the *vector* path rather than the
scalar one gives a third, separate stream, so this is not "the reference and one
broken build" but three mutually inconsistent conforming builds, with the
divergent message set depending on which code path was contracted. Two compilers, native against emulated floating
point, and the reference's separate AVX2 implementation of the FFT all agree
bit for bit — so the hazard is one flag, not build variation in general. (We
note in passing that the emulated build agreeing with the native one is itself
worth recording, since `config.h` recommends emulation precisely because native
FPUs "may yield slight discrepancies that could affect determinism".)

**What is protecting the reference today is an accident.** The native build
wraps `double` in `typedef struct { double v; } fpr;`, and `fpr.h` states the
reason: so the compiler complains if raw arithmetic operators are used on the
type. It has a second effect that is documented nowhere. GCC will not contract
through the wrapper at any setting we tried, including its default
`-ffp-contract=fast` and `-Ofast`; clang will, but only at
`-ffp-contract=fast`, which `-Ofast` and `-ffast-math` both imply. A
reimplementation that does not use the wrapper — which is most of them, since
the wrapper is a C type-safety device with no analogue in many languages —
loses that protection without any way to know it had it.

**The reference offers the switch and under-rates it.** Beyond the compiler
flag, the reference ships an explicit `FALCON_FMA` build option (AVX2 path), and
`config.h` (lines 124-134) assesses its risk verbatim as: signatures "might
theoretically change, but only with low probability, less than 2^(-40); produced
signatures are still safe and interoperable." We built the reference with and
without `FALCON_FMA=1` on AVX2+FMA hardware and measured the signature-change
rate at **1.09e-4 = 2^-13.2** over 274944 messages -- about 2^27 (eight decimal
orders) higher than the stated bound -- and the divergent pairs recover the key
and forge. "Safe and interoperable" holds for a signature in isolation; it fails
when both builds' signatures on one message are observed, because their
difference is the key.

**No shipped build enables contraction today** (we surveyed liboqs, PQClean, the
Rust bindings, Bouncy Castle, the distributions, and the libraries that do not
yet ship Falcon at all). The deployed ecosystem is protected by integer
emulation, by hand-written two-rounding AVX2 intrinsics, and by the reference's
own warning -- none of which is in the specification. The threat is therefore
latent, and now is the time to make the prohibition normative, before FIPS 206's
move toward native-FP, bit-exact implementations makes it active.

**Suggestion.** State normatively that contraction must be disabled (`#pragma
STDC FP_CONTRACT OFF`, or an equivalent prohibition expressed over the affected
expressions), and that FMA must not be enabled in signing (the reference's
`FALCON_FMA` included), rather than leaving either to a build configuration. A
prohibition on "using FMA" that an implementer reads as "do not call `fma()`"
does not cover the compiler-contraction case, in which nobody wrote `fma`
anywhere; and a build option whose own documentation calls it "safe and
interoperable" will be enabled unless the standard forbids it.

### 2. The first two sampler calls are not harmless

Four of the twelve divergent pairs above are not recoverable by section 5.1 of
ePrint 2024/1709, because their difference is not 2-sparse in `z0`. They are
divergences at the *first* two sampler calls, which that paper's section 5 sets
aside: the resulting difference vector "is a short lattice vector, but it is not
expected to be short enough to make key recovery feasible."

That reasoning is about the difference vector, and there is a route at those
positions that does not read it. In `sample_preimage` the second target
component is `t1 = (-c_hat * b)/q` with `b = -f`, so as a ring element
`t1 = (c*f)/q` exactly, where `c` is the public hashed message. The first two
sampler centres are therefore coefficients of `(c*f)/q`, which we verified
against exact `Int128` arithmetic. An integer centre at the first call is
exactly the condition

    (c*f)[n/2 - 1] = 0   (mod q)

— one `F_q`-linear equation on the secret, with coefficients the adversary
computes from the public message. We collected `n-1 = 511` such events and
recovered the private key by Gaussian elimination over `F_q`: no difference
vectors, no lattice reduction, 49 seconds. The kernel is one-dimensional and
`||f|| << q`, so the lift is a scan over `q-1` scalars.

**Suggestion.** This is not a request for a change to the algorithm; it is a
request that the draft not justify a partial countermeasure by the harmlessness
of the first two positions. In particular, if FIPS 206 adopts the rounding
sampler of ePrint 2024/1709 section 7.1, it should mandate *both* parts of that
countermeasure in the same clause: with `round` in place of `floor` the
sensitive centres become the half-integers, and since `q` is odd a first-two
centre `N/q` can never be a half-integer, so part 1 alone closes the first two
calls unconditionally — while leaving the last two, which are the cheapest
key-recovery positions, exactly as exposed as `floor` did.

### 3. Two algebraically equal formulations that the text does not choose between

Two places where the specification's formulas and the reference's code are
algebraically identical produce different intermediate values, and those
differences reach the signature.

**(a) Complex division.** The reference computes the reciprocal explicitly:

    m = 1 / (b_re^2 + b_im^2);  b' = (b_re*m, -b_im*m);  d = a * b'

(`FPC_DIV` in `fft.c`). The specification gives no formula, and an
implementation that uses its language's complex division — Python's, or
Julia's, both of which use Smith's algorithm — gets a different last bit.
Smith's algorithm is the numerically better choice; the reference's is the one
the KATs encode.

**(b) The `D11` entry of LDL\*.** Algorithm 8 of the specification, and the
Python reference (`ffsampling.py`), compute

    L10 = G10 / G00;  D11 = G11 - L10 * adj(L10) * G00

while the C reference (`Zf(poly_LDL_fft)`) computes

    mu  = G01 / G00;  D11 = G11 - mu * adj(G01);  L10 = adj(mu)

Three multiplications against one. Equal in exact arithmetic because `G00` is
Hermitian and therefore real.

**These change the signature.** Over 400000 signatures at n = 512 with
identical key, message and PRNG state, and with every other computation
including the hand-unrolled bottom levels held at the reference's spelling,
three diverged — a rate of 7.5e-6 (Poisson 95%: 1.5e-6 to 2.2e-5). When one
diverges, roughly 470 of its 512 coefficients differ, because the byte stream
desynchronises.

All three are the mechanism of ePrint 2024/1709 (Lin, Tibouchi, Yu, Zhang,
EUROCRYPT 2025), Lemma 1:

    key 70 sig  534   call 1023 of 1024   mu -267.00000000000006  vs -267
    key 65 sig 1239   call 1023 of 1024   mu  337                 vs  336.99999999999994
    key 68 sig 1885   call 1024 of 1024   mu  285.99999999999989  vs  286

`SamplerZ` begins `s = floor(mu)`, and `floor` is discontinuous, so a 1e-13
discrepancy becomes a difference of 1 in `s`. In all three the leaf standard
deviation at that node is bit-identical between the two runs, so the
divergence is carried entirely by the centre. All three are at call 2n-1 or
2n, which is where section 5 of that paper recovers the entire private key
from a single discrepant pair.

**Why this is not already covered.** That paper's section 6.1 identifies a
different instance of the same mechanism — the reference's `sign_dyn` and
`sign_tree` modes order the bottom two levels of the tree traversal
differently — and its section 7.2 gives a countermeasure: make `sign_tree`
follow `sign_dyn`'s ordering. We reproduce that result as a positive control
(5 in 100000, against the ~3e-5 of its Table 2). But **(a) and (b) above are
present in both signing modes**, so aligning the two modes with each other
does not remove them. They are a difference between the specification and the
reference, not between two entry points of the reference.

**Suggestion.** If FIPS 206 requires bit-exact KAT agreement, the text must
give these two computations explicitly. We had thought a second option was
available — to say that any algebraically equivalent formulation is acceptable
— and the measurement removes it: in the presence of the sampler's `floor`,
algebraically equivalent formulations are not interchangeable.

We would also suggest adopting the countermeasure of that paper's section 7.1
(sample the centre with `round` rather than `floor`, and require
`||(g,-f)||^2` odd), which removes the sensitivity itself rather than any one
instance of it — subject to the caveat in section 2 above that both of its
parts must be mandated together. We reproduce its observation that the C
reference generates only keys with `||(g,-f)||^2` even (all 100 of ours are),
and we can confirm that this is a fixable property rather than an obstacle: the
reason is that `gen_poly_cdt` forces the coefficient sums of *both* `f` and `g`
odd, so `||(g,-f)||^2 = f(1) + g(1) = 0 (mod 2)`; requiring only one of the two
to be odd leaves the Pornin-Prest descent working (40 of 40 solver successes)
and every resulting key has an odd norm (0 of 40 in the control). The
key-generation change is therefore a normative change to the standard's key
generation, not an implementation note, and it needs to be written down.

### 4. `sigma` and `sigma_min` come from different `epsilon`

Table 3.3 of the round-3 specification lists both. Section 2.6 gives
`sigma` by equation (2.13) with `eps <= 1/sqrt(Q_s * lambda)`, `Q_s = 2^64`,
`lambda = 128` (Level I) or `256` (Level V). Evaluating that reproduces the
printed `sigma` to a relative error of 2e-16.

`sigma_min` is printed as a value only, with no derivation. Inverting the
smoothing parameter of Z,

    eta_eps(Z) = (1/pi) * sqrt( (1/2) * ln(2*(1 + 1/eps)) )

recovers `eps = 1/sqrt(2^64 * n^3)`, which reproduces both listed values to
3e-13. That is *not* the `epsilon` of (2.13): evaluating the smoothing
parameter at (2.13)'s `epsilon` is 11% off, and evaluating (2.13) at this one
is 10% off.

**Suggestion.** State the `epsilon` used for `sigma_min`, or give its
derivation. A reader who assumes the two share an `epsilon` — which the text
does not warn against — gets a value that is 10% wrong and that no test will
catch, because `sigma_min` only enters as a lower bound.

### 5. `fpr_inv_sigma` is not the correctly rounded reciprocal of `sigma`

The reference stores `1/sigma` per degree in `fpr_inv_sigma[]` (`fpr.h`) and
multiplies the ffLDL tree leaves by it (`ffLDL_binary_normalize` in `sign.c`:
"we actually store in the tree leaf the inverse of the value mandated by the
specification"). We could not reproduce those entries from any derivation the
specification supports. Computed at 400 bits and compared on raw bit patterns:

| derivation of `1/sigma`                          | logn = 9 | logn = 10 |
|:-------------------------------------------------|---------:|----------:|
| correctly rounded reciprocal of Table 3.3's printed sigma | +941 ulp | +18545 ulp |
| correctly rounded reciprocal of the Float64 sigma the Python reference carries | +1 ulp | +2 ulp |
| Float64 division `1.0 / sigma`                    |   +1 ulp |    +2 ulp |
| correctly rounded reciprocal of (2.13) evaluated at 400 bits | +1 ulp | +1 ulp |
| (2.13) evaluated in Float64, then `1.0 /` it      |   +2 ulp |    +1 ulp |

The table is above every one of them. The 35-digit decimal in `fpr.h` is also
not the 35-digit expansion of `1/sigma`: it agrees to about sixteen
significant figures and then diverges, so it appears to be some Float64
printed at length rather than a high-precision constant.

**How much this matters.** An implementation that computes the reciprocal
instead of transcribing the table builds a different expanded private key --
452 of 512 tree leaves differ at `logn = 9`. It does not build different
signatures, and this is settled by theory rather than by our measurement:
Lemma 2 of ePrint 2024/1709 bounds the probability of a divergent execution by
160*eps for a perturbation eps of the standard deviation, which at eps ~ 2^-52
is about 1e-11 per signature. Our 100000-signature run observed none, but its
resolution is 3e-5, six orders of magnitude too coarse to have seen anything
either way.

The contrast with section 3 is still the useful part, and it is the contrast
between those two lemmas: this constant perturbs the sampler's *width*, which
enters `BerExp`'s comparison smoothly, while section 3's differences perturb
its *centre*, which passes through `floor`. Whether an underdetermined choice
reaches the signature depends on which quantity it reaches, not on how large
it is.

**Suggestion.** If FIPS 206 keeps a precomputed reciprocal, publish the
constant itself at full precision rather than leaving it to be derived from a
rounded `sigma`, and say whether intermediate values are within the scope of
the KAT requirement.

**Related.** `inner.h` line 749 documents `fpr_sigma_min[]` as `1/sigma_min`.
The array holds `sigma_min` itself (1.2778... at `logn = 9`, which is Table
3.3's value, not its reciprocal). The code uses it correctly; only the comment
is wrong.

### 6. The SamplerZ test vectors' byte-consumption convention

Table 3.2 of the round-3 specification gives sixteen `SamplerZ` test vectors
as `randombytes` hex strings. Section 3.9.3 explains which bytes go to which
step ("the first 9 random bytes are used by BaseSampler, the next one by line
5 and the last one(s) by BerExp") but not the order in which those nine are
assembled into the 72-bit word. Read left to right, four of the sixteen
reproduce; with each chunk reversed, all sixteen do. The convention appears
only in the reference test harness (`KAT_randbytes` in `test.py`, which
applies `bytes.fromhex(oc)[::-1]`).

**Suggestion.** State the convention alongside the vectors.

### Artifact

[URL of the public repository]

- `scripts/cmp_cref_fp.jl` — the stage-by-stage bit comparison
- `scripts/cref_shim.c` — exposes the reference's sampler and ffSampling with
  a directly seeded PRNG, which is what made the comparison possible
- `test/vectors/cref_sign_kat.jl` — the reference's own signing output,
  recorded so the reproduction is checkable without a C compiler
- `scripts/cref_contract.sh` — section 1 end to end: the contraction probe, the
  three builds, the 100000-signature comparison, the eight-configuration
  negative control, the key recovery, and the FALCON-1024 replication, in one
  script (2 minutes 41 seconds)
- `scripts/cref_contract_recover.jl` — the section 5.1 recovery run on the real
  cross-build pairs, from the two signatures and the public key alone
- `scripts/event_solve.jl` — section 2: the key recovered from first-two events
  by Gaussian elimination over `F_q`
- `scripts/window.jl` — the two-sided condition under which an implementation
  difference is an oracle for those events, with both ends measured
- `scripts/divergence_rate.jl` — the 100000-signature measurement of both
  section 3 and section 5, run at the same sample size
- `scripts/first_divergence.jl` — for each divergent pair, the first sampler
  call at which the two disagree, with both centres printed
- `scripts/inv_sigma_audit.jl` — the ulp table in section 3
- `docs/debug_log.md` #046–#071 — the
  measurements behind each claim, including the hypotheses we tested and
  rejected. #054 records that we asserted the opposite of section 1's
  conclusion three times on a 480-signature sample before measuring at a
  sample size that could see it; #055 records that our first version of
  section 1 was a rediscovery of ePrint 2024/1709 section 6.1, because a
  single toggle we had not read carefully was moving only one of the three
  formulations we had named. #071 records that our first statement of section 1
  above — that GCC and clang *defaults* diverge — was wrong, and that building
  the reference rather than modelling it is what produced the correct and
  narrower claim.

We would be glad to supply anything further that is useful.

---

## 日本語メモ（投稿時には削る）

- **§1（FP_CONTRACT）が最も強く、最も具体的で、最も直しやすい。**
  「参照実装そのものが、適合する 2 通りのビルドで違う署名を出し、
  12 件中 8 件で鍵が出る」は、仕様書の書き方の話ではなく**実測された事故**である。
  IPD が FMA を禁ずるにしても、「`fma()` を呼ぶな」ではなく
  **「縮約を無効にせよ」**でなければ今回の事例を覆えない ―
  誰も `fma` と書いていないのだから。ここを強調すること。
- **陰性対照を必ず一緒に書く。** 8 構成がバイト一致することを書かずに
  「ビルドで署名が変わる」とだけ言うと、警告として強すぎて行動に移せない。
  「旗ひとつ」に絞れるのは対照があるからである。
- **§2 は「アルゴリズムを変えろ」ではない。**「最初の2回は無害」という
  <em>根拠</em>に寄りかかった対策を書かないでほしい、という要請にとどめる。
  こちらのオラクルは実装差であって、plain Falcon への実働攻撃は主張しない。
- **偶数ノルムの件は #067 で撤回済み。** 「配備できない」と書かないこと。
  正しくは「片方のパリティを外せば 40/40 で解け、全部奇数ノルムになる」＝
  **鍵生成の規範的変更として書き下す必要がある**。
- 論点は 2024/1709（Sleeping Falcon）と**同じ機構**である。
  違うのは摂動源で、あちらは **FMA** と **同一実装の 2 入口**
  （`sign_dyn`/`sign_tree`）、こちらは**仕様書と参照実装の差**。
  「違うビルドは要らない」は**誤り**だったので書かないこと
  （`sign_dyn`/`sign_tree` はビルドの差ではない）。
  言えるのは「**§7.2 の対策では塞がらない**」の一点。
  §1 の「Why this is not already covered」がそれ。
- §1 と §3 の**対照**は書くが、**発見として書かない**。
  中心（`floor` を通る）は Lemma 1、幅（`BerExp` の比較を通る）は Lemma 2 で、
  どちらも 2024/1709 が証明済み。こちらの数値はその追認にすぎない。
  §1 で新しいのは **(a)(b) が両 signing mode に共通で、§7.2 で塞がらない**点だけ。
- 3 番（`fpr_inv_sigma`）が最も具体的で、最も直しやすく、
  最も「言われないと気づかない」。かつ §1 の対照群として効く。
- 4 番は round-3 の話なので、FIPS 206 が KAT 形式を変えるなら
  自動的に解消する。IPD を見てから残すか決める。
- 分量は長すぎる。IPD を見て、既に解決している項を削れば適量になるはず。
