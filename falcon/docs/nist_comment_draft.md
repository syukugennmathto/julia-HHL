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

---

## 本文（英語・草案）

Subject: FIPS 206 (FN-DSA) — two formulations the round-3 specification leaves
open, which change the signature, plus two under-specified constants

Dear NIST PQC team,

We have written an independent implementation of Falcon (round-3
specification v1.2, 01/10/2020) in Julia, from the specification, and then
made it reproduce the C reference implementation's signatures bit for bit.
Since the FIPS 206 status update (Perlner, September 2025) indicates that the
standard will fix the order of floating-point operations, forbid fused
multiply-add, and require implementations to match KATs exactly, we think the
following observations are relevant to the draft. All of them are reproducible
from the artifact linked at the end.

**Summary.** Two places where the specification's formulas and the reference's
code are algebraically identical produce different intermediate values, and
those differences reach the signature at a rate of about 7.5e-6, at the
positions where ePrint 2024/1709 turns a discrepant pair into full key
recovery. They are not the instance of that mechanism which that paper's
section 6.1 identifies, and its section 7.2 countermeasure does not remove
them. Two further items concern constants the specification does not let a
reader derive. If the standard requires bit-exact agreement, all four need to
be pinned in the text.

### 1. Two algebraically equal formulations that the text does not choose between

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
instance of it. We note in passing that we reproduce its observation that the
C reference generates only keys with `||(g,-f)||^2` even, which blocks that
countermeasure until the key generation's parity condition is relaxed: all 100
of our independently generated keys are even.

### 2. `sigma` and `sigma_min` come from different `epsilon`

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

### 3. `fpr_inv_sigma` is not the correctly rounded reciprocal of `sigma`

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

The contrast with section 1 is still the useful part, and it is the contrast
between those two lemmas: this constant perturbs the sampler's *width*, which
enters `BerExp`'s comparison smoothly, while section 1's differences perturb
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

### 4. The SamplerZ test vectors' byte-consumption convention

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
- `scripts/divergence_rate.jl` — the 100000-signature measurement of both
  section 1 and section 3, run at the same sample size
- `scripts/first_divergence.jl` — for each divergent pair, the first sampler
  call at which the two disagree, with both centres printed
- `scripts/inv_sigma_audit.jl` — the ulp table in section 3
- `docs/debug_log.md` #046, #048, #050, #053, #054, #055, #056 — the
  measurements behind each claim, including the hypotheses we tested and
  rejected. #054 records that we asserted the opposite of section 1's
  conclusion three times on a 480-signature sample before measuring at a
  sample size that could see it; #055 records that our first version of
  section 1 was a rediscovery of ePrint 2024/1709 section 6.1, because a
  single toggle we had not read carefully was moving only one of the three
  formulations we had named.

We would be glad to supply anything further that is useful.

---

## 日本語メモ（投稿時には削る）

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
