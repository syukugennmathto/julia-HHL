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

Subject: FIPS 206 (FN-DSA) — three places where the round-3 specification and
the reference implementation are algebraically equal but round differently

Dear NIST PQC team,

We have written an independent implementation of Falcon (round-3
specification v1.2, 01/10/2020) in Julia, from the specification, and then
made it reproduce the C reference implementation's signatures bit for bit.
Since the FIPS 206 status update (Perlner, September 2025) indicates that the
standard will fix the order of floating-point operations, forbid fused
multiply-add, and require implementations to match KATs exactly, we think the
following observations are relevant to the draft. All of them are reproducible
from the artifact linked at the end.

**Summary.** Reproducing the reference's signatures did not require matching
its floating-point arithmetic; it required access to its sampler's
randomness. But three places where the specification's formulas and the
reference's code are algebraically identical do produce different intermediate
values, and if the standard requires bit-exact agreement, those places need to
be pinned in the text.

### 1. Three algebraically equal, differently rounded formulations

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

**(c) The bottom two levels of ffSampling.** The reference hand-unrolls
`logn == 2` and folds the split/merge twiddles into the constants `1/sqrt(2)`
and `1/sqrt(8)`, doing one multiplication where the generic routine does two.

Our measurements: with every other floating-point primitive already agreeing
bit for bit (`FFT`, `iFFT`, `poly_split_fft`, `poly_merge_fft`, `poly_add`,
`poly_sub`, `poly_mul_fft`, `poly_adj_fft`, `poly_muladj_fft`,
`poly_mulselfadj_fft`, and the whole `B0` matrix of `expand_privkey`), these
three are what remain. Respelling them makes the expanded private key agree in
all 5120 doubles at `logn = 9`.

**They do not, however, change the signature.** Over 480 signatures at
n = 512 with identical key, message and PRNG state, the two spellings did not
differ in a single one of 245760 coefficients. The reason is margin:
`BerExp` compares a fixed-point exponential against random bytes, and an ulp
of slack in its argument flips that comparison with probability on the order
of 2^-52 per draw.

**Suggestion.** If FIPS 206 requires bit-exact KAT agreement, the text should
either give these three computations explicitly, or state that any
algebraically equivalent formulation is acceptable and that the KATs are not
to be read as pinning them. Either is fine; leaving it implicit means
implementers will discover the difference the way we did.

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
444 of 512 tree leaves differ at `logn = 9`. It builds the *same signatures*:
over 120 signatures with identical key, message and PRNG state, none of 61440
coefficients changed. So this affects a KAT on intermediate values or on the
expanded-key format, and not a KAT on signatures. We state that limit
explicitly because the opposite reading would be the alarming one and it is
not what we measured.

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
- `docs/debug_log.md` #046, #048, #050 — the measurements behind each claim,
  including the hypotheses we tested and rejected

We would be glad to supply anything further that is useful.

---

## 日本語メモ（投稿時には削る）

- 「浮動小数点が原因ではなかった」を**先に**書くのが要点。
  NIST 側も 2024/1709（Sleeping Falcon）を意識しているので、
  「丸めの差が署名を変える」という話と混同されると論点がぼける。
  こちらの主張は逆で、**bit は違うが署名は変わらない**である。
- 3 番（`fpr_inv_sigma`）が最も具体的で、最も直しやすく、
  最も「言われないと気づかない」。ここが本命かもしれない。
- 4 番は round-3 の話なので、FIPS 206 が KAT 形式を変えるなら
  自動的に解消する。IPD を見てから残すか決める。
- 分量は長すぎる。IPD を見て、既に解決している項を削れば適量になるはず。
