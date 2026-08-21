# NIST への意見（下書き）― FIPS 206 の IPD 公開時に投稿する用

**まだ投げない。** FIPS 206 の initial public draft は 2026-08-21 時点で
未発行である（`docs/refs.md`）。IPD が出た時点で、
本文の節番号を確認したうえで pqc-forum に投稿する。

> ## ⚠ 2026-08-21 夕 ― §1 と §3 は大幅に後退させる必要がある（#055）
>
> ePrint 2024/1709 の原典を読んだ結果:
>
> - **§1 の (c)**（ffSampling 最下段の手展開）は同論文 §6.1 そのもので、
>   §7.2 に対策がある。**投稿から落とすか、既知として引用に留める。**
>   そして 10 万署名で測った 5 件の発散は、この (c) によるものだった。
> - **§1 の (a)(b)**（複素除算・`D11`）は同論文に出てこず、
>   まだ伝播するか測っていない。ここだけが残る。
> - **§3** の「署名は変わらない」は同論文 **Lemma 2** が定理として
>   述べており、しかも予測率 1e-11 に対しこちらの実験の分解能は 3e-5。
>   「測って確かめた」と書ける水準にない。
>
> A2 の測定結果を待って書き直すこと。それまで投稿しない。

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

**They do change the signature, at a measurable rate.** Over 100000
signatures at n = 512 with identical key, message and PRNG state, five
diverged -- a rate of 5e-5. When one diverges, roughly 470 of its 512
coefficients differ, because the byte stream desynchronises.

The mechanism is the one ePrint 2024/1709 identifies. All five divergences are
nearly-integer centres:

    call 1023   mu = 221.00000000000006   vs   220.99999999999997
    call 1024   mu = 436.99999999999994   vs   437.00000000000011
    call 1024   mu = 444.0                vs   443.99999999999977
    call 1023   mu = -180.99999999999997  vs  -181.0000000000002
    call 1023   mu =  60.000000000000007  vs   59.999999999999957

`SamplerZ` begins `s = floor(mu)`, and `floor` is discontinuous, so a 1e-13
discrepancy becomes a difference of 1 in `s`. All five are at call 1023 or
1024 of 1024 -- positions 2n-2 and 2n-1, which is where that paper predicts
the centres concentrate near integers.

What we would add to it: 2024/1709 obtains its discrepancies by changing the
*build* (FMA on or off, emulated against native against AVX2, dynamic against
tree mode) -- different compilations of the same source. **A different build
is not required.** Two implementations that both follow the specification
suffice, because the specification does not say which of two algebraically
equal formulas to use.

**Suggestion.** If FIPS 206 requires bit-exact KAT agreement, the text must
give these three computations explicitly. We had thought a second option was
available -- to say that any algebraically equivalent formulation is
acceptable -- and the measurement removed it: the formulations are not
interchangeable in the presence of the sampler's `floor`.

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
444 of 512 tree leaves differ at `logn = 9`. Over 100000 signatures, measured
the same way as section 1 above, it builds the *same signatures*: none of
51200000 coefficients changed (95% upper bound on the rate, 3e-5).

That contrast is the useful part. This constant perturbs the sampler's
*width*, which enters `BerExp`'s comparison smoothly; section 1's differences
perturb its *centre*, which passes through `floor`. Whether an
underdetermined choice reaches the signature depends on which quantity it
reaches, not on how large it is.

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
- `docs/debug_log.md` #046, #048, #050, #053, #054 — the measurements behind
  each claim, including the hypotheses we tested and rejected. #054 records
  that we asserted the opposite of section 1's conclusion three times on a
  480-signature sample before measuring at a sample size that could see it.

We would be glad to supply anything further that is useful.

---

## 日本語メモ（投稿時には削る）

- 論点は 2024/1709（Sleeping Falcon）と**同じ機構**である。
  違うのは摂動源で、あちらは**ビルド構成**（FMA の有無、fpemu 対
  native 対 AVX2）、こちらは**仕様書に忠実な 2 つの実装**。
  「違うビルドは要らない」が言いたいこと。混同されないよう
  §1 の最後で明示的に差分を書いてある。
- §1 と §3 の**対照**が本体。中心（`floor` を通る）は 5e-5 で伝播し、
  幅（`BerExp` の比較を通る）は 0。「差の大きさ」ではなく
  「差がどの量に届くか」で決まる、という一文が結論。
- 3 番（`fpr_inv_sigma`）が最も具体的で、最も直しやすく、
  最も「言われないと気づかない」。かつ §1 の対照群として効く。
- 4 番は round-3 の話なので、FIPS 206 が KAT 形式を変えるなら
  自動的に解消する。IPD を見てから残すか決める。
- 分量は長すぎる。IPD を見て、既に解決している項を削れば適量になるはず。
