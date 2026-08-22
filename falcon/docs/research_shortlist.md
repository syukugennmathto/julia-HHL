# Research shortlist (2026-08-22)

Produced by a 10-lens fan-out of 100 candidate directions, three independent
judge passes (novelty vs ePrint 2024/1709, feasibility on this machine, impact
if true), and a synthesis pass. The raw 80 first-round cards are in
`docs/research_cards.md`; this file is the ranked survivor list.

Status of directions already acted on is in `docs/debug_log.md` (#058-#067).

---

# FALCON / FN-DSA — Final Research Shortlist

Prepared from 100 candidate cards, three judge passes, and direct verification against `/home/user/julia-HHL/falcon` (paper §6.7, §7.1, §7.2, §7.3; `scripts/rounding_mode.jl`, `scripts/rounding_attack.jl`, `scripts/first_two_probe.jl`, `scripts/snapping.jl`), plus external anchors confirmed on disk (`/home/user/pornin/c-fn-dsa/kgen_ntru.c:1927`, `/home/user/algorand/falcon/deterministic.c`, `scripts/pyref/{encoding,falcon}.py`).

**Three facts that reshaped the ranking and must be respected by every card below:**

1. **§7.3 already owns the rounding-mode result** (1800/1800 under FE_UPWARD/DOWNWARD/TOWARDZERO, 0/1800 control). Any card whose novelty is "the FP environment is unpinned" is dead.
2. **§7.3 also established a ceiling**: key recovery via 2024/1709 §5.1 needs the *exact* centre at one of the last two calls to be an integer — probability `2/‖(g,−f)‖² ≈ 1.2×10⁻⁴` — and *no perturbation, however large, exceeds it*. Tail-window directed rounding gave 0 recoveries in 48 000 pairs. Any card promising "amplify the perturbation → amplify key recovery at the last two calls" is refuted in-house.
3. **§6.7 already publishes the first-two linear form** (`t1 = (c·f)/q` exactly; event ⇔ `(c·f)_{n/2−1} ≡ 0 mod q`; `n` events determine `f`) **and explicitly names the missing piece: an oracle for the event.** So the first-two direction is not new — *supplying and bounding the oracle is.* This distinction is the single most important novelty-risk line in the document.

---

## TOP 12

### 1. The oracle window: exactly which perturbations turn the first-two channel into a working attack

**Merged from cards 21, 22, 23, 24, 27, 30.**

**Hypothesis.** §6.7 establishes the channel and concedes it has no oracle. The missing piece is a *two-sided* condition, and it is a theorem, not a search. A perturbation of magnitude `|Δμ|` between two implementations is a usable event-oracle at calls 1,2 iff

```
   |η|  ≪  |Δμ|  ≪  (2/q) / 2n
   1.1e-13        1.6e-4 / 1024 = 1.6e-7
```

Lower bound: `η = μ̂ − (c·f)_{n/2−1}/q` is the *shared* rounding offset, measured at ≤ 1.4×10⁻⁹/q ≈ 1.1×10⁻¹³ (§7.1 ii). Below it, the two implementations sit on the same side of the integer — which is exactly why `first_two_probe.jl` got 0 straddles on 10 manufactured integer centres. Upper bound: the interior has **no atoms** (§6.4: 21 near-integer events in 10⁷ interior draws, exactly the uniform expectation), so interior straddles cost `2n·|Δμ|`, and once that exceeds `2/q` the boolean stops carrying information. Inside the window, `P(straddle | integer centre) ≈ 1/2` and `P(false positive) ≈ 2n·|Δμ|/(2/q) ≲ 10⁻³` — a near-noise-free oracle yielding one exact `F_q`-linear equation per event, and full key recovery by Gaussian elimination with **no lattice reduction and no signature values**.

The window also *explains and bounds* §7.3: directed rounding empirically gives rate 1, so it sits **above** the window (`|Δμ| ≳ 5×10⁻⁴`) — it destroys the oracle rather than creating it. That is the reason §7.3 measured 0 recoveries, and stating it converts a negative result into a design law.

**First experiment on this machine.**
(a) Free, five minutes: re-read the first-flip position histogram already emitted by `scripts/rounding_mode.jl` (it records `firstpos` bucketed into calls 1–2 / 2n−1..2n / elsewhere). If "elsewhere" dominates, `|Δμ|` for directed rounding is above the window — the window theorem's first empirical confirmation, at zero cost.
(b) Extend `scripts/precision_law.jl` to inject a **sign-biased** offset `+δ` at every centre, sweeping `δ ∈ {1e-14, 1e-13, 1e-12, 1e-10, 1e-8, 1e-6}`. For each `δ` measure, over the mod-q-preselected message set (below), `P(divergence | integer centre at call 1 or 2)`, and separately over 10⁵ generic messages the false-positive rate. Predict a plateau at ≈ 1/2 across `[1e-12, 1e-8]` with FP rate rising as `2n·δ`.
(c) Event generation with **no signing**: for one key, compute `(c·f)_{255}` and `(c·f)_{511} mod q` by NTT over ~4×10⁶ candidate messages (~minutes) — expect ~650 hits at `2/q`. Sign only those.
(d) Run the solve: 512 events → 512×512 Gaussian elimination mod 12289 (seconds) → centre-lift → check `h·f ≡ g`. Report the number of events at which the kernel becomes 1-dimensional (expect 511–515).
(e) Labeling (card 22): report `|Δμ|` at call 1 and call 2 **separately** (`position_profile.jl` currently pools them) — the two centres are `(c_re+c_im)·INVSQRT8` and `(c_im−c_re)·INVSQRT8` in `_ffsampling_c4`, so asymmetry is plausible and would remove the 2^k labeling ambiguity.

**Kill criterion.** Dead if (i) no `δ` produces `P(straddle|integer) ≥ 0.3` while keeping the false-positive fraction below 10⁻²  — i.e. the window is empty because `η` is much larger than measured; or (ii) the 512 collected equations do not give a 1-dimensional kernel whose centre-lift reproduces `h` with short `g`; or (iii) the events are unlabelable *and* no ambiguity-tolerant solver exists, in which case the honest output is an open problem, not an attack.

**Compute.** (a) free. (b) ~6 × 10⁵ signatures ≈ 25 min. (c) minutes, no signing. (d) seconds. (e) ~1.6 × 10⁶ signatures ≈ 60 min across arms. **Total: 2–3 CPU-hours.**

**What is most likely to go wrong.** *Novelty, not feasibility.* §6.7 already publishes the identity, the theorem, and the attack shape — a referee who reads the manuscript will say "you already had this." The genuinely new content is exactly three things: the two-sided window, the demonstration that directed rounding falls *outside* it (which reinterprets §7.3), and the executed solve. Frame it that way or it reads as self-plagiarism. Second risk: the window's realism. An *injected* biased δ is not an implementation; the card is only a full attack if some *realistic* difference lands in `[1e-12, 1e-8]` — the candidates are a radix-4/Stockham FFT, FMA in the FFT butterflies, x87 80-bit intermediates, and a fixed-point signer, none yet measured. Say so.

---

### 2. The r = 0 single-trace oracle: breaking *non-derandomized* FALCON

**Merged from cards 28, 64.**

**Hypothesis.** Every attack in this literature — 2024/1709's, this project's, and card 1 above — needs the same `(message, salt)` evaluated twice. The §6.7 equation does not: it needs only the boolean "call 1 (or 2) had an exact-integer centre," which is a property of **one** signature. When `μ` is exactly integral, `r = μ − ⌊μ⌋ = 0.0` bit-exactly, which is an extreme point of the rejection sampler. Two complementary claims: **(a)** timing is provably blind — by Poisson summation `ρ_σ(Z−r) = σ√(2π)(1 + O(e^{−2π²σ²}))` with `e^{−2π²σ_min²} ≈ 10⁻¹⁴`, and `ccs = σ_min/σ` cancels the width, so `E[iterations]` is independent of both `r` and `σ` to relative 10⁻¹⁴ (this *is* the isochrony property, and it doubles as a perfect mask); **(b)** a value-level probe is not blind — at `r = 0` with `(z0, b) = (0,0)` the exponent `x` is **bit-exactly 0.0**, so `trunc(x·2⁶³) = 0` and every Horner step of `approxexp` becomes a zero-operand multiply, the loudest possible Hamming-weight event, occurring with probability ≈ 0.18 per rejection round. A single-trace SPA/EM adversary who can locate the *first* sampler call gets one equation per signature **with a fresh salt** — breaking plain, unmodified, randomized FALCON in `n·(q/2) ≈ 3×10⁶` traces.

**First experiment.** Pure statistics, entirely on `samplerz`, no full signing. Over 10⁷ calls at the two leaf widths that calls 1 and 2 actually use, tabulate the joint distribution of (rejection-loop iterations, total bytes consumed, BerExp inner rounds, byte index at which BerExp exits) conditioned on `r = 0.0` exactly vs `r ~ U(0,1)`. Compute total-variation distance and optimal Bayes error → traces-per-equation → traces-to-key. Confirm `x === 0.0` bit-exactly in the `(0,0)` branch. Then run `scripts/dudect.jl` (existing controls, measured noise floor) on fixed-`r=0` vs random-`r` to confirm the timing null.

**Kill criterion.** Dead if TV distance between the `r=0` and generic-`r` distributions of (iterations, bytes, BerExp path) is below ~2⁻¹⁰ per call — a noiseless trace could not flag the event. Also dead if the advantage survives only under averaging across signatures, since each signature has a different `c` and the decision must be per-signature. If dudect *does* find timing above 3× the noise floor, that is the stronger result and the card flips.

**Compute.** 10⁷ SamplerZ calls ≈ 10–30 s. dudect ≈ minutes. **Under an hour.**

**What is most likely to go wrong.** The physical half cannot be run — no oscilloscope, no target board — so the honest deliverable is a **leakage-model quantification plus a separation theorem**, with the EM measurement marked future work. Reviewers at CHES will want the trace. Second: Falcon sampler side channels are a crowded field (Guerreau et al. TCHES 2022; Karabulut–Aysu; Fouque et al. EUROCRYPT 2020) — check that no one has already published an `r ≈ 0` distinguisher before claiming it. Third: the project has already recorded "timing oracle from BerExp byte consumption: judged not externally observable" as a negative; part (a) *converts that into a theorem*, which is good, but the write-up must not read as re-litigating it.

---

### 3. Adversarial spelling: a KAT-passing, spec-defensible, key-leaking implementation

**Card 43, constrained by the §7.3 ceiling and the card-1 window.**

**Hypothesis.** Spelling choice is an *attack parameter*, not an accident. The structural ceiling is `2/‖(g,−f)‖² ≈ 1.2×10⁻⁴`; A2 sits 16× below it only because `|Δμ| ≈ 10⁻¹⁴` is comparable to the shared rounding offset. An adversary who *chooses* the spellings drives `|Δμ|` into the same window as card 1 — `[10⁻¹², 10⁻⁸]` — where `P(straddle | integer centre)` saturates near 1/2 while the generic rate `2n·|Δμ| ≤ 10⁻⁵` stays invisible. The result is an implementation in which every individual choice is defensible (Smith's division, pairwise vs sequential FFT summation, a valid alternative twiddle ordering, `hypot`, compiler-contracted FMA) that leaks the whole key from one repeated-syndrome pair at ~10⁻⁴ and passes a 100-vector KAT with ≈ 99% probability. **An algorithm-substitution attack whose payload is deniable as a rounding preference.**

**First experiment.** (1) Use the cheap surrogate — mean `|Δμ|` at positions 2n−2, 2n−1 over 10⁴ signatures via `scripts/position_profile.jl` (~20 s per candidate) — as the objective, not end-to-end rates. (2) Build a menu of individually-defensible spellings per primitive in `src/fft.jl` / `src/ffsampling.jl` (≥5 complex-division algorithms, ≥3 summation orders, LDL variants, split/merge twiddle foldings, FMA on/off) behind the existing `with_spec_spelling` toggles; greedily maximise last-two `|Δμ|` subject to interior `|Δμ| < 10⁻¹⁰`. (3) Confirm the winner end-to-end over 2×10⁵ signatures and run `scripts/key_recovery.jl` on a resulting pair. (4) Plot KAT-catch probability vs leak rate. (5) State the defence: card 1's window makes the backdoor detectable by construction, and §7.2 snapping neutralises it outright.

**Kill criterion.** Dead if no combination of defensible spellings pushes last-two `|Δμ|` above ~3×10⁻¹³ (rate stays pinned near A2's 7.5×10⁻⁶, no amplification available), **or** if every amplifying combination also raises interior `|Δμ|` enough that `2n·|Δμ|` exceeds the structural rate — which would make it catchable by ordinary KATs and destroy deniability.

**Compute.** ~200 candidates × 10⁴ signatures ≈ 1 CPU-hour; winner confirmation 2×10⁵ ≈ 8 min. **~2 hours.**

**What is most likely to go wrong.** The ceiling is hard: the *best possible* leak is 1.2×10⁻⁴ per same-syndrome pair, i.e. ~10⁴ pairs to a key. That is a real attack but it is not "one pair"; do not overclaim. Second, ASA/kleptography is a well-populated literature (Bellare–Paterson–Rogaway; Young–Yung) — the novelty is the *deniable-numerics* payload, not the ASA framing. Third, the search may find nothing in the window, since the shared rounding offset is what suppresses A2 and it suppresses most respellings equally.

---

### 4. Bias-and-average: key recovery from ordinary, verifying, randomized signatures

**Card 62.**

**Hypothesis.** Adding an *exact integer* `a` to the centre of one bottom leaf leaves `r = μ − ⌊μ⌋` bit-identical, so the sampler consumes identical randomness and returns `z + a` deterministically. Every signature then becomes `s' = s + a·x^j·(g,−f)`, still satisfying `s1' + s2'h = c` and still under `β²` (for `a = 8`, ≈ +1.0×10⁶ against a mean `‖s‖² = 2.81×10⁷` and `β² = 3.40×10⁷`). Since `E[s] = 0` over random salts, the coefficientwise mean of `N` faulted signatures converges to `a·x^j·(g,−f)`; with per-coefficient noise `165.7/√N` needing to fall below `a/2`, `a = 8` gives `N ≈ 2×10³`. **This removes the same-syndrome requirement, the derandomization requirement, the chosen-message requirement, and the second implementation** — leaving one persistent fault and a few thousand ordinary signatures that all verify.

**First experiment.** Add a flag to `_ffsampling_c4` in `src/ffsampling.jl` adding constant `a` to the last leaf's `μ`. Sign one key 3000 times with fresh random salts (~7 s). Assert every signature verifies; record the norm-rejection rate against the unfaulted baseline. Average `s1`, `s2` coefficientwise, round, compare against `a·x^j·(g,−f)` for the stored key. Divide out with an extended `scripts/key_recovery.jl` search (allow `a·x^j` as well as `a + b·x^{n/2}`) and check the recovered `(f,g)` reproduces `h mod q`. Sweep `a ∈ {1,2,4,8,16}` and confirm `N ∝ a⁻²`. Then the realistic variant: scale one deepest-node `l10` entry by `(1+ε)` so `a` varies per signature, and confirm `E[a] ≠ 0` still suffices (only the direction is needed).

**Kill criterion.** Dead if faulted signatures fail the norm/compression check more than ~5% of the time (the fault is loud), **or** if the coefficientwise mean does not converge to an integer multiple of `(g,−f)` (a single-leaf offset does not map to a fixed short `dz0`), **or** if the smallest silent `a` still needs `N > 10⁷`.

**Compute.** ~3 s of signing per point; whole sweep **under 10 minutes.** Cheapest strong card on the board.

**What is most likely to go wrong.** *Novelty.* Falcon fault attacks are published (BEARZ; loop-abort faults on lattice signatures, Espitau et al.) and "average many faulted signatures to recover a fixed offset" is close to standard statistical fault analysis. **Check the fault-attack literature before writing a line.** The genuinely new observation is the `r`-invariance — that an *integer* offset leaves the rejection sampler's randomness consumption bit-identical, so the fault is perfectly silent and the offset is exactly fixed. If that specific point is unpublished, the card holds; if not, it collapses to an instantiation. Second risk: the fault model is *assumed*, not demonstrated — there is no fault-injection rig here.

---

### 5. Key generation is not a function of its seed: the acceptance-threshold band

**Merged from cards 77, 84, 79.**

**Hypothesis.** All prior work — 2024/1709 and this project — studies *signing*. Key generation contains a floating-point **rejection** decision, and a rejected candidate advances the RNG stream, so one flipped decision produces an **entirely unrelated key pair from the same seed**. This is not a 10⁻⁵ statistical leak; it is a total, visible interop break that FIPS 206's announced keyGen KATs cannot tolerate. And the disagreement band is not one ulp: verified on disk at `/home/user/pornin/c-fn-dsa/kgen_ntru.c:1927` —

> `/* Constant is (0.999*1.17*sqrt(q))^2 … (Standard says 0.9999, we are slightly more restrictive just to be sure that our keys are always compliant despite using fixed-point approximations for the orthogonalized norm.) */`

— a **deliberate 1.8×10⁻³ relative disagreement**, a thousand times larger than any FP effect anywhere in this literature, in the FN-DSA-tracking implementation by the reference's own author. Two further sites amplify: the `max|F|` encodability cliff against `max_fg_bits` (whose input is FP-reduced, and whose failure re-rolls the whole stream), and FIPS 206's *announced* new clause requiring ffLDL leaves in `[σ_min, σ_max]` — which converts this project's Class II constant (`fpr_inv_sigma`, 452 of 512 leaves differ, certified harmless by Lemma 2) into a key-generation rejection decision.

**First experiment.** No NTRUSolve needed. (1) Sample ~10⁶ candidate `(f,g)` with `gen_poly` and compute the orthogonalized norm in Float64 (BigFloat only for the near-boundary tail); histogram `bnorm/(1.17²q)` near 1; report the mass in `[0.999², 0.9999²]` (the c-fn-dsa vs standard band) and in a ±10⁻¹³ relative band (the Julia vs C-ref band). Restrict to candidates where the orthogonalized part is the max, since the `ℓ2` part is integer-exact and can never straddle. (2) Evaluate *both* thresholds on the same candidate stream — no c-fn-dsa build required. (3) Run Julia and C keygen on identical seeds (`scripts/cref_keygen_instr.c` exists) for as many seeds as the boundary density predicts, and exhibit one seed producing two different public keys. (4) Cheap one-run kill for the FIPS-206 leaf clause: min/max leaf σ over 10⁴ expanded keys (`leaf_sigmas` already exposed) — margin to `σ_min` in ulps.

**Kill criterion.** Dead if fewer than 1 in 10⁴ candidates lands within 2×10⁻³ relative of the bound (the 0.999 fudge is harmless) **and** fewer than 1 in 10⁸ within 10⁻¹³ (the two-spelling flip is unreachable). Also dead if every implementation performs this test in exact integer or fixed-point arithmetic — **check c-fn-dsa's `fxr` path, the round-3 `fpr` path and PQClean *first*.**

**Compute.** 10⁶ candidates ≈ 30 min in Float64; leaf-σ census ≈ 100 s of keygen; C/Julia seed diff ≈ 30 min. **~1.5 hours.**

**What is most likely to go wrong.** The 0.999 constant is *documented in a comment*, so a referee may call it known-and-intentional. The defensible claim is not "Pornin made a mistake" — it is "seed → key is therefore not a function across conforming implementations, which breaks seed-based derivation, HSM key backup, and keyGen KATs, and nobody has said so." Second: c-fn-dsa uses fixed point here, so this is a *threshold* disagreement, not an FP-nondeterminism one — do not conflate them. Third: `#064` already found 0/400 key divergences under the `div_fft` respelling, so the Babai-rounding sub-story is probably a null; lead with the threshold, not with Babai.

---

### 6. The ecosystem already disagrees: a third spelling, and a per-ISA one

**Merged from cards 76, 78, 50.**

**Hypothesis.** The underdetermination is not hypothetical for FIPS 206 — the implementations that will *define* its conformance already spell the sensitive operations differently. Verified: `pornin/c-fn-dsa` computes the LDL step as `inv_g00_re = 1/g00_re; mu_re = g01_re*inv_g00_re; mu_im = g01_im*inv_g00_re`, replicated across its generic/SSE2/NEON/RV64D backends, whereas the round-3 reference forms `1/(g00_re²+g00_im²)` then multiplies, and Algorithm 8 uses three multiplications. **Three values of the same quantity.** Separately, GCC defaults to `-ffp-contract=fast` on AArch64/RISC-V and c-fn-dsa's Makefile is plain `-O2`, so *the same source with the same flags gives different signatures per ISA* — which makes CMVP vendor-affirmation across operational environments a **correctness** claim, not a performance one.

**First experiment.** (1) Ten-minute analytic pre-check: `poly_mulselfadj_fft` stores an exact zero imaginary part, so FPC_DIV computes `g01_re·g00_re·(1/g00_re²)` while Pornin computes `g01_re·(1/g00_re)` — different roundings, so arm A3 is *not* a no-op. Confirm by bit-comparison on real Gram entries. (2) Add arm A3 to the toggle framework and measure its divergence rate against the C spelling over ~4×10⁵ signatures with `divergence_rate.jl`, plus the first-divergence position profile. (3) Cross-compile `sign_fpoly.c` and the reference `fft.c` with `clang --target=aarch64-linux-gnu` (verified working here) and grep the emitted assembly for `fmadd`/`fmsub`/`vfma` — this settles the per-ISA claim from the artifact, in minutes, with no ARM hardware. Reproduce the ARM numerics on x86 with Julia's hardware `fma()` at exactly the identified sites. (4) Static survey (grep for library complex division / `__divdc3`, the 3-mul D11, the unrolled `logn==2` block) across PQClean, liboqs, `pornin/rust-fn-dsa`, Bouncy Castle, pqm4 — all clone through the session git proxy; Go, Rust, Java and clang toolchains are installed.

**Kill criterion.** Dead for A3 if it is bit-identical to FPC_DIV on real inputs (kill in ten minutes, before any rate run). Dead for the ISA half if the AArch64 assembly contains no FMA for these expression shapes, or if the Rust and C builds are demonstrably bit-identical on ARM. Weakened if the survey shows every reachable implementation is a C transliteration.

**Compute.** A3 rate: 4×10⁵ signatures ≈ 15–30 min. Cross-compile + grep: minutes. Survey: 2–3 hours of reading. **~3 hours.**

**What is most likely to go wrong.** A3 is another instance of the established Class I(b), not a new mechanism — the contribution is "the ecosystem already disagrees," which is evidence, not theory. FMA divergence itself is 2024/1709 §6.2, so the ISA half must be framed as a *deployment/validation* argument, not a new perturbation. And the survey's likely outcome is "everything is a transliteration."

---

### 7. The exact-shadow certificate and the boundary inventory: a completeness theorem for §7.2

**Merged from cards 93, 94, 55, 32.**

**Hypothesis.** Two halves that need each other. **(a) The certificate:** run the signer alongside an exact shadow of the centre computation; at each of the `2n` leaves record the margin `m_k = dist(μ̂_k, Z)` and the forward error `e_k = |μ̂_k − μ_exact,k|`. The predicate `REACHABLE := min_k(m_k − C·e_k) ≤ 0` decides, **per signature**, whether *any* conforming respelling could change the output — quantifying over the whole equivalence class, not the two arms someone implemented. This turns a 10⁻⁵ rare-event measurement into a deterministic test with an unbiased low-variance estimator, and generalises verbatim to any FP cryptographic implementation with a discrete decision on an FP value. **(b) The inventory:** FALCON contains a small, enumerable set of FP→discrete boundaries — `fpr_floor(μ)` (`sign.c:1369`), `fpr_trunc(x·1/ln2)` in `fpr_expm`, BerExp's 64-bit byte comparison, `fpr_rint(t0)/(t1)` (`sign.c:854/872/1061/1079`), `fpr_rint(k)` in Babai (`keygen.c:3249`), `keygen.c:3682/3706/3707/3934`, and `fpr_lt(bnorm, fpr_bnorm_max)` (`keygen.c:4269`). Claim: **only `fpr_floor(μ)` is structural** (has an atom), because `rint(t0)/(t1)` evaluate exact integers with margin 1/2 always, BerExp's third input is Lipschitz with a measurable constant (`rate ≤ 2n·k·e`, an input neither Lemma 1 nor Lemma 2 covers), `1.17²q = 16822.47…` sits strictly between integers, and `exp(−x)·ccs` is never dyadic. **Corollary — the payoff:** §7.2's rational snapping is not a patch but **complete**: with an exact centre split, bit-exactness of FN-DSA is equivalent to nothing else.

**First experiment.** (1) Add a shadow mode to `src/ffsampling.jl` carrying a parallel `BigFloat(128)` copy of `t0/t1` down the tree, logging `(m_k, e_k)` per leaf. **Validate on the 13 divergent pairs already on disk** (A2 key 70/534, 65/1239, 68/1885; A1 4/105; the three `countermeasure_break` pairs; the n=1024 pair): sensitivity must be 1.0. Specificity: run over 10⁴ ordinary signatures and check the predicted rate covers the measured 1.9×10⁻⁵. (2) Histogram the margin at each inventory site over 10⁴ signatures / 3×10² keygens and test for an atom. (3) Constructive completeness test: with snapping **on**, re-run `precision_law.jl`'s uniform-ε injection — the law must change from `2n·ε` (measured 2.4e-7→4.2e-4, 9.8e-4→9.2e-1) to a much shallower slope, and that slope **is** the residual constant `C`. (4) Headline demo: degrade the whole FFT to `k` mantissa bits with snapping on and find the smallest `k` at which 1200 signatures stay bit-identical.

**Kill criterion.** Dead if the certificate fires on ≫10⁻⁴ of signatures (near-boundary margins common, crossings rare → certifies nothing) or fails to fire on even one known divergent pair (unsound). Dead as a theorem if any site other than `fpr_floor(μ)` shows an atom — in particular if `keygen.c:3249`'s `rint` argument turns out to be a small-denominator rational, which would be a *bigger* result than the theorem. Dead if, with snapping on, the injected-ε law does not flatten.

**Compute.** BigFloat shadow ≈ 15 min for 10⁴ signatures; margin histograms < 10 min; injection sweep ≈ 20 min. **~1 hour**, plus a day of code (the concrete `ComplexF64` annotations in `fft.jl`/`ffsampling.jl` must be loosened before a shadow type will dispatch).

**What is most likely to go wrong.** Shadow-value execution is standard numerical-debugging technique, so the novelty is entirely in the *completeness corollary*, not the instrument. And the soundness constant `C` needs a rigorous forward-error bound that §7.2 already flags as missing — without it the certificate is heuristic, and a referee will say so. Budget the derivation, not the code.

---

### 8. The malleability atlas: FN-DSA's verifier accept-set is underdetermined too

**Card 92.**

**Hypothesis.** Every result in this project is about the *signer's output*. The **verifier's accept-set** is also underdetermined — in integer arithmetic, with no floating point, and it is worse, because the attacker needs no key and no rare event. Verified on disk: `scripts/pyref/encoding.py`'s `decompress` strips trailing zero bits from the whole bitstring and then stops after the `n`-th coefficient, never checking the remainder; and `Falcon.verify` (`falcon.py:425–455`) checks **neither the header byte nor the signature length**. The C reference checks all three (`falcon.c:878–895`, `codec.c:465–472`). So for a verifier written from the specification's own reference, every FALCON signature has ~2^(8·slack) equivalent byte strings that all verify — **breaking SUF-CMA** and breaking any protocol treating a signature as a unique handle (transaction ids, Merkle-committed state proofs, signature-hashing aggregators). Independently: the C API with `sig_type = 0` accepts both `0x30` and `0x50` headers and `falcon_verify_start` hashes only `sig[1..40]`, so the header byte is unauthenticated and a COMPRESSED signature can be re-encoded into CT format and still verify. This is live standards business (the pqc-forum "which signature formats should FIPS 206 support" thread asserts canonicity is "enforced by the code").

**First experiment.** Build a mutation generator over the genuine C-reference signature in `test/vectors/cref_kat.jl`: (i) set each unused padding bit; (ii) randomise each padding byte; (iii) alter the header byte; (iv) re-encode `s2` in fixed-width CT format with a `0x50` header. Run every mutant through three verifiers already present: `src/falcon.jl`, the buildable `scripts/cref` (with `sig_type=0` and each explicit type), and `pyref`. Tabulate accept/reject. Measure malleability entropy over 10⁴ signatures: `slack = 8·(sig_bytes − 41) − used_bits`. Then clone PQClean, liboqs and `pornin/rust-fn-dsa` and run the same corpus for a cross-implementation compatibility matrix. (Note: pyref *does* reject `−0`, so drop that mutation class.)

**Kill criterion.** Dead if every third-party implementation rejects every mutant **and** the round-3 spec text (Algorithm Decompress, §3.11.2) mandates the trailing-bit, header and length checks — then this is a Python-reference bug report, not an underdetermination result. Weakened to incremental if only the Python reference is permissive and nothing deployed derives from it.

**Compute.** Minutes. The mutation corpus is tiny and there is no signing at volume. **Under an hour**, plus clone/build time.

**What is most likely to go wrong.** The C codec guards exactly these cases, which strongly suggests implementers knew — so parts may be folklore that was simply never written down for FN-DSA. Check the round-3 spec text *first*; if it mandates the checks, downgrade to a compatibility-matrix contribution and a standards comment.

---

### 9. Verifier-side floating point: the economics inversion

**Merged from cards 98, 47.**

**Hypothesis.** Everything above needs a 10⁻⁵ accident because the *signer's* inputs are honest. A scheme whose **verifier** uses floating point inverts the economics: the adversary supplies the verifier's input and can **grind** for a signature whose accept/reject decision sits within FP error of the boundary — no key, no repeated syndrome, no derandomization, just work. This project's established result 8 located exactly such a site and stopped there: HAWK's signing is integer-only but its verifier rebuilds `s0` by floating-point FFT (`RebuildS0`). If HAWK's accept decision compares an FP quantity rather than an exactly-recovered integer, two conforming verifiers differing only in summation order or FMA can be made to disagree on an attacker-chosen signature — a **non-repudiation failure and a consensus split**, the exact class that forced ZIP-215 for Ed25519, and strictly more damaging than a signer-side leak. It also **attacks this project's own two-condition characterisation**, which was derived only for signature *production*.

**First experiment.** Clone the HAWK C reference and `mjosaarinen/lil-hawk-py`. **Read first, grind second:** locate the accept decision downstream of `RebuildS0` and determine whether the compared quantity is (i) an exactly-recovered integer with a 1/2 margin — safe, the same structure as Falcon's `fpr_rint(t0)` — or (ii) a genuine FP comparison. If (ii): instrument the verifier to log `|value − bound|` over ~10⁴–10⁵ attacker-generated candidates (candidates need no key, so this is cheap), fit the margin density near the boundary, and extrapolate grinding work. Because the map from `s1` to `s0` is affine and the attacker chooses `s1` freely, **solve** for an `s1` driving one coefficient's margin below FP error rather than waiting. Then build two verifiers differing only in norm summation order, grind, and exhibit one byte string one accepts and the other rejects.

**Kill criterion.** Dead if the verdict is taken on exactly-recovered integers with a large margin (no boundary to grind) — in which case the two-condition characterisation simply extends to verification, which is a clean positive for the paper. Also dead if the margin density is so low that grinding exceeds 2⁶⁴.

**Compute.** Reading: hours. Instrumented verification of 10⁵ candidates: use the C reference, not the Python one (Python at 10⁶ is hours-to-days). **~1 day.**

**What is most likely to go wrong.** The likely finding is that `s0` is uniquely determined mod 2 by `h0` with a reconstruction error bounded far below 1/2 by design — i.e. the card dies, but as a *sharpening* of the characterisation into signing-side vs verification-side conditions, which is still worth a subsection. Building the HAWK C reference here is unverified.

---

### 10. How much of the tree is live: is "pin the operation order" even possible?

**Merged from cards 59, 49, 96, 15.**

**Hypothesis.** FIPS 206's announced strategy — fix the order of FP operations, forbid FMA, require bit-exact KATs — and **this project's own §13 recommendation** both assume the sensitive surface is small enough to pin textually. The claim is the opposite: a single 1-ulp change at essentially *any* of the 9728 doubles of the LDL tree propagates to `|Δμ|` at the last two calls exceeding the shared rounding floor, hence yields the full `≈1/‖(g,−f)‖²` rate — because every `z0`-branch centre passes through the top-level `l10` correction (the depth-accumulation mechanism §6.3 already established). If most of the tree is live, bit-exactness cannot be achieved by pinning a list of formulas, and only a structural fix (§7.2) works. **This converts a policy preference into a measurable theorem, and forces the project to change its own recommendation.**

**First experiment.** Extend `scripts/position_profile.jl` with a hook applying `nextafter` to exactly one tree double (chosen by index) per run, then measure `|Δμ|` at the last two calls over 10³ signatures. Sample ~200 of the 9728 indices stratified by tree level. Produce the per-level histogram of `|Δμ|` and the fraction of indices whose `|Δμ|` exceeds the measured shared-offset floor (`≈1.1×10⁻¹³`); convert to a predicted rate per index via the straddle law. Report the fraction of the 9728 doubles that are **live**. Then the completeness count: how many distinct arithmetic *sites* (not doubles) a standard would have to pin.

**Kill criterion.** Dead if a large majority (>70%) of single-ulp perturbations give last-two `|Δμ|` **below** the shared rounding floor — the surface is small and localisable, pinning is viable, and FIPS 206's plan is vindicated (a publishable outcome in its own right, but it kills the argument).

**Compute.** 200 indices × 10³ signatures = 2×10⁵ signatures ≈ **8 minutes.** Run card 1(a)/the straddle-law baseline first to fix the floor.

**What is most likely to go wrong.** The measurement is easy; the *inference* is the risk. Showing that many tree doubles are live does not by itself show that many *realistic spellings* move them — a referee will demand the bridge from "1-ulp anywhere" to "an implementer would actually write that." Pair it with the (cheap, mechanically enumerable) rewrite census and deduplicate variants that are bit-identical on random inputs, or the claim is unfalsifiable in the wrong direction.

---

### 11. Odd-norm keys: why they are unreachable, and attacking the half of §7.1 we conceded

**Merged from cards 81, 3, 82, 88.**

**Hypothesis.** Three linked claims. **(a) Mechanism:** 2024/1709 calls the reference's even `‖(g,−f)‖²` "an idiosyncrasy of the C implementation itself, that is easily fixable." It is load-bearing: `zint_bezout` requires both inputs odd (its return is literally `(gcd==1) & x[0] & y[0]`), `keygen.c:4190` states "we require that Res(f,φ) and Res(g,φ) are both odd (the NTRU equation solver requires it)", `Res(f,φ) ≡ f(1)^n`, and `t ≡ f(1)+g(1) (mod 2)` — so the **constant-time binary GCD forces the even norm**. That directly substantiates this project's headline §7.1 break against its most obvious referee objection. **(b) Constructive:** this implementation may be the only one that can build §7.1-conformant keys today, because `src/ntrugen.jl:403` uses a generic `xgcd_floor` that tolerates an even operand. **(c) The attack:** with odd-norm keys in hand, test the half of §7.1 the paper concedes is sound. `m₂ = t²−2u²` and `m₃ = t³−2t(u²+v²+w²)+2u(v−w)²` read as toy-degree derivations; at n=512 the last six leaves sit nine levels down a tree whose `l10` correction mixes in other nodes' rationals. If any of the twelve in-precision denominators acquires a factor of 2 that `t`'s parity does not control, **§7.1 is broken outright, not merely undeployable.** **(d)** And a small non-auditability theorem: whether a round-based signer protects a key depends on the parity of `‖(g,−f)‖²`, which no party without the private key can check (`q` odd ⇒ `Λ + 2Z^{2n} = Z^{2n}`; `φ(G,−F) ≡ 1` forced by `fG − gF = q`) — so a subverted generator silently keeps victims exploitable.

**First experiment.** Patch `ntru_gen` with `parity = :one_odd` so exactly one of `f,g` has odd coefficient sum. Generate ~500 keys and check per key: `t` odd; `ntru_equation_holds`; `gs_norm_ok`; `max|F|,|G|` within `max_fg_bits`; encode/decode round-trip; sign+verify. Record `NTRUSolveFailure` rate and wall clock vs standard parity. Then recompute the `2n` centres for one signature at `BigFloat(200)` and recover each of the twelve in-precision denominators by **best-rational-approximation** (theory-independent — do not trust the published `m₂`, `m₃` formulas), and tabulate 2-adic valuations for ~50 odd-norm and ~50 even-norm keys. Finally re-run `scripts/countermeasure_break.jl` on odd-norm keys and scan ≥10⁵ signatures for a half-integer centre at calls 2n−5..2n.

**Kill criterion.** Dead as a break if all twelve recovered denominators are odd for every odd-norm key **and** the NewSamplerZ campaign produces zero half-integer centres at the last six calls (against ~6 expected under even norm). Dead as a *deployment* claim if a ≤20-line constant-time `zint_bezout` patch (`gcd(x,y) = gcd(x, x+y)` with `x` odd, back-transforming `x(u'−v') − y v' = 1`) handles one even input with no measurable cost — in which case the honest move is to **retract** the paper's "undeployable" framing.

**Compute.** 500 keygens: seconds to minutes. BigFloat(200) centres: ~0.1 s each, so 100 keys × a few signatures is trivial. NewSamplerZ campaign 10⁵ signatures ≈ 4–8 min. **~2 hours**, plus generic-type work on the tree code.

**What is most likely to go wrong.** Either outcome of (a)/(b) is a footnote-sized result on its own — this card only pays if (c) lands, and (c) is a long shot because the paper has already reasoned (correctly, on its face) that odd `t` ⇒ odd `m₂`, `m₃`. Also: the honest outcome may *weaken* an established claim of ours. That is fine — reviewers weight self-corrections heavily — but budget for it.

---

### 12. ⌊−0⌋: the two supported backends of the normative implementation disagree about what `floor` means

**Card 91.**

**Hypothesis.** FALCON's determinism rests on one operation, `s = fpr_floor(μ)` (`sign.c:1369`), and the reference ships **two supported floating-point backends with different semantics on the −0 bit pattern**: emulated `fpr.h` returns −1 (its own comment calls the choice "debatable"); native `fpr.h` (`r = (int64_t)x.v; r − (x.v < (double)r)`) returns 0. The comment excuses this as unreachable because "the other functions normalize zero to +0" — and that excuse is questionable: `fpr_neg` is `x ^= 1<<63`, `fpr_add` documents `(−0)+(−0) → −0`, `FPR()` preserves the sign when clamping an underflow, and `sign.c:815/992` apply `fpr_neg(ni)` to the entire `t1` branch. If a sampler centre is ever exactly zero with the sign bit set, two builds of the **normative** implementation emit different signatures — deterministically, with no perturbation involved. That is a **fourth, semantic class**: not how an operation is spelled, but a disagreement about what `⌊·⌋` *means*, and the specification just writes `⌊μ⌋`.

**First experiment.** Step 1 (minutes, decisive on the premise): compile a ~25-line harness against `scripts/cref/fpr.h` twice — once emulated, once `FALCON_FPNATIVE` — and print `fpr_floor(fpr_neg(fpr_zero))` and `fpr_floor(0x8000000000000000)`. Expect −1 vs 0. Step 2 (seconds): instrument `src/samplerz.jl` to count, per leaf call, `μ == 0.0` and `signbit(μ)`, bucketed by position, over ~2000 signatures (≈2×10⁶ leaf calls). Step 3: audit reachability of −0 by replaying one signing run through a Julia model of `fpr_add`/`fpr_mul`/`fpr_neg` sign rules, logging every zero intermediate and its sign bit. Step 4: if a negatively-signed zero reaches a centre, replay under both floor semantics and diff the bytes; then confirm §7.2 snapping removes it (predicted: `round(g·(−0.0)) = 0`, so yes).

**Kill criterion.** Dead if over ~10⁷ leaf calls no centre is ever exactly ±0.0 **and** the emulated-backend replay finds no route by which a negatively-signed zero reaches `t0`/`t1` at a leaf. Reduced to a documentation defect if both backends in fact return 0.

**Compute.** **Under an hour, total.** Highest feasibility on the board; self-kills or self-confirms before lunch.

**What is most likely to go wrong.** Reachability. §6.2 measured computed centres sitting 1–3 ulp off the exact rational, i.e. they do not land exactly on integers, let alone on −0. The realistic landing is "a documentation/robustness defect in the reference plus a normative gap in the spec's `⌊·⌋`," which is worth a paragraph in the NIST comment, not a paper. Run it anyway because it costs an hour.

---

## RUNNERS-UP (next 15, one line each)

13. **ACVP vs FIPS 206 contradiction (75+74)** — ACVP's randomized `sigGen` supplies `rnd` to the module, so every validated FN-DSA module must expose the derandomization hook the standard intends to forbid; `usnistgov/ACVP` clones here, and PCT power for a randomized signer is *exactly zero* because divergent signatures verify.
14. **Arithmetic weak keys (53+86)** — factor the integer-centre condition by CRT over the prime powers of `t`; the established weak-key negative excluded only the *size* of `‖(g,−f)‖²`, never its factorisation, and §7.1's 11-vs-6.1 half-integer excess (p≈0.04) is a live hint.
15. **The AM/HM tree theorem (85)** — the `2n` leaves are the arithmetic/harmonic-mean binary iterates of the Gram embeddings, all-AM leaf `= ‖(g,−f)‖²` exactly and all-HM `= q²/that`, which would explain `σ_min = σ/(1.17√q)` to eleven digits and give the whole exposure profile from the key with no signing.
16. **`round`'s tie rule is itself underdetermined (66)** — `floor` is tie-free and mode-independent; `round` is ties-to-even in Julia/Python, ties-away in C, MXCSR-controlled via `rint` — a second, arithmetic break axis of §7.1 that compounds with the parity one.
17. **The last-six key-recovery window (5+37+57)** — the merge tower gives *every* leaf 2-sparse support `x^s(a + b x^{n/2})`, so §5.1 generalises with a known shift; force the divergence rather than waiting for a 10⁻⁹ event.
18. **The precision paradox (33)** — two implementations at *different* precision share no rounding offset, so `P(straddle|integer) → 1/2` and a 113-bit tree should disagree with double ~20× more often than two double spellings do, killing "compute it in higher precision to be safe."
19. **Reduction quality drives the leak rate (89)** — adding `k·(f,g)` to `(F,G)` leaves `h`, the GS norms and the sampled point invariant but scales `|L10|`, hence `|Δμ|`; invisible in the public key, chosen at keygen, and a normative requirement on Reduce that nobody has stated.
20. **The FFT radix is underdetermined too (39)** — radix-4/split-radix/Stockham perturbs the *target* `t`, not only the tree, so unlike A2 it should reach calls 1,2; ten-minute no-op check before any rate run, and it is a prime candidate for card 1's oracle window.
21. **The precision floor `p*` (36)** — locate the mantissa width (~32 bits) at which the four structured positions give way to `2n` generic ones; a number FIPS 206 needs for the FPGA/Cortex-M regime and does not have.
22. **The straddle law, measured not assumed (54)** — replace the hand-set "2 ulp spread" with the directly computable shared offset `δ = (gμ̂ − round(gμ̂))/g`, validate against four measured rates, and make every future arm cost minutes instead of days.
23. **Lemma 2 as a fault-targeting rule (65)** — width faults do reach the signature at large ε but land at *uniformly random* positions where §5.1 recovery fails, while far smaller centre faults recover keys; a hardening-budget rule in both directions.
24. **The fault-susceptibility map of the expanded key (63)** — find single bit flips in the long-lived tree whose `|Δμ|` lands in the silent band `[3×10⁻¹⁴, 1/2n]`: reaching, confined to the last two calls, and invisible to the operator.
25. **Calibrated ulp-ladder KAT vectors (42)** — a ~40-vector suite labelled by the ulp threshold each detects, turning pass/fail conformance into a measurement; honest per-vector power is ~0.1–0.3, not 1.
26. **Spec-conformant complex division is less constant-time (67)** — Julia's/Python's robust division branches on secret Gram entries and does two divides where `FPC_DIV` is branch-free with one, so the spec's silence is normative for side-channel behaviour, not only determinism.
27. **What a first-end divergence actually leaks (25+60)** — settle 2024/1709 §5's unexperimented dismissal with a byte-consumption count and a `D_{Λ,√2σ}` moment test; the *negative* is the useful outcome, because it makes the event channel (card 1) the only channel.

---

## KILLED — the negative space, with reasons

**Already established (in this project or 2024/1709) — do not re-run:**
- **10, 35, 45, 51** — all four are verbatim §7.2 rational snapping (`n = round(g·μ̂)`, `s = fld(n,g)`, first six and last six calls), already implemented in `scripts/snapping.jl`, margins already measured (1.4e-9 / 2.4e-8 vs 0.5), all four known divergent pairs already replayed to AGREE, 0 of 1200 ordinary signatures moved.
- **58, and the identity half of 21** — §6.7/`scripts/linear_form.jl` already proves `t1 = (c·f)/q` exactly against Int128 negacyclic arithmetic and already states Heuristic 1 becomes a theorem at calls 1,2 with density exactly `1/q`.
- **16** — store-vs-rebuild *is* `sign_tree` vs `sign_dyn` (2024/1709 §6.1, reproduced as the A1 control), and the proposed measurement is verbatim the established A2 arm.
- **11, 13, 17, 19, 99(partly)** — hedged signing, IBE extraction, VM snapshot replay and load-balanced fleets are all "derandomized FALCON", the threat model 2024/1709 explicitly scopes to; no new mechanism, position or rate. (Keep only the *artifact*: `/home/user/algorand/falcon/deterministic.c` is a shipped derandomized Falcon-1024 — cite it, don't build a card on it.)
- **14, 34, 61** — the rounding-mode result is §7.3 (1800/1800 vs 0/1800 control), and the amplification claim is precisely what §7.3 refuted (0 recoveries in 48 000 tail-window pairs). Residual worth one paragraph: FTZ/DAZ via `-Ofast`/`crtfastmath.o`, the subnormal census, and FPEMU-vs-native diverging by construction.
- **4, 18, 29, 40, 71, 72, 97** — the KAT-power calculation (`≈3/ρ` vectors) is one line of arithmetic on an already-measured rate, restated seven times, and is **subsumed** by §7.3's far stronger rate-1 result: bit-exact KATs are unachievable at *any* suite size until the FP environment is normative.
- **6** — "round protects iff the denominator is odd" is verbatim §7.1(iii); and the interior atlas has no consumer (§6.4: interior denominators exceed 2⁵³, zero integer centres in 1.02×10⁸ draws).
- **48** — same as 6; the interior-denominator negative is established at a level of rigour this would only restate.
- **30** — wholly conditional on the established §6.7, and its parity half is already proved in §7.1(iii).
- **96, 15** — "the rate is set by the integer-centre density, not by the perturbation" is §6.6; the maintenance catalogue's entries are the known A1/A2 set. (The *completeness count* survives inside card 10.)
- **52, 83, 90, 95** — keygen Babai rounding: `#064` already ran the differential (identical seeds, both `div_fft` spellings, 0/400, live positive control), and the Babai ratio's denominator is a resultant, so the margin is astronomically safe. (The *rejection-threshold* sites survive as card 5.)
- **12** — 2024/1709 §5's infeasibility claim is better overturned by §6.7's exact `F_q` solve than by BKZ on desynchronised difference vectors.
- **20** — SNARK completeness reduces to the already-measured divergence rate with a relabelled second arm; keep the framing sentence, drop the card.
- **70** — `s = (t − z)B`, so a post-descent change to one `z` coefficient changes `s` by one basis row times a monomial: true by construction, and it restates the structure the working recovery already depends on.

**Premise factually refuted:**
- **1, 2, 8** — all three assume 2024/1709 §7.1 proposes `floor → round` as a drop-in on Algorithm 15. It does not: Algorithm 4 replaces the base sampler too, with `x = (y−r)²/2σ² − (y₊²−y₊)/2σ_max²`, under which `x ≥ 0` for every `r ∈ [−1/2,1/2)`. §7.1(i) already implemented exactly this and verified it distributionally sound (χ²/df 0.56–1.02 at five `(μ,σ)`, 200 000 samples each). No acceptance clamp ⇒ no correctness break (1), no per-leaf bias to regress on GS norms (2), no isochrony kink (8).

**Infeasible on this machine:**
- **26** — belief propagation over `F_q` with q=12289 is ~10¹¹ message-vector ops per sweep; "exhaustive posterior at toy n=16" is `12289¹⁶` states; and `src/params.jl` has no toy degrees. *Salvage the five-line information-theoretic floor (≈3×10⁶ observations) and nothing else.*
- **6 (as written)** — a `Rational{BigInt}` lockstep run of ffSampling is impossible: the FFT domain is over complex roots of unity, not `Q`. Use continued-fraction reconstruction (card 27-runner-up / card 11's BigFloat(200) route) instead.
- **80** — a fixed-point (Int128, 2⁻⁸⁰) FFT and ffSampling tower from scratch plus a 200-bit oracle over 10⁵ signatures is days, not hours. *Salvage the exactness sub-theorem: how many fractional bits make `floor` provably exact at the six in-precision positions.*
- **12, 27 (n=512 endpoint)** — dimension-1024 BKZ is out of reach; `fpylll` is installed and n≤128 sweeps run, but the headline must be reported as extrapolation.
- **100 (OpenFHE half)** — extracting `DiscreteGaussianGeneratorGeneric` into a standalone compiling harness is the classic C++ header yak-shave; scope to Lattigo plus one extracted file, and expect the null (FHE sampler decisions are integer- or CDT-indexed, so condition (b) of the two-condition characterisation fails).

**Real but low impact / dominated by a better card:**
- **9** — DAG-diff tooling whose validation target is the divergence already found by hand; dominated by card 7's certificate.
- **22, 27, 54** — supporting measurements, not results; fold into cards 1 and the runners-up.
- **38, 44, 68** — three more instances of the sensitivity-instrument cluster; card 7 (certificate + inventory) has the cleanest soundness story and the completeness payoff.
- **46** — the precision-independence lower bound is already §6.6; the sound upper bound will very likely go vacuous through ten levels with `div_fft` in the loop (test at n=16 in minutes before investing).
- **47** — dominated by card 9's framing of the same HAWK question.
- **49** — dominated by card 10; e-graph saturation on the FFT butterfly is the likely blow-up.
- **50** — survives inside card 6 as the static survey; standalone it is survey work whose likely outcome is "everything is a transliteration."
- **55, 32** — absorbed into card 7; the BerExp Lipschitz lemma is a lemma, not a paper.
- **69** — the FPU-latency channel is very likely a clean negative on x86-64 (no subnormals at these magnitudes, fixed-latency `divsd`), and this is a Firecracker guest so sub-cycle timing claims are untrustworthy. Report the negative in one sentence.
- **73** — the trilemma probably collapses to "publish the sampled `z`", which is the card's own kill test and takes an afternoon.
- **74** — one table and one paragraph; folded into runner-up 13.
- **84** — shares its premise with card 5, which has a 1000× larger and already-documented disagreement band.
- **87** — the acceptance window caps the achievable rate boost at ~1.5–2×, and shorter secrets pay for it in lattice security.
- **88** — a small advisory theorem; folded into card 11.
- **93/94** — merged into card 7.
- **98** — merged into card 9.

---

## If you can only run ONE experiment next

**Run card 1(a)+(b): the oracle window.** It costs almost nothing — part (a) is *re-reading a histogram `scripts/rounding_mode.jl` already computes* — and it is the only experiment that simultaneously (i) tests the project's most valuable unfinished claim, since §6.7 explicitly names the missing oracle and this supplies the exact two-sided condition `|η| ≪ |Δμ| ≪ (2/q)/2n` under which one exists, (ii) *reinterprets* the §7.3 negative result rather than repeating it — directed rounding fails not because the leak is bounded but because it lands **above** the window, which turns "rate 1 and 0 recoveries" from an embarrassment into a design law, (iii) tells you immediately whether cards 3, 18, 19, 20 and 24 are worth running, because every one of them is a search for a realistic perturbation inside that window, and (iv) has a symmetric, publishable outcome either way: if the window is non-empty you get an end-to-end key recovery by Gaussian elimination that contradicts 2024/1709 §5's infeasibility claim; if it is empty you get a theorem bounding the entire first-two channel, which is the cleanest possible closure of a direction the paper currently leaves open. The one discipline to impose before writing anything: state plainly, in the first paragraph, that §6.7 already publishes the linear form and the attack shape, and that the new content is the window, the reinterpretation of directed rounding, and the executed solve — the project has been burned three times by claiming something already in 2024/1709, and the fourth time would be claiming something already in its own manuscript.