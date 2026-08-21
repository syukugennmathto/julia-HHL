#!/usr/bin/env julia
#
# dudect.jl -- timing-leakage detection by Welch's t-test.
#
#     julia --project=falcon falcon/scripts/dudect.jl [measurements]
#
# Implements the "dude, is my code constant time?" procedure (Reparaz, Balasch,
# Verbauwhede, DATE 2017) against this implementation.  The method needs no
# model of the machine: it times the routine under two input classes and asks
# whether the two timing distributions are distinguishable.  If they are, the
# routine leaks; if they are not, nothing has been proved, but the experiment
# has bounded what a timing attacker at this sample size could have seen.
#
# ---------------------------------------------------------------------------
# WHY THE CONTROLS ARE NOT OPTIONAL
# ---------------------------------------------------------------------------
#
# A leakage test that reports "no leak" is worthless unless it is known to
# report "leak" when there is one.  Every run therefore includes:
#
#   * a POSITIVE control -- a function whose running time is a blatant
#     function of its secret input.  If this does not trip, the harness is
#     broken and every other line of the report is meaningless.
#   * a NEGATIVE control -- a function whose running time cannot depend on its
#     input.  If this trips, the harness is producing false positives, and the
#     real results are noise.
#
# This is the same discipline as docs/debug_log.md #035 (the profiler that
# hid the cost it was meant to measure) and #050 (the baseline that already
# agreed): an instrument is not to be trusted until it has been shown to move
# in both directions.
#
# ---------------------------------------------------------------------------
# METHOD, AND THE CHOICES INSIDE IT
# ---------------------------------------------------------------------------
#
#   - The two classes are INTERLEAVED, one measurement each, alternating.  A
#     block of class A followed by a block of class B would let CPU frequency
#     drift, cache state or a garbage collection wander in and be attributed
#     to the class.  Interleaving makes any such drift common-mode.
#
#   - The class is chosen by a coin flip per iteration, not by strict
#     alternation, so that a periodic disturbance at the alternation frequency
#     cannot alias onto the class label.
#
#   - THE EFFECT SIZE IS REPORTED AS A DIFFERENCE OF MEDIANS, not of means.
#     `compress_sig` allocates its output buffer, so its timing distribution
#     has a collector-driven tail that no amount of cropping fully removes;
#     reported as a mean difference it came out at 6.35 microseconds, larger
#     than the routine's entire median runtime, which is incoherent.  The
#     median spread between the same inputs measured directly was 90 ns.  The
#     t statistic stays as the *detection* statistic -- it is what answers "are
#     these distinguishable" -- and the median difference answers "by how
#     much".  Reporting only t invites reading it as a magnitude, which it is
#     not: a routine made faster and less variable will show a LARGER t for
#     the SAME absolute leak (docs/debug_log.md #051).
#
#   - Measurements are cropped at a range of upper percentiles before the
#     t-test, as in the original dudect.  The raw distribution has a heavy
#     right tail from the collector and the scheduler, which is not the signal
#     and which inflates the variance enough to hide one.  The reported
#     statistic is the largest |t| over the crops, which is the standard (and
#     deliberately conservative-in-the-wrong-direction) choice: it multiplies
#     the number of hypotheses, so a threshold of 4.5 rather than 1.96 is used.
#
#   - GC is NOT disabled.  Allocation is a real cost of this implementation
#     (docs/debug_log.md #035), and a collection whose timing depends on how
#     much a secret-dependent path allocated is a real leak, not an artefact.
#     What the cropping removes is the tail, not the mean shift.
#
#   - A warm-up pass runs before timing.  In Julia the first call is
#     compilation; mixing it in would swamp everything.
#
# THRESHOLD.  |t| > 4.5 is dudect's convention and corresponds to a false
# positive probability around 1e-5 for a single test.  It is a detection
# threshold, not a security bound: |t| < 4.5 at N measurements means "no
# leakage large enough to see with N measurements", and the corresponding
# statement in the report is deliberately worded that way.

using Falcon
using Printf
using Random
using Statistics

const THRESHOLD = 4.5

"""
    NOISE_FLOOR

The largest |t| produced in this run by a routine that *cannot* leak, filled in
by the controls before any real result is reported.

dudect's fixed |t| > 4.5 is a per-test false-positive threshold under the
assumption that the timing noise is independent across measurements.  On this
setup that assumption fails as `n` grows: at n = 100000 the input-free
`basesampler` self-check reached **t = 11.61**, with no input to leak.  Slow
drift -- CPU frequency, page cache, the allocator's state -- is common to both
classes in the mean but not independent measurement to measurement, and the
t statistic grows as sqrt(n) against it.

So a run reports its own floor and judges against it.  A result counts as
leakage only if it clears both 4.5 and three times whatever the controls
managed on this machine, on this day, at this `n`.  Publishing the floor
alongside is the point: it is what makes "not detected" a quantitative claim
instead of an absence of evidence.
"""
const NOISE_FLOOR = Ref(0.0)

# ---------------------------------------------------------------------------
# Welch's t-test
# ---------------------------------------------------------------------------

"""
    welch_t(a, b) -> Float64

Welch's t statistic for two samples of unequal variance.  Returns 0.0 when
either sample is too small or degenerate, so that a crop that removed
everything cannot masquerade as a clean result.
"""
function welch_t(a::AbstractVector{Float64}, b::AbstractVector{Float64})
    (length(a) < 2 || length(b) < 2) && return 0.0
    ma, mb = mean(a), mean(b)
    va, vb = var(a), var(b)
    denom = sqrt(va / length(a) + vb / length(b))
    denom == 0 && return 0.0
    return (ma - mb) / denom
end

"""
    max_abs_t(a, b) -> (t, crop)

The largest |t| over a family of upper-percentile crops of the *pooled*
distribution, and which crop produced it.

Cropping at a percentile of the pooled sample rather than of each class
separately matters: cropping each class at its own percentile would remove a
genuine mean difference by construction.
"""
function max_abs_t(a::Vector{Float64}, b::Vector{Float64})
    best = welch_t(a, b)
    bestcrop = 1.0
    # The effect size is measured on the UNCROPPED samples, and deliberately
    # not at whichever crop maximised |t|.  Taking it from the winning crop
    # gave -6465 ns for `compress_sig`, larger than the routine's entire
    # 5.2 microsecond median runtime, because at a crop that separates the
    # classes one of them keeps only the few samples nearest the cut.  Timed
    # directly, the same inputs differ by about 90 ns.  The crop is a device
    # for making the *test* sensitive; it must not be used to state the
    # *magnitude* (docs/debug_log.md #051).
    dmed = median(a) - median(b)
    pooled = sort(vcat(a, b))
    for p in (0.999, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
        cut = pooled[max(1, min(length(pooled), round(Int, p * length(pooled))))]
        ac = filter(<=(cut), a)
        bc = filter(<=(cut), b)
        t = welch_t(ac, bc)
        if abs(t) > abs(best)
            best = t
            bestcrop = p
        end
    end
    return (best, bestcrop, dmed)
end

# ---------------------------------------------------------------------------
# The measurement loop
# ---------------------------------------------------------------------------

"""
    measure(prepare_fixed, prepare_random, run; n, warmup) -> (t, crop, na, nb)

Time `run(x)` under two input classes and return the t statistic.

## Both classes' inputs are materialised BEFORE the loop

The first version called `prepare_fixed()` / `prepare_random()` inside the
loop, just outside the timed region.  The negative control then reported
`t = -5.46` -- a false positive, and the harness said so.  The cause was in
the harness, not the code under test: `prepare_random()` does work
(`rand`, an allocation) that `prepare_fixed()` does not, so the two classes
entered the timed region with the cache and the branch predictor in
systematically different states.  The asymmetry was outside the measurement
and still inside the experiment.

Pre-building both input vectors and indexing them removes it: every iteration
now does one array read before `time_ns()`, whichever class it is.

That is the same failure as docs/debug_log.md #035, one level up -- there the
instrument moved a cost out of the timed region, here it moved one in.
"""
function measure(prepare_fixed, prepare_random, run; n::Int = 20_000,
                 warmup::Int = 200, rng = MersenneTwister(0xC0FFEE))
    half = n ÷ 2 + 8
    inputs_a = [prepare_fixed() for _ in 1:half]
    inputs_b = [prepare_random() for _ in 1:half]

    for i in 1:warmup
        run(inputs_a[mod1(i, half)])
        run(inputs_b[mod1(i, half)])
    end

    ta = Float64[]
    tb = Float64[]
    sizehint!(ta, half)
    sizehint!(tb, half)
    ia = 1
    ib = 1

    for _ in 1:n
        fixedclass = rand(rng, Bool)
        x = fixedclass ? inputs_a[ia] : inputs_b[ib]
        t0 = time_ns()
        run(x)
        dt = Float64(time_ns() - t0)
        if fixedclass
            push!(ta, dt); ia = mod1(ia + 1, half)
        else
            push!(tb, dt); ib = mod1(ib + 1, half)
        end
    end

    t, crop, dmu = max_abs_t(ta, tb)
    return (t, crop, length(ta), length(tb), dmu)
end

"Does `t` clear both the fixed threshold and the run's measured noise floor?"
significant(t) = abs(t) > THRESHOLD && abs(t) > 3 * NOISE_FLOOR[]

function report(name, t, crop, na, nb, dmu = NaN; expect::Symbol = :unknown)
    verdict = significant(t) ? "LEAK" : "  --"
    ratio = NOISE_FLOOR[] > 0 ? abs(t) / NOISE_FLOOR[] : NaN
    mark = if expect === :leak
        significant(t) ? "ok" : "HARNESS BROKEN"
    elseif expect === :clean
        significant(t) ? "FALSE POSITIVE" : "ok"
    elseif expect === :floor
        "-> noise floor"
    else
        ""
    end
    @printf("%-44s  t = %8.2f  %5.1fx  dmed = %+8.1f ns  n = %d/%d  %s  %s\n",
            name, t, ratio, dmu, na, nb, verdict, mark)
    return significant(t)
end

"Record a control that cannot leak, raising the run's noise floor."
function note_floor!(name, t, crop, na, nb, dmu = NaN)
    NOISE_FLOOR[] = max(NOISE_FLOOR[], abs(t))
    @printf("%-44s  t = %8.2f   ----  dmed = %+8.1f ns  n = %d/%d  -> noise floor\n",
            name, t, dmu, na, nb)
    return false
end

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

"A long byte string the sampler can draw from without the source itself varying."
const POOL = shake256(codeunits("dudect/pool"), 1 << 22)

"""
    distinct(make, k) -> Vector

`k` separately allocated objects, each built by `make()`.

The fixed class must not be a single object that every iteration re-reads.
The first version made class A one array and class B one of several, and
`sqnorm` came out at `t = -16.3` -- which is not a property of `sqnorm` at
all: class A read the same cache-hot buffer every time while class B walked a
working set.  The values were what differed by design; the *memory* differed
by accident, and the accident is larger.

So both classes now cycle through the same number of distinct allocations.
Only the contents differ.
"""
distinct(make, k::Int) = [make() for _ in 1:k]

"Return the next element of `pool`, advancing the cursor.  A function rather
than an inline `pool[(i[] = mod1(i[]+1, k))]` because that is a keyword
argument as far as Julia's parser is concerned."
function nextof!(pool::Vector, cursor::Base.RefValue{Int})
    cursor[] = mod1(cursor[] + 1, length(pool))
    return pool[cursor[]]
end

"""
A sink for results that must not be optimised away.

The first version of the positive control was a `for i in 1:x; s += i; end`
loop whose result was discarded.  It did not trip, and the harness correctly
reported itself broken: LLVM turns that loop into `x*(x+1)/2` and then deletes
it, so the "secret-dependent" running time was constant.  Touching memory the
compiler cannot reason about, and storing the result somewhere observable,
removes both escapes.
"""
const SINK = Ref(UInt64(0))

mutable struct PoolSource
    pos::Int
end
function (p::PoolSource)(k::Integer)
    p.pos + k > length(POOL) && (p.pos = 0)
    out = POOL[(p.pos + 1):(p.pos + k)]
    p.pos += k
    return out
end

function main()
    n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20_000
    rng = MersenneTwister(20260821)
    p = FALCON_512
    println("# dudect.jl  julia=", VERSION, "  measurements=", n,
            "  threshold |t| > ", THRESHOLD)
    println("#")
    println("# class A = fixed secret input, class B = random secret input.")
    println("# 'LEAK' means the two timing distributions are distinguishable.")
    println()

    broken = false

    # --- controls ----------------------------------------------------------
    # The floor is established FIRST, from routines that cannot leak, and
    # every later verdict is relative to it.
    println("## controls (the two negative ones set this run's noise floor)")
    let
        src = PoolSource(0)
        t, c, na, nb, dm = measure(() -> nothing, () -> nothing,
                              _ -> basesampler(src); n = n)
        note_floor!("negative control (basesampler, no input at all)", t, c, na, nb, dm)
    end
    let
        # POSITIVE: running time is a plain function of the secret, and the
        # loop reads memory so it cannot be closed-formed away.
        t, c, na, nb, dm = measure(() -> 64,
                              () -> rand(rng, 1:4096),
                              x -> begin
                                  s = UInt64(0)
                                  @inbounds for i in 1:x
                                      s = xor(s * 0x9e3779b97f4a7c15, UInt64(POOL[i]))
                                  end
                                  SINK[] = s
                              end; n = n)
        broken |= !report("positive control (loop count = secret)", t, c, na, nb, dm;
                          expect = :leak)
    end
    let
        # NEGATIVE: same work regardless of input.
        t, c, na, nb, dm = measure(() -> 0x5a5a5a5a5a5a5a5a,
                              () -> rand(rng, UInt64),
                              x -> begin
                                  s = x
                                  @inbounds for i in 1:4096
                                      s = xor(s * 0x9e3779b97f4a7c15, UInt64(POOL[i]))
                                  end
                                  SINK[] = s
                              end; n = n)
        note_floor!("negative control (fixed work)", t, c, na, nb, dm)
    end
    println()

    # --- the implementation ------------------------------------------------
    println("## verification -- a second control, out of this codebase")
    let
        # Verification handles no secret, so it does not owe constant time.
        # It earns its place here as a control made of real code: an invalid
        # signature is rejected on a different path from a valid one, so the
        # two classes MUST be distinguishable.  A harness that cannot see this
        # cannot be trusted to see anything.
        r = chacha20(collect(UInt8, 0x00:0x37))
        sk, pk = falcon_keygen(512, k -> randombytes!(r, k))
        msg = collect(codeunits("dudect verification message"))
        sig = falcon_sign(sk, msg, k -> randombytes!(r, k))
        bad = copy(sig); bad[end] = xor(bad[end], 0xff)
        broken |= !report("positive control: verify, valid vs tampered",
                          measure(() -> sig, () -> bad,
                                  x -> falcon_verify(pk, msg, x); n = n)...;
                          expect = :leak)

        # ...and with both classes on the accepting path, which is the
        # comparison that would matter if verification did handle a secret.
        msgs = [collect(codeunits("dudect verify $i")) for i in 1:8]
        sigs = [falcon_sign(sk, m, k -> randombytes!(r, k)) for m in msgs]
        idx = Ref(1)
        t, c, na, nb, dm = measure(() -> 1,
                              () -> rand(rng, 2:8),
                              i -> falcon_verify(pk, msgs[i], sigs[i]); n = n)
        report("verify, fixed vs varying valid signature", t, c, na, nb, dm)
    end
    println()

    println("## the sampler (the part that must be isochronous)")
    let
        sigmin = p.sigma_min
        sigma = 1.5
        src = PoolSource(0)
        # sigma fixed, mu fixed vs mu random.  `mu` is the secret: in signing
        # it is a coordinate of the target expressed in the secret basis.
        t, c, na, nb, dm = measure(() -> 0.375,
                              () -> rand(rng) * 100 - 50,
                              mu -> samplerz(mu, sigma, sigmin, src); n = n)
        report("samplerz, fixed vs random centre (sigma fixed)", t, c, na, nb, dm)
    end
    let
        sigmin = p.sigma_min
        src = PoolSource(0)
        mu = 0.375
        # mu fixed, sigma fixed vs sigma random in the admissible range.
        t, c, na, nb, dm = measure(() -> 1.5,
                              () -> sigmin + rand(rng) * (p.sigma_max - sigmin) * 0.98,
                              s -> samplerz(mu, s, sigmin, src); n = n)
        report("samplerz, fixed vs random width (centre fixed)", t, c, na, nb, dm)
    end
    println()

    println("## signing")
    let
        r = chacha20(collect(UInt8, 0x00:0x37))
        sk1, _ = falcon_keygen(512, k -> randombytes!(r, k))
        sk2, _ = falcon_keygen(512, k -> randombytes!(r, k))
        keys = [sk1, sk2]
        for i in 3:6
            push!(keys, falcon_keygen(512, k -> randombytes!(r, k))[1])
        end
        msg = collect(codeunits("dudect signing message"))
        src = chacha20(collect(UInt8, 0x40:0x77))
        sb = k -> randombytes!(src, k)
        # fixed key vs a key drawn from a small pool -- what a timing attacker
        # watching a signer would be trying to separate.
        nsig = min(n, 20_000)
        NK = 6
        pool_a = distinct(() -> expand_privkey(keys[1].f, keys[1].g, keys[1].F, keys[1].G, p), NK)
        pool_b = distinct(() -> (k = keys[rand(rng, 2:length(keys))];
                                 expand_privkey(k.f, k.g, k.F, k.G, p)), NK)
        ia = Ref(0); ib = Ref(0)
        report("falcon_sign (whole), fixed vs varying key",
               measure(() -> nextof!(pool_a, ia),
                       () -> nextof!(pool_b, ib),
                       k -> falcon_sign(k, msg, sb); n = nsig)...)

        # --- bisected -------------------------------------------------------
        # Where in signing does the difference live?  Each stage below is
        # driven with the same two key classes, so the numbers are comparable
        # to the whole-signature line above.
        pt = hash_to_point(msg, shake256(codeunits("dudect salt"), SALT_LEN),
                           512; q = p.q)
        report("  sample_preimage (ffSampling + basis)",
               measure(() -> nextof!(pool_a, ia),
                       () -> nextof!(pool_b, ib),
                       k -> Falcon.sample_preimage(k, pt, sb); n = nsig)...)

        # Both classes cycle through NPOOL distinct allocations; only the
        # contents differ.  See `distinct`.
        NPOOL = 8
        pre = [Falcon.sample_preimage(k, pt, sb) for k in keys]
        s2_a = distinct(() -> Int.(pre[1][2]), NPOOL)
        s2_b = distinct(() -> Int.(pre[rand(rng, 2:length(pre))][2]), NPOOL)
        s1_a = distinct(() -> Int.(pre[1][1]), NPOOL)
        s1_b = distinct(() -> Int.(pre[rand(rng, 2:length(pre))][1]), NPOOL)
        ka = Ref(0); kb = Ref(0)
        report("  sqnorm (generic, scans coefficients)",
               measure(() -> (ka[] = mod1(ka[] + 1, NPOOL); (s1_a[ka[]], s2_a[ka[]])),
                       () -> (kb[] = mod1(kb[] + 1, NPOOL); (s1_b[kb[]], s2_b[kb[]])),
                       x -> sqnorm(x[1], x[2]); n = n)...)
        ka[] = 0; kb[] = 0
        report("  sqnorm_machine (what signing uses now)",
               measure(() -> (ka[] = mod1(ka[] + 1, NPOOL); (s1_a[ka[]], s2_a[ka[]])),
                       () -> (kb[] = mod1(kb[] + 1, NPOOL); (s1_b[kb[]], s2_b[kb[]])),
                       x -> sqnorm_machine(x[1], x[2]); n = n)...)
        salt = shake256(codeunits("s"), SALT_LEN)
        report("  encode_signature",
               measure(() -> nextof!(s2_a, ka),
                       () -> nextof!(s2_b, kb),
                       x -> encode_signature(salt, x, p.logn, p.sig_bytes); n = n)...)
        report("  compress_sig only",
               measure(() -> nextof!(s2_a, ka),
                       () -> nextof!(s2_b, kb),
                       x -> compress_sig(x, p.sig_bytes - 41); n = n)...)
        keypool_a = distinct(() -> expand_privkey(keys[1].f, keys[1].g, keys[1].F, keys[1].G, p), NPOOL)
        keypool_b = distinct(() -> (k = keys[rand(rng, 2:length(keys))];
                                    expand_privkey(k.f, k.g, k.F, k.G, p)), NPOOL)
        report("  expand_privkey (key setup, not per-signature)",
               measure(() -> keys[1], () -> keys[rand(rng, 2:length(keys))],
                       k -> expand_privkey(k.f, k.g, k.F, k.G, p); n = min(n, 5_000))...)
    end
    println()

    if broken
        println("!! a control failed: the numbers above are not evidence.")
        exit(1)
    end
end

main()
