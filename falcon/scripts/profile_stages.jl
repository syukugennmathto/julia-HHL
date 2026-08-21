#!/usr/bin/env julia
#
# Where does the time actually go?
#
#     julia --project=falcon falcon/scripts/profile_stages.jl [logn]
#
# `bench.jl` answers "how long does keygen/sign/verify take".  This answers
# "and which line inside it is responsible", by timing the stages separately on
# the same inputs and reporting allocation alongside time.  Allocation is
# reported because in this implementation it is usually the answer: nothing
# here is written to reuse a buffer, and Julia's `BigInt` is a heap object with
# a finalizer, so a loop over BigInt coefficients allocates once per operation
# and then pays the collector.
#
# Read the output as a hierarchy: an indented line is a component of the line
# above it.  The components will not sum exactly to their parent -- some of the
# parent is glue, and the timer itself costs something -- but if they sum to
# far less than the parent, that gap is itself a finding.
#
# METHOD NOTES (the ways this measurement could lie, and what is done about
# them):
#
#   - Everything is run once before timing.  In Julia the first call is
#     compilation.  This is the single most common way to publish a wrong
#     benchmark in this language.
#
#   - The median of many runs is reported, not the mean, and not one run.
#     Key generation retries a random number of times, and the GC fires when it
#     fires; both put a long right tail on the distribution.
#
#   - GC time is reported separately (`gc%`).  If a stage's gc% is high, the
#     stage is not slow -- allocation somewhere is, and the collector merely
#     presented the bill to whoever was running at the time.
#
#   - **The collector is NOT forced before each timed call.**  An earlier
#     version called `GC.gc(false)` there, reasoning that each run should start
#     from a comparable heap.  It does -- and it also moves the collection an
#     allocation-heavy stage would have triggered *outside* the timed region,
#     so that stage is billed for its allocation but not for its collection.
#     Measured, that understated signing by 42% (1.745 ms reported as 1.012)
#     and verification by 18%.  Every gc% of 0.0 in the old output was an
#     artefact of the method (docs/debug_log.md #035).  Set
#     FALCON_PROFILE_FORCE_GC=1 to get the old behaviour back for comparison;
#     the honest number is the default.
#
#   - Stages that consume randomness are given a *fresh, identical* byte source
#     each time, so run k and run k+1 do the same work.  Without that, a
#     rejection-sampling stage would drift into different amounts of work as
#     the stream advances and the timings would not be comparable.

using Falcon
using Printf

# `ReplayBytes` is already a `randombytes`-style callable -- it has a
# `(rb::ReplayBytes)(k)` method -- so no wrapper is needed.  Named anyway, so
# that the call sites read the same as the library's own.
bytesource_of(s) = s

# ---------------------------------------------------------------------------
# measurement plumbing
# ---------------------------------------------------------------------------

struct Row
    label::String
    depth::Int
    ms::Float64
    mean_ms::Float64
    alloc_mb::Float64
    gc_pct::Float64
    n::Int
end

const ROWS = Row[]

"""
Whether to run the collector before each timed call.  Default `false`: forcing
it hides the collection that an allocation-heavy stage causes, which is exactly
what this script is trying to find.  See the header.
"""
const _FORCE_GC = Ref(get(ENV, "FALCON_PROFILE_FORCE_GC", "0") == "1")

"""
Time `f` `iters` times and record the median and mean wall time, the allocation
of a single call, and the share of wall time spent in GC.

`setup` runs before each timed call and its result is passed to `f`; use it to
hand each run an identical rewound input.

The GC share is **aggregated over all runs**, not taken from the median run.
Taking it from the median run reported 0.0% everywhere, and that was structural
rather than lucky: collections land in the tail, so the median run is by
definition one that did not collect.  The mean column is there for the same
reason -- for an allocation-heavy stage the mean sits well above the median and
that gap *is* the collector.
"""
function measure!(label, depth, f; iters = 20, setup = () -> nothing)
    f(setup())                                   # compile, and warm any cache
    ts = Float64[]
    total_gc = 0.0
    local alloc
    for _ in 1:iters
        arg = setup()
        _FORCE_GC[] && GC.gc(false)              # off by default; see the header
        st = Base.gc_num()
        t0 = time_ns()
        f(arg)
        dt = (time_ns() - t0) / 1e6
        d = Base.GC_Diff(Base.gc_num(), st)
        push!(ts, dt)
        total_gc += d.total_time / 1e6
        alloc = d.allocd
    end
    sort!(ts)
    wall = sum(ts)
    push!(ROWS, Row(label, depth, ts[length(ts) ÷ 2 + 1], wall / length(ts),
                    alloc / 2^20, wall > 0 ? 100 * total_gc / wall : 0.0, iters))
    return nothing
end

"""
Report a whole operation twice: as it behaves in steady state, and with the
collector run immediately before each call.

The gap between the two is not "GC pause time".  Forcing a collection also
hands the timed call a *compacted* heap, which allocation is much faster into
and which no real workload ever has.  So the pair brackets the truth: the
steady-state figure is what a caller sees, the fresh-heap figure is what the
code would cost if allocation were free, and the difference is the price of
allocating into a live heap -- pauses included.

This is reported explicitly because the `gc%` column reads 0.0 at these
iteration counts and that is honest rather than broken: 20 calls allocating
~1.9 MB each will often not trigger a single collection, so no pause is
attributed even though allocation is dominating.  The finding lives in the
gap, not in the column (docs/debug_log.md #035).
"""
function measure_both!(label, f; iters = 200, setup = () -> nothing)
    was = _FORCE_GC[]
    _FORCE_GC[] = false; measure!(label * " [steady state]", 1, f; iters, setup)
    _FORCE_GC[] = true;  measure!(label * " [after a forced GC]", 1, f; iters, setup)
    _FORCE_GC[] = was
    return nothing
end

function print_rows()
    @printf("%-42s %10s %10s %10s %6s %6s\n",
            "stage", "median ms", "mean ms", "alloc MB", "gc%", "iters")
    println("-"^90)
    for r in ROWS
        @printf("%-42s %10.4f %10.4f %10.3f %6.1f %6d\n",
                "  "^r.depth * r.label, r.ms, r.mean_ms, r.alloc_mb, r.gc_pct, r.n)
    end
    println()
end

"""
A deterministic byte source that is *rewound* before each run, so every run does
equal work.

The bytes are generated once and the same `ReplayBytes` is reused with `pos`
reset.  An earlier version built a new `ReplayBytes` per run, which allocated a
fresh megabyte outside every timed region -- and that turned out to bias the
result the same way an explicit `GC.gc()` does, by making the collector fire
during setup instead of during the call being measured.  Signing came out at
1.34 ms that way against 2.13 ms honestly measured.  Two different mechanisms,
one mistake: *anything* that empties the heap between runs hides the cost of
filling it (docs/debug_log.md #035).
"""
function fresh_source(seed::AbstractString, nbytes::Integer)
    rb = ReplayBytes(shake256(codeunits(seed), nbytes))
    return function ()
        rb.pos = 0
        return rb
    end
end

"""
Wraps a `randombytes` source and counts calls and bytes.  The *number of calls*
is the interesting quantity: each one allocates a `Vector{UInt8}`, so a stage
that draws 20 kB in 5000 calls costs very differently from one that draws it in
three.
"""
mutable struct CountingSource
    inner::Any
    calls::Int
    bytes::Int
end
CountingSource(inner) = CountingSource(inner, 0, 0)
(c::CountingSource)(k::Integer) = (c.calls += 1; c.bytes += k; c.inner(k))

# ---------------------------------------------------------------------------
# the stages
# ---------------------------------------------------------------------------

function main(args)
    logn = length(args) >= 1 ? parse(Int, args[1]) : 9
    n = 1 << logn
    p = params(n)

    println("# profile_stages.jl  julia=", VERSION, "  n=", n)
    println()

    # Set FALCON_PROFILE_SKIP_KEYGEN=1 to measure only signing and verification.
    # One accepted (f, g) is still needed -- there is no signing without a key --
    # but the repeated ntru_solve timings below are what make this script take
    # minutes rather than seconds, and once recorded they rarely need redoing.
    skip_keygen = get(ENV, "FALCON_PROFILE_SKIP_KEYGEN", "0") == "1"

    # ---- key generation -----------------------------------------------------
    # gen_poly draws 4096 samplerz values regardless of n, so its cost is
    # constant in n; ntru_solve's is not.  Timing them apart is the whole point.
    src = fresh_source("falcon-jl/profile/keygen", 1 << 22)

    skip_keygen || measure!("gen_poly (4096 samplerz draws)", 1,
                            s -> gen_poly(n, bytesource_of(s)); iters = 5, setup = src)

    # One accepted (f, g) pair, reused for every stage below.
    f, g = let s = src()
        bs = bytesource_of(s)
        local ff, gg
        while true
            ff = gen_poly(n, bs)
            gg = gen_poly(n, bs)
            gs_norm(Float64.(ff), Float64.(gg); q = p.q) > gram_schmidt_quality()^2 * p.q &&
                continue
            is_invertible_zq(Int.(ff)) || continue
            try
                ntru_solve(ff, gg; q = p.q)
                break
            catch e
                e isa NTRUSolveFailure || rethrow()
            end
        end
        (ff, gg)
    end

    if !skip_keygen
    measure!("gs_norm (the quality rejection)", 1,
             _ -> gs_norm(Float64.(f), Float64.(g); q = p.q); iters = 50)
    measure!("is_invertible_zq", 1, _ -> is_invertible_zq(Int.(f)); iters = 50)
    measure!("ntru_solve (whole descent)", 1,
             _ -> ntru_solve(f, g; q = p.q); iters = 5)

    # ---- inside ntru_solve --------------------------------------------------
    # The descent is  field_norm -> recurse -> lift/karamul -> babai_reduce.
    # Time one level of each at the top degree; the recursion halves the degree
    # each time, so the top level is where the big coefficients live.
    fb = BigInt.(f); gb = BigInt.(g)
    measure!("field_norm (one level, n=$n)", 2, _ -> field_norm(fb); iters = 20)
    measure!("karamul (BigInt, n=$n)", 2, _ -> karamul(fb, gb); iters = 20)
    measure!("galois_conjugate", 2, _ -> galois_conjugate(gb); iters = 50)

    # babai_reduce needs an unreduced (F, G): reproduce what the descent hands
    # it at the top level.
    Fp, Gp = ntru_solve(field_norm(fb), field_norm(gb); q = p.q)
    Fraw = karamul(lift(Fp), galois_conjugate(gb))
    Graw = karamul(lift(Gp), galois_conjugate(fb))
    measure!("babai_reduce (top level only)", 2,
             _ -> Falcon.babai_reduce(fb, gb, Fraw, Graw); iters = 5)
    @printf("# top-level (F,G) before reduction: %d bits max; after: %d bits max\n",
            maximum(Falcon.bitsize, Fraw),
            maximum(Falcon.bitsize, first(Falcon.babai_reduce(fb, gb, Fraw, Graw))))

    # ---- the same multiply, three ways -------------------------------------
    # karamul is the workhorse of the descent and it is written over BigInt.
    # How much of its cost is Karatsuba, and how much is BigInt itself?
    fi = Int.(f); gi = Int.(g)
    measure!("karamul on BigInt (again, for contrast)", 2, _ -> karamul(fb, gb); iters = 20)
    measure!("polymul on Int (schoolbook, same ring)", 2, _ -> polymul(fi, gi); iters = 20)
    measure!("polymulq on Int (schoolbook mod q)", 2, _ -> Falcon.polymulq(fi, gi); iters = 20)

    end  # !skip_keygen

    F, G = ntru_solve(f, g; q = p.q)
    measure!("expand_privkey (FFT basis + ffLDL tree)", 1,
             _ -> expand_privkey(Int.(f), Int.(g), Int.(F), Int.(G), p); iters = 10)

    sk = expand_privkey(Int.(f), Int.(g), Int.(F), Int.(G), p)
    pk = public_key(sk)
    print_rows(); empty!(ROWS)

    # ---- signing ------------------------------------------------------------
    println("# signing")
    println()
    msg = collect(b"falcon-jl profile message")
    ssrc = fresh_source("falcon-jl/profile/sign", 1 << 20)
    measure_both!("falcon_sign (whole)", s -> falcon_sign(sk, msg, bytesource_of(s));
                  iters = 200, setup = ssrc)

    # How much of signing is the `randombytes` interface itself?  Every call
    # returns a freshly allocated `Vector{UInt8}`, and the sampler makes a great
    # many small ones -- 9 bytes for the base sampler, 1 for the sign bit, 1 per
    # rejection round inside berexp.  Count them, then price that many calls at
    # the same size mix, so the figure is the interface's cost and not a guess.
    cnt = CountingSource(ssrc())
    falcon_sign(sk, msg, cnt)
    @printf("# one signature makes %d randombytes calls for %d bytes (mean %.1f per call)\n",
            cnt.calls, cnt.bytes, cnt.bytes / cnt.calls)
    ncalls = cnt.calls
    measure!("randombytes interface alone (that many calls)", 2,
             s -> (for i in 1:ncalls; s(i % 3 == 1 ? 9 : 1); end);
             iters = 20, setup = fresh_source("falcon-jl/profile/src", 1 << 20))

    salt = shake256(codeunits("salt"), SALT_LEN)
    measure!("hash_to_point", 2, _ -> hash_to_point(msg, salt, n; q = p.q); iters = 50)
    point = hash_to_point(msg, salt, n; q = p.q)
    measure!("sample_preimage (ffSampling etc.)", 2,
             s -> Falcon.sample_preimage(sk, point, bytesource_of(s));
             iters = 20, setup = ssrc)
    s1, s2 = Falcon.sample_preimage(sk, point, bytesource_of(ssrc()))
    measure!("encode_signature", 2,
             _ -> encode_signature(salt, s2, p.logn, p.sig_bytes); iters = 50)
    measure!("sqnorm(s1, s2) (BigInt)", 2, _ -> sqnorm(s1, s2); iters = 50)

    # inside sample_preimage
    measure!("fft (one forward transform)", 2, _ -> fft(Float64.(point)); iters = 50)
    tf = fft(Float64.(point))
    measure!("ifft (one inverse transform)", 2, _ -> ifft(tf); iters = 50)
    measure!("mul_fft (one FFT-domain product)", 2, _ -> mul_fft(tf, tf); iters = 50)
    # samplerz's precondition is 1 < sigmin < sigma < 1.8205, and the sigma it
    # is actually called with is a *leaf* sigma of the normalised tree -- NOT
    # the key's sigma (165.7), which is a width in the lattice, not in Z.
    # Passing p.sigma here throws inside approxexp; that is the sampler
    # correctly refusing an out-of-range input, not a defect.
    leafsig = first(leaf_sigmas(sk.tree))
    measure!("samplerz x n (the sampler alone)", 2,
             s -> (for _ in 1:n; samplerz(0.5, leafsig, p.sigma_min, s); end);
             iters = 10, setup = ssrc)
    print_rows(); empty!(ROWS)

    # ---- verification -------------------------------------------------------
    println("# verification")
    println()
    sig = falcon_sign(sk, msg, bytesource_of(ssrc()))
    measure_both!("falcon_verify (whole)", _ -> falcon_verify(pk, msg, sig); iters = 400)
    measure!("decode_signature", 2, _ -> decode_signature(sig); iters = 200)
    _, vsalt, vs2 = decode_signature(sig)
    measure!("hash_to_point", 2, _ -> hash_to_point(msg, vsalt, n; q = p.q); iters = 200)
    s2q = Int[mod(c, p.q) for c in vs2]
    measure!("polymulq (schoolbook)", 2, _ -> Falcon.polymulq(s2q, pk.h, p.q); iters = 200)
    measure!("polymulq_ntt", 2, _ -> polymulq_ntt(s2q, pk.h); iters = 200)
    measure!("ntt (one forward transform)", 2, _ -> ntt(s2q); iters = 200)
    vprod = Falcon.polymulq(s2q, pk.h, p.q)
    vs1 = Falcon.centered(Falcon.polysubq(hash_to_point(msg, vsalt, n; q = p.q), vprod, p.q), p.q)
    measure!("sqnorm (BigInt)", 2, _ -> sqnorm(vs1, vs2); iters = 200)
    measure!("sqnorm equivalent in Int64", 2,
             _ -> (a = 0; for c in vs1; a += c * c; end; for c in vs2; a += c * c; end; a);
             iters = 200)
    print_rows()

    return nothing
end

main(ARGS)
