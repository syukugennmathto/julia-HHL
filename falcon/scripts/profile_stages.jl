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
    alloc_mb::Float64
    gc_pct::Float64
    n::Int
end

const ROWS = Row[]

"""
Time `f` `iters` times and record the median wall time, the allocation of a
single call, and the fraction of wall time spent in GC.

`setup` runs before each timed call and its result is passed to `f`; use it to
hand each run an identical fresh input (a rewound byte source, say).
"""
function measure!(label, depth, f; iters = 20, setup = () -> nothing)
    f(setup())                                   # compile, and warm any cache
    ts = Float64[]
    gcs = Float64[]
    local alloc
    for _ in 1:iters
        arg = setup()
        GC.gc(false)                             # start from a comparable heap
        st = Base.gc_num()
        t0 = time_ns()
        f(arg)
        dt = (time_ns() - t0) / 1e6
        d = Base.GC_Diff(Base.gc_num(), st)
        push!(ts, dt)
        push!(gcs, d.total_time / 1e6)
        alloc = d.allocd
    end
    p = sortperm(ts)
    mid = p[length(p) ÷ 2 + 1]
    push!(ROWS, Row(label, depth, ts[mid], alloc / 2^20,
                    ts[mid] > 0 ? 100 * gcs[mid] / ts[mid] : 0.0, iters))
    return nothing
end

function print_rows()
    @printf("%-44s %12s %12s %7s %6s\n", "stage", "median ms", "alloc MB", "gc%", "iters")
    println("-"^86)
    for r in ROWS
        @printf("%-44s %12.4f %12.3f %7.1f %6d\n",
                "  "^r.depth * r.label, r.ms, r.alloc_mb, r.gc_pct, r.n)
    end
    println()
end

"A deterministic byte source that can be rewound, so every run does equal work."
function fresh_source(seed::AbstractString, nbytes::Integer)
    return () -> ReplayBytes(shake256(codeunits(seed), nbytes))
end

# ---------------------------------------------------------------------------
# the stages
# ---------------------------------------------------------------------------

function main(args)
    logn = length(args) >= 1 ? parse(Int, args[1]) : 9
    n = 1 << logn
    p = params(n)

    println("# profile_stages.jl  julia=", VERSION, "  n=", n)
    println()

    # ---- key generation -----------------------------------------------------
    # gen_poly draws 4096 samplerz values regardless of n, so its cost is
    # constant in n; ntru_solve's is not.  Timing them apart is the whole point.
    src = fresh_source("falcon-jl/profile/keygen", 1 << 22)

    measure!("gen_poly (4096 samplerz draws)", 1, s -> gen_poly(n, bytesource_of(s));
             iters = 5, setup = src)

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
    measure!("falcon_sign (whole)", 1, s -> falcon_sign(sk, msg, bytesource_of(s));
             iters = 20, setup = ssrc)

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
    measure!("samplerz x n (the sampler alone)", 2,
             s -> (bs = bytesource_of(s); for _ in 1:n; samplerz(0.5, p.sigma, p.sigma_min, bs); end);
             iters = 10, setup = ssrc)
    print_rows(); empty!(ROWS)

    # ---- verification -------------------------------------------------------
    println("# verification")
    println()
    sig = falcon_sign(sk, msg, bytesource_of(ssrc()))
    measure!("falcon_verify (whole)", 1, _ -> falcon_verify(pk, msg, sig); iters = 200)
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
