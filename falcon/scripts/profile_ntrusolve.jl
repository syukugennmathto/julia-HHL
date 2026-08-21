#!/usr/bin/env julia
#
# Where inside ntru_solve does the time go?
#
#     julia --project=falcon falcon/scripts/profile_ntrusolve.jl [logn]
#
# scripts/profile_stages.jl establishes that ntru_solve is essentially all of
# key generation.  That is not yet an answer: "the descent is slow" could mean
# either of two very different things.
#
#   (a) The coefficients genuinely get huge near the bottom of the descent, so
#       the work is irreducibly big-integer work and the only real fix is to
#       change representation the way the C reference does (a residue number
#       system over 31-bit limbs, i.e. machine words throughout).
#
#   (b) Most of the time is spent at the *top* levels, where the coefficients
#       are small enough to fit a machine word and BigInt is being paid for
#       nothing.
#
# These have opposite implications for what to do next, and the two are
# distinguished by measurement, not by argument.  So this walks the same
# recursion `ntru_solve` walks, timing each level and recording how many bits
# the coefficients actually occupy there.
#
# The recursion below MUST stay a faithful copy of src/ntrugen.jl's.  It calls
# the library's own field_norm / karamul / lift / galois_conjugate /
# babai_reduce, so only the control flow is duplicated -- but if ntru_solve
# changes, this changes with it or it starts measuring a different algorithm.
# The check at the end (that the instrumented descent produces a solution
# satisfying the NTRU equation) is what keeps that honest.

using Falcon
using Printf

const Q_DEFAULT = 12289

mutable struct Level
    n::Int
    calls::Int
    t_field_norm::Float64
    t_karamul_lift::Float64
    t_babai::Float64
    alloc_field_norm::Int
    alloc_karamul_lift::Int
    alloc_babai::Int
    bits_in::Int          # max bitsize of f, g arriving at this level
    bits_raw::Int         # max bitsize of (F, G) before Babai reduction
    bits_out::Int         # ... and after
end

Level(n) = Level(n, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0)

const LEVELS = Dict{Int,Level}()
level(n) = get!(() -> Level(n), LEVELS, n)

"Run `f`, returning (value, elapsed_ms, bytes_allocated)."
function tracked(f)
    st = Base.gc_num()
    t0 = time_ns()
    v = f()
    dt = (time_ns() - t0) / 1e6
    d = Base.GC_Diff(Base.gc_num(), st)
    return (v, dt, d.allocd)
end

maxbits(v) = isempty(v) ? 0 : maximum(Falcon.bitsize, v)

function solve_instrumented(f::Vector{BigInt}, g::Vector{BigInt}, q::Integer)
    n = length(f)
    L = level(n)
    L.calls += 1
    L.bits_in = max(L.bits_in, maxbits(f), maxbits(g))

    if n == 1
        d, u, v = Falcon.xgcd_floor(f[1], g[1])
        d == 1 || throw(Falcon.NTRUSolveFailure("gcd != 1 at the bottom"))
        return (BigInt[-q * v], BigInt[q * u])
    end

    (fp, t1, a1) = tracked(() -> Falcon.field_norm(f))
    (gp, t2, a2) = tracked(() -> Falcon.field_norm(g))
    L.t_field_norm += t1 + t2
    L.alloc_field_norm += a1 + a2

    Fp, Gp = solve_instrumented(fp, gp, q)

    (FG, t3, a3) = tracked(function ()
        F = Falcon.karamul(Falcon.lift(Fp), Falcon.galois_conjugate(g))
        G = Falcon.karamul(Falcon.lift(Gp), Falcon.galois_conjugate(f))
        return (F, G)
    end)
    L.t_karamul_lift += t3
    L.alloc_karamul_lift += a3
    L.bits_raw = max(L.bits_raw, maxbits(FG[1]), maxbits(FG[2]))

    (out, t4, a4) = tracked(() -> Falcon.babai_reduce(f, g, FG[1], FG[2]))
    L.t_babai += t4
    L.alloc_babai += a4
    L.bits_out = max(L.bits_out, maxbits(out[1]), maxbits(out[2]))

    return out
end

function main(args)
    logn = length(args) >= 1 ? parse(Int, args[1]) : 9
    n = 1 << logn
    q = Q_DEFAULT

    println("# profile_ntrusolve.jl  julia=", VERSION, "  n=", n)
    println()

    # An (f, g) that ntru_solve accepts.  Found with the library's own
    # acceptance test so that the descent measured is a real one.
    src = ReplayBytes(shake256(codeunits("falcon-jl/profile/ntrusolve"), 1 << 23))
    local f, g
    attempts = 0
    while true
        attempts += 1
        f = gen_poly(n, src)
        g = gen_poly(n, src)
        gs_norm(Float64.(f), Float64.(g); q = q) > gram_schmidt_quality()^2 * q && continue
        is_invertible_zq(Int.(f)) || continue
        try
            ntru_solve(f, g; q = q)
            break
        catch e
            e isa Falcon.NTRUSolveFailure || rethrow()
        end
    end
    println("# accepted (f, g) after ", attempts, " attempt(s); ",
            "max bitsize of f: ", maxbits(BigInt.(f)))
    println()

    # Warm up: compile everything before the timed descent.
    solve_instrumented(BigInt.(f), BigInt.(g), q)
    empty!(LEVELS)

    (FG, total, alloc) = tracked(() -> solve_instrumented(BigInt.(f), BigInt.(g), q))
    F, G = FG

    Falcon.ntru_equation_holds(BigInt.(f), BigInt.(g), F, G; q = q) ||
        error("the instrumented descent did not solve the NTRU equation -- " *
              "it has drifted from src/ntrugen.jl and is measuring something else")

    @printf("%6s %6s %10s %10s %10s %10s %8s %8s %8s\n",
            "n", "calls", "fieldnrm", "karamul", "babai", "level ms", "bits_in",
            "raw", "out")
    println("-"^92)
    tot_fn = tot_km = tot_bb = 0.0
    for k in sort(collect(keys(LEVELS)), rev = true)
        L = LEVELS[k]
        lv = L.t_field_norm + L.t_karamul_lift + L.t_babai
        tot_fn += L.t_field_norm; tot_km += L.t_karamul_lift; tot_bb += L.t_babai
        @printf("%6d %6d %10.2f %10.2f %10.2f %10.2f %8d %8d %8d\n",
                L.n, L.calls, L.t_field_norm, L.t_karamul_lift, L.t_babai, lv,
                L.bits_in, L.bits_raw, L.bits_out)
    end
    println("-"^92)
    @printf("%6s %6s %10.2f %10.2f %10.2f %10.2f\n",
            "total", "", tot_fn, tot_km, tot_bb, tot_fn + tot_km + tot_bb)
    @printf("\nwhole descent: %.1f ms, %.2f GB allocated\n", total, alloc / 2^30)

    println("\n# allocation by stage (GB)")
    @printf("%6s %12s %12s %12s\n", "n", "fieldnrm", "karamul", "babai")
    for k in sort(collect(keys(LEVELS)), rev = true)
        L = LEVELS[k]
        @printf("%6d %12.3f %12.3f %12.3f\n", L.n,
                L.alloc_field_norm / 2^30, L.alloc_karamul_lift / 2^30,
                L.alloc_babai / 2^30)
    end

    # --- the decisive contrast ---------------------------------------------
    # How much of a BigInt operation's cost is the arbitrary precision, and how
    # much is the mere fact that BigInt is a heap object?  Compare the same
    # arithmetic on values that fit a machine word.
    println("\n# BigInt versus machine integers, on values that fit either")
    small = rand(Int64(-1000):Int64(1000), 1 << 16)
    sb = BigInt.(small)
    for (name, v) in (("Int64", small), ("Int128", Int128.(small)), ("BigInt", sb))
        f2 = () -> (a = zero(eltype(v)); @inbounds for i in eachindex(v); a += v[i] * v[i]; end; a)
        f2()
        (_, t, al) = tracked(f2)
        @printf("  sum of squares over %d %-7s %8.3f ms  %8.2f MB\n",
                length(v), name * ":", t, al / 2^20)
    end
    return nothing
end

main(ARGS)
