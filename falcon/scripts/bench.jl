#!/usr/bin/env julia
#
# Time this implementation's keygen / sign / verify, in the same shape as
# scripts/cref_bench.c times the C reference's, so the two can be put side by
# side.
#
# Usage:
#     julia --project=falcon falcon/scripts/bench.jl [logn] [keygen] [sign] [verify]
#
# Output: one "op n median_ms mean_ms min_ms iters" line per operation, the
# same format cref_bench.c emits, so `scripts/bench_compare.jl` can read both.
#
# What is and is not being measured
# ---------------------------------
# This is a *specification-conformant* implementation, not an optimised one.
# Nothing here is vectorised, `polymul` is schoolbook, `ntru_solve` carries its
# intermediate values in `BigInt`, and no attempt has been made to avoid
# allocation.  The C reference, by contrast, has had years of work put into it
# by its author, keeps everything in a single caller-supplied scratch buffer,
# and does its big-integer work in a hand-rolled RNS with 31-bit limbs.  A
# large ratio is the expected outcome and is not by itself a bug.
#
# What the ratio is *useful* for is finding the places where the gap is much
# bigger than the general factor -- those are real algorithmic differences
# rather than constant factors, and they are what a later optimisation pass
# should go after.
#
# Julia-specific care:
#
#   - The first call to anything is compilation, not execution.  Every
#     operation is run once before timing starts.  Forgetting this is the
#     classic way to publish a Julia benchmark that is wrong by 100x.
#
#   - The GC is not disabled.  Allocation *is* part of this implementation's
#     cost, and hiding it would flatter the numbers.
#
#   - `median` is quoted rather than `mean` for the same reason as in the C
#     driver: keygen retries, so its distribution has a long right tail.

using Falcon
using Printf
using Statistics

"A deterministic byte source, so two runs of this script are comparable."
function make_rng(seed::AbstractString)
    buf = UInt8[]
    pos = Ref(0)
    ctr = Ref(0)
    return function (k::Integer)
        while pos[] + k > length(buf)
            append!(buf, shake256(codeunits(seed * "/" * string(ctr[])), 1 << 16))
            ctr[] += 1
        end
        out = buf[(pos[] + 1):(pos[] + k)]
        pos[] += k
        return out
    end
end

function report(op::AbstractString, n::Integer, t::Vector{Float64})
    @printf("%s %d %.6f %.6f %.6f %d\n",
            op, n, median(t), mean(t), minimum(t), length(t))
    flush(stdout)
end

function main(args)
    logn = length(args) >= 1 ? parse(Int, args[1]) : 9
    nkg  = length(args) >= 2 ? parse(Int, args[2]) : 5
    nsg  = length(args) >= 3 ? parse(Int, args[3]) : 50
    nvf  = length(args) >= 4 ? parse(Int, args[4]) : 200
    n = 1 << logn

    println("# bench.jl julia=", VERSION, " n=", n)
    println("# op n median_ms mean_ms min_ms iters")

    rng = make_rng("falcon-jl/bench")
    msg = "falcon-jl benchmark message"

    # --- warm-up: compile everything before any clock starts -----------------
    let (sk0, pk0) = falcon_keygen(n, rng)
        sig0 = falcon_sign(sk0, msg, rng)
        falcon_verify(pk0, msg, sig0)
        expand_privkey(sk0.f, sk0.g, sk0.F, sk0.G, sk0.params)
    end

    # --- keygen --------------------------------------------------------------
    t = Float64[]
    local sk, pk
    for _ in 1:nkg
        t0 = time_ns()
        sk, pk = falcon_keygen(n, rng)
        push!(t, (time_ns() - t0) / 1e6)
    end
    report("keygen", n, t)

    # --- expand_privkey ------------------------------------------------------
    # Split out for the same reason the C driver splits it out: our
    # `falcon_sign` takes an already-expanded key, so `sign` below is the
    # counterpart of C's `sign_tree`, and `expand + sign` is the counterpart
    # of C's `sign_dyn`.
    t = Float64[]
    for _ in 1:min(nsg, 20)
        t0 = time_ns()
        expand_privkey(sk.f, sk.g, sk.F, sk.G, sk.params)
        push!(t, (time_ns() - t0) / 1e6)
    end
    report("expand_privkey", n, t)

    # --- sign ----------------------------------------------------------------
    t = Float64[]
    local sig
    for _ in 1:nsg
        t0 = time_ns()
        sig = falcon_sign(sk, msg, rng)
        push!(t, (time_ns() - t0) / 1e6)
    end
    report("sign_tree", n, t)

    # --- verify --------------------------------------------------------------
    # Measured **from bytes**, i.e. `pubkey_from_bytes` + `falcon_verify`, because
    # C's `falcon_verify` takes the encoded public key and decodes it on every
    # call.  Timing our parsed-key form against that would be comparing two
    # different amounts of work and would flatter us by the cost of a decode.
    pkb = pubkey_bytes(pk)
    @assert falcon_verify(pubkey_from_bytes(pkb), msg, sig)
    t = Float64[]
    for _ in 1:nvf
        t0 = time_ns()
        falcon_verify(pubkey_from_bytes(pkb), msg, sig)
        push!(t, (time_ns() - t0) / 1e6)
    end
    report("verify", n, t)

    return nothing
end

main(ARGS)
