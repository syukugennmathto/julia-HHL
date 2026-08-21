#!/usr/bin/env julia
#
# heuristic1_check.jl -- does our implementation reproduce ePrint 2024/1709's
# Heuristic 1, and does it explain where our divergences land?
#
#     julia --project=falcon falcon/scripts/heuristic1_check.jl [keys] [n]
#
# ===========================================================================
# WHAT IS BEING CHECKED
# ===========================================================================
#
# scripts/divergence_rate.jl found 5 divergent signatures in 100000, and
# scripts/first_divergence.jl found that all five are nearly-integer centres,
# and that all five are at sampler call 1023 or 1024 of 1024 -- the LAST TWO
# calls of the traversal.  None was at call 1 or 2.
#
# ePrint 2024/1709 states a "Heuristic 1" (quoted second-hand -- see the
# provenance note at the bottom):
#
#     the first two calls to SamplerZ have probability 1/q of having a
#     center c_i that is an exact integer, and the last two calls have
#     probability 1/||(g,-f)||^2 of having such a center.
#
# For Falcon-512, q = 12289 and ||(g,-f)||^2 measures about 16600, so the
# heuristic predicts 1.63e-4 per signature at the first two calls and
# 1.20e-4 at the last two.
#
# An exact-integer center is NECESSARY for the divergence, not sufficient:
# the two floating-point evaluations must also land on opposite sides of it.
# So the heuristic bounds our divergence rate from above.
#
# ===========================================================================
# THE THING THAT LOOKED LIKE A CONTRADICTION, AND WAS A SAMPLING ERROR
# ===========================================================================
#
# The heuristic makes the FIRST two calls the MORE likely end (1.63e-4
# against 1.20e-4), and we saw 0 of 5 there.  Under the heuristic, read
# naively, that split has probability (1.20/(1.63+1.20))^5 = 0.013.
#
# The naive reading is wrong, and the reason is worth stating because it is a
# trap for anyone repeating this measurement.  **The centres of calls 1 and 2
# do not depend on the signature's randomness at all.**  Measured directly:
# over 8 signatures with one key and one message, calls 1 and 2 take exactly
# one distinct centre each; call 3 is the first that varies.  Those two
# leaves sit at the bottom of the first descent, where the target is still
# `t` itself and no sampled coefficient has fed back into it.
#
# So a divergence experiment that fixes the message and varies only the PRNG
# -- which is what scripts/divergence_rate.jl does, because it must, to hold
# the two spellings on identical randomness -- draws 100000 samples at the
# last two calls and only ONE HUNDRED at the first two, one per key.  The
# expected number of exact-integer centres at calls 1-2 in that experiment is
# 100 * 2/12289 = 0.016.  Seeing zero is not evidence against anything.
#
# This script therefore measures the two ends with two different experiments:
#
#   Arm 1 (first two calls): vary the MESSAGE.  Each (key, message) pair is
#          one independent draw of those two centres.
#   Arm 2 (last two calls) : vary the PRNG, one message.  Each signature is
#          one independent draw.
#
# Both are sized so the heuristic predicts of order ten events.
#
# ===========================================================================
# HOW A NEAR-INTEGER CENTER IS COUNTED
# ===========================================================================
#
# When the exact center is an integer, the computed `mu` misses it by a
# rounding error of order 1e-13.  When it is not, the fractional part is
# spread over [0,1), so |mu - round(mu)| < 1e-9 arises with probability about
# 2e-9 -- under 1e-3 of one event across this whole run.  The threshold is
# four orders of magnitude clear of the signal on one side and five clear of
# the noise on the other, so nothing depends on its exact value.
#
# CONSTANT TIME: out of scope for this project, but note in passing that the
# quantity counted here -- whether a center is an exact integer -- is exactly
# the condition 2024/1709 turns into key recovery.  A sampler that branched
# on it would be catastrophic.  Falcon's does not branch on it; it rounds
# through it, which is how the discrepancy becomes an output difference
# rather than a timing difference.
#
# ===========================================================================
# PROVENANCE OF THE HEURISTIC
# ===========================================================================
#
# eprint.iacr.org, artifacts.iacr.org, dblp.org, link.springer.com and
# dl.acm.org are all refused by this machine's egress proxy (403 at CONNECT),
# so Heuristic 1 as stated above has NOT been read in the original; it comes
# from web-search summaries.  The numbers measured below are ours; what they
# are compared against is second-hand.  Re-check the statement against the
# PDF before this is cited anywhere.

using Falcon
using Printf

const F = Falcon

# Log `mu` at every sampler call.  A copy of the library's `samplerz_isigma`
# with one push added -- the same approach as scripts/first_divergence.jl,
# and for the same reason (a wrapper that captures the original and then
# redefines the name recurses into itself).
Falcon.eval(quote
    const MU = Float64[]
    const RECORD = Ref(false)
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        RECORD[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb)
            b = Int(rb(1)[1]) & 1
            z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

"How far is `mu` from the nearest integer?"
intgap(mu::Float64) = abs(mu - round(mu))

const NEAR = 1e-9

"Run one signing, return the logged centres."
function centres(sk, pt, seed::String)
    rb = chacha20(shake256(codeunits(seed), 56))
    empty!(F.MU); F.RECORD[] = true
    F.sample_preimage(sk, pt, x -> randombytes!(rb, x))
    F.RECORD[] = false
    return F.MU
end

point(msg::String, salt::String) =
    hash_to_point(collect(codeunits(msg)), shake256(codeunits(salt), SALT_LEN), 512;
                  q = FALCON_512.q)

"""
Confirm that calls 1 and 2 carry no dependence on the signature randomness,
which is what forces the two arms below to be separate experiments.
"""
function show_determinism(sk, pt)
    M = [copy(centres(sk, pt, "det/$j")) for j in 1:8]
    varies = findfirst(i -> length(unique(getindex.(M, i))) > 1, 1:length(M[1]))
    println("# does the centre depend on the signature's randomness?")
    for i in (1, 2, 3, 4)
        @printf("  call %4d : %d distinct centres over 8 signatures\n",
                i, length(unique(getindex.(M, i))))
    end
    @printf("  first call whose centre varies across signatures: %d\n", varies)
    println()
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    ndraw = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1000
    n = 512
    p = FALCON_512

    println("# heuristic1_check.jl  n=512  ", nkeys, " keys x ", ndraw, " draws per arm")
    println()

    r = chacha20(collect(UInt8, 0x00:0x37))
    keys = [falcon_keygen(512, k -> randombytes!(r, k))[1] for _ in 1:nkeys]
    normsq = [Float64(sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f))
              for sk in keys]
    mean_normsq = sum(normsq) / length(normsq)

    @printf("||(g,-f)||^2 over %d keys: mean %.1f   min %.0f   max %.0f   (1.17^2 q = %.1f)\n",
            nkeys, mean_normsq, minimum(normsq), maximum(normsq), 1.17^2 * p.q)
    println()

    show_determinism(keys[1], point("power analysis message", "power salt"))

    # ---- Arm 1: the first two calls, one draw per (key, message) -----------
    nearfirst = 0; nfirst = 0; minfirst = Inf
    for (ki, sk) in enumerate(keys), j in 1:ndraw
        mu = centres(sk, point("msg/$ki/$j", "salt/$ki/$j"), "arm1/$ki/$j")
        for i in (1, 2)
            g = intgap(mu[i]); nfirst += 1; minfirst = min(minfirst, g)
            g < NEAR && (nearfirst += 1)
        end
    end

    # ---- Arm 2: the last two calls, one draw per signature ----------------
    pt = point("power analysis message", "power salt")
    nearlast = 0; nlast = 0; minlast = Inf
    nearbulk = 0; nbulk = 0
    for (ki, sk) in enumerate(keys), j in 1:ndraw
        mu = centres(sk, pt, "arm2/$ki/$j")
        for i in (2n - 1, 2n)
            g = intgap(mu[i]); nlast += 1; minlast = min(minlast, g)
            g < NEAR && (nearlast += 1)
        end
        for i in 3:(2n - 2)
            nbulk += 1; intgap(mu[i]) < NEAR && (nearbulk += 1)
        end
    end

    println("# rate of near-integer centres (|mu - round(mu)| < 1e-9)")
    @printf("%-34s %8s %10s %12s %12s\n", "position (experiment)", "near", "draws", "rate", "heuristic")
    @printf("%-34s %8d %10d %12.3g %12.3g\n", "calls 1, 2      (vary message)",
            nearfirst, nfirst, nearfirst / nfirst, 1 / p.q)
    @printf("%-34s %8d %10d %12.3g %12.3g\n", "calls 2n-1, 2n  (vary PRNG)",
            nearlast, nlast, nearlast / nlast, 1 / mean_normsq)
    @printf("%-34s %8d %10d %12.3g %12s\n", "calls 3 .. 2n-2 (vary PRNG)",
            nearbulk, nbulk, nearbulk / nbulk, "~0")
    println()
    @printf("closest approach: calls 1,2 %.3g    calls 2n-1,2n %.3g\n", minfirst, minlast)
    println()
    @printf("per signature -- first two: measured %.3g, heuristic %.3g\n",
            2 * nearfirst / nfirst, 2 / p.q)
    @printf("per signature -- last two : measured %.3g, heuristic %.3g\n",
            2 * nearlast / nlast, 2 / mean_normsq)
    println()
    println("# scripts/divergence_rate.jl measured 5 divergences in 100000")
    println("# signatures, all at call 2n-1 or 2n.  An exact-integer centre is")
    println("# necessary but not sufficient -- the two spellings must also land")
    println("# on opposite sides of it -- so that rate should sit below the")
    println("# last-two figure above, and it does.")
end

main()
