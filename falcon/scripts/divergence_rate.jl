#!/usr/bin/env julia
#
# divergence_rate.jl -- do the specification's underdetermined choices reach
# the signature, and how often?
#
#     julia --project=falcon falcon/scripts/divergence_rate.jl [keys] [sigs_per_key]
#
# ===========================================================================
# WHAT THIS MEASURES, AND WHY THE SAMPLE SIZE IS THE WHOLE POINT
# ===========================================================================
#
# Two implementations of Falcon can both follow the specification and still
# compute different intermediate values, because the specification does not
# determine them.  This project found two independent instances
# (docs/debug_log.md #048, #053):
#
#   A. HOW AN OPERATION IS SPELLED.  Complex division, the `D11` entry of
#      LDL*, and the bottom two levels of ffSampling are written one way in
#      the specification and the Python reference, and an algebraically
#      identical other way in the C reference.  These perturb the sampler's
#      CENTRE, `mu`.
#
#   B. A CONSTANT THAT CANNOT BE DERIVED.  `fpr_inv_sigma[logn]` is not the
#      correctly rounded reciprocal of the sigma Table 3.3 publishes
#      (scripts/inv_sigma_audit.jl).  This perturbs the sampler's WIDTH.
#
# The first measurement of A used 480 signatures and found zero divergences.
# That was **under-powered and the conclusion drawn from it was wrong**.
# ePrint 2024/1709 reports that Falcon's sampler diverges at nearly-integer
# centres with probability on the order of 1e-4 per signature; 480 signatures
# cannot distinguish 1e-4 from 0.  At 100000 signatures the answer changes:
#
#   A: 5 of 100000 signatures diverge (5e-5).  When one does, about 470 of
#      512 coefficients differ -- the byte stream desynchronises and the rest
#      of the signature is unrelated.
#   B: 0 of 100000 (95% upper bound 3e-5), even though 444 of 512 tree leaves
#      differ.
#
# ===========================================================================
# THE MECHANISM, AND WHY A AND B DIFFER
# ===========================================================================
#
# Every one of the five divergences under A is a **nearly-integer centre**:
#
#     call 1023   mu = 221.00000000000006   vs   220.99999999999997
#     call 1024   mu = 436.99999999999994   vs   437.00000000000011
#     call 1024   mu = 444.0                vs   443.99999999999977
#     call 1023   mu = -180.99999999999997  vs  -181.0000000000002
#     call 1023   mu =  60.000000000000007  vs   59.999999999999957
#
# `SamplerZ` begins `s = floor(mu)`.  `floor` is discontinuous, so a
# discrepancy of 1e-13 in `mu` becomes a difference of 1 in `s` whenever the
# two land on opposite sides of an integer.  That changes `z`, changes how
# many bytes `BerExp` consumes, and desynchronises the PRNG.
#
# All five are at call 1023 or 1024 of 1024 -- positions 2n-2 and 2n-1.  That
# is exactly where ePrint 2024/1709 predicts the centres concentrate near
# integers.  This is their mechanism, reproduced from a different cause.
#
# **B does not fire because it perturbs the width, not the centre.**  The
# width enters through `dss = 1/(2*sigma^2)` and `ccs = sigma_min/sigma`,
# which feed `BerExp`'s comparison of a fixed-point exponential against random
# bytes.  That comparison is smooth: a relative perturbation of 1e-16 flips it
# with probability of the same order, not with the ~1e-8 per draw that a floor
# straddle achieves.
#
# So the boundary is not "big differences propagate and small ones do not".
# It is **which quantity the difference reaches**: the centre goes through a
# floor, the width does not.
#
# ===========================================================================
# WHAT THIS ADDS TO ePrint 2024/1709
# ===========================================================================
#
# They perturbed by changing the BUILD: FMA on or off, fpemu against fpnative
# against AVX2, dynamic against tree signing -- different compilations of the
# same source.  This shows the perturbation does not need a different build.
# **Two implementations that both follow the specification are enough**,
# because the specification does not say which of two algebraically equal
# formulas to use.  For a standard that intends to require bit-exact KATs and
# forbids the derandomized settings where this becomes key recovery, that is
# the difference between "implementations should agree" and "the text must say
# which formula".

using Falcon
using Printf
using Statistics

const F = Falcon

"""
Rate at which the two spellings of the operations (A) produce different
signatures, given identical key, message and PRNG state.
"""
function rate_spelling(nkeys::Int, nsig::Int)
    p = FALCON_512
    r = chacha20(collect(UInt8, 0x00:0x37))
    keys = [falcon_keygen(512, k -> randombytes!(r, k))[1] for _ in 1:nkeys]
    msg = collect(codeunits("power analysis message"))
    pt = hash_to_point(msg, shake256(codeunits("power salt"), SALT_LEN), 512; q = p.q)

    n = 0; ndiff = 0; ncoef = 0; ncoefdiff = 0
    divergent = Tuple{Int,Int,Int}[]
    for (ki, sk) in enumerate(keys), j in 1:nsig
        st = shake256(codeunits("prng/$ki/$j"), 56)
        r1 = chacha20(st); r2 = chacha20(st)
        _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
        _, b = with_spec_ffsampling() do
            F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
        end
        n += 1
        d = count(Int.(a) .!= Int.(b))
        ncoef += length(a); ncoefdiff += d
        d > 0 && (ndiff += 1; push!(divergent, (ki, j, d)))
    end
    return (n, ndiff, ncoef, ncoefdiff, divergent)
end

"""
Rate at which using `1/sigma` instead of the reference's `fpr_inv_sigma`
table (B) produces different signatures.
"""
function rate_constant(nkeys::Int, nsig::Int)
    p = FALCON_512
    recip = 1 / p.sigma
    r = chacha20(collect(UInt8, 0x00:0x37))
    msg = collect(codeunits("power analysis message"))
    pt = hash_to_point(msg, shake256(codeunits("power salt"), SALT_LEN), 512; q = p.q)

    n = 0; ndiff = 0; ncoef = 0; ncoefdiff = 0; nleaf = 0
    for ki in 1:nkeys
        sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
        t2 = F.ffldl_fft(F.gram_fft(sk.B0_fft))
        F.normalize_tree!(t2, p.sigma)
        patch!(t) = t isa F.FFLDLNode ? (patch!(t.left); patch!(t.right)) :
                    (t.isigma = sqrt(real(t.value[1])) * recip)
        patch!(t2)
        nleaf += count(zip(leafisig(sk.tree), leafisig(t2))) do (x, y)
            reinterpret(UInt64, x) != reinterpret(UInt64, y)
        end
        sk2 = F.FalconPrivateKey(p, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t2)
        for j in 1:nsig
            st = shake256(codeunits("prngB/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            _, a = F.sample_preimage(sk,  pt, x -> randombytes!(r1, x))
            _, b = F.sample_preimage(sk2, pt, x -> randombytes!(r2, x))
            n += 1
            d = count(Int.(a) .!= Int.(b))
            ncoef += length(a); ncoefdiff += d
            d > 0 && (ndiff += 1)
        end
    end
    return (n, ndiff, ncoef, ncoefdiff, nleaf / nkeys)
end

leafisig(t::F.FFLDLLeaf) = [t.isigma]
leafisig(t::F.FFLDLNode) = vcat(leafisig(t.left), leafisig(t.right))

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1000
    println("# divergence_rate.jl  n=512  ", nkeys, " keys x ", nsig, " signatures")
    println()

    n, nd, nc, ncd, div = rate_spelling(nkeys, nsig)
    println("## A -- how the operations are spelled (perturbs the CENTRE)")
    @printf("  signatures differing   : %d of %d   (rate %.3g)\n", nd, n, nd / n)
    @printf("  coefficients differing : %d of %d\n", ncd, nc)
    if nd == 0
        @printf("  zero events -> 95%% upper bound on the rate: %.3g\n", 3 / n)
    else
        @printf("  when one differs, %.0f of 512 coefficients do\n", ncd / nd)
        println("  divergent cases (key, signature, coefficients):")
        for d in div
            @printf("    %s\n", string(d))
        end
    end
    println()

    n2, nd2, nc2, ncd2, leaves = rate_constant(nkeys, nsig)
    println("## B -- fpr_inv_sigma vs 1/sigma (perturbs the WIDTH)")
    @printf("  tree leaves differing  : %.1f of 512 per key\n", leaves)
    @printf("  signatures differing   : %d of %d   (rate %.3g)\n", nd2, n2, nd2 / n2)
    @printf("  coefficients differing : %d of %d\n", ncd2, nc2)
    nd2 == 0 && @printf("  zero events -> 95%% upper bound on the rate: %.3g\n", 3 / n2)
    println()
    println("# ePrint 2024/1709 reports ~1e-4 per signature for its own")
    println("# perturbation class (build configuration), by the same mechanism.")
end

main()
