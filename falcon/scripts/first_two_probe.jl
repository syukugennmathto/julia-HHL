#!/usr/bin/env julia
#
# first_two_probe.jl -- construct FIRST-two-call divergences on demand and look
# at what they actually give the attacker.
#
#     julia --project=falcon falcon/scripts/first_two_probe.jl [keys] [msgs]
#
# ePrint 2024/1709 section 5.1 dismisses first-two divergences in one paragraph,
# without an experiment: "the entire remainder of the ffSampling computation is
# affected, which yields a large, somewhat unstructured difference ds between
# the two signatures.  This vector ds is a short lattice vector, but it is not
# expected to be short enough to make key recovery feasible."
#
# That is a claim about STRUCTURE, and it has a testable alternative.  A
# divergence at call 0 changes z(0) by 1.  If the sampler happens to consume the
# SAME number of random bytes in both runs, the PRNG stays synchronised, and
# every later call sees the same random tape with only a deterministically
# shifted centre -- in which case most later z's could still agree and dz would
# be SPARSE, making ds short and structured, i.e. exploitable.  If instead the
# byte consumption differs, the PRNG desynchronises and everything after is
# unrelated, which is the paper's assumption.  Nobody has measured which happens.
#
# Finding a first-two divergence is cheap if you go looking: the centre of call 0
# is n/q with q = 12289, so it is an exact integer for about one message in
# 12289 (Heuristic 1), and the message is the attacker's to choose.  We scan
# messages, and whenever a first-two centre lands on an integer we run both
# spellings and report what the difference looks like.

using Falcon
using Printf
using Statistics

const F = Falcon

Falcon.eval(quote
    const MU = Float64[]
    const USED = Ref(0)
    const ON = Ref(false)
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); USED[] += 9
            b = Int(rb(1)[1]) & 1; USED[] += 1
            z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            z0 = basesampler(rb); USED[] += 9
            b = Int(rb(1)[1]) & 1; USED[] += 1
            z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

gap(x) = abs(x - round(x))

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 40
    nmsg  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1500
    p = FALCON_512; n = p.n; L = 2n; q = p.q
    r = chacha20(collect(UInt8, 0x00:0x37))

    println("# first_two_probe.jl  n=", n, "  scanning ", nkeys*nmsg, " (key,message) pairs")
    println("# looking for an EXACT integer centre at sampler call 1 or 2 (density ~2/q)")
    println()

    hits = 0; diverged = 0
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        for j in 1:nmsg
            pt = hash_to_point(collect(codeunits("ft/$ki/$j")),
                               shake256(codeunits("fts/$ki/$j"), SALT_LEN), n; q = q)
            rb = chacha20(shake256(codeunits("ftp/$ki/$j"), 56))
            empty!(F.MU); F.ON[] = true
            s1a, s2a = F.sample_preimage(sk, pt, x -> randombytes!(rb, x))
            F.ON[] = false
            length(F.MU) == L || continue
            (gap(F.MU[1]) < 1e-9 || gap(F.MU[2]) < 1e-9) || continue

            hits += 1
            which = gap(F.MU[1]) < 1e-9 ? 1 : 2
            muhit = F.MU[which]

            # same key, same message, same PRNG state -- the two spellings
            rb1 = chacha20(shake256(codeunits("ftp/$ki/$j"), 56))
            rb2 = chacha20(shake256(codeunits("ftp/$ki/$j"), 56))
            F.USED[] = 0; empty!(F.MU); F.ON[] = true
            b1, b2 = F.sample_preimage(sk, pt, x -> randombytes!(rb1, x))
            F.ON[] = false; u1 = F.USED[]; MUA = copy(F.MU)
            F.USED[] = 0; empty!(F.MU); F.ON[] = true
            c1, c2 = with_spec_ffsampling() do
                F.sample_preimage(sk, pt, x -> randombytes!(rb2, x))
            end
            F.ON[] = false; u2 = F.USED[]; MUB = copy(F.MU)
            dmu = length(MUB) >= which ? abs(MUA[which] - MUB[which]) : NaN
            @printf("   centres at that call: C %.17g  spec %.17g   |dmu| = %.3g%s\n",
                    MUA[which], length(MUB) >= which ? MUB[which] : NaN, dmu,
                    dmu == 0.0 ? "   <- EXACTLY equal: no straddle possible" : "")

            ds1 = Int.(b1) .- Int.(c1); ds2 = Int.(b2) .- Int.(c2)
            nd = count(!=(0), ds1) + count(!=(0), ds2)
            @printf("hit %d: key %d msg %d, integer centre at call %d (mu = %.1f)\n",
                    hits, ki, j, which, muhit)
            if nd == 0
                println("   -> the two spellings AGREE here (perturbation did not straddle)")
                continue
            end
            diverged += 1
            dsn = sqrt(sum(x -> Float64(x)^2, ds1) + sum(x -> Float64(x)^2, ds2))
            sgn = sqrt(sum(x -> Float64(x)^2, Int.(b1)) + sum(x -> Float64(x)^2, Int.(b2)))
            @printf("   -> DIVERGED. %d of %d signature coefficients differ\n", nd, 2n)
            @printf("      PRNG bytes consumed: %d vs %d  -> %s\n", u1, u2,
                    u1 == u2 ? "SYNCHRONISED (structured difference possible)" :
                               "DESYNCHRONISED (difference is unrelated)")
            @printf("      ||ds|| = %.1f   ||s|| = %.1f   ratio = %.3f\n", dsn, sgn, dsn/sgn)
            @printf("      beta bound sqrt = %.1f ; a key-recovery difference would be O(1)\n",
                    sqrt(Float64(p.sig_bound)))
        end
    end
    println()
    @printf("# %d integer-centre hits at calls 1-2 in %d pairs (rate %.3g); %d diverged\n",
            hits, nkeys*nmsg, hits/(nkeys*nmsg), diverged)
    println("# Expected hit rate from Heuristic 1: 2/q = ", @sprintf("%.3g", 2/q))
end

main()
