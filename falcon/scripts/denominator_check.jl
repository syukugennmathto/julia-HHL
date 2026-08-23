#!/usr/bin/env julia
#
# denominator_check.jl -- verify, directly, the arithmetic structure that both
# the attack of ePrint 2024/1709 and its section 7.1 countermeasure rest on.
#
#     julia --project=falcon falcon/scripts/denominator_check.jl [keys] [sigs]
#
# Theorem 1 / Heuristic 1 of that paper say the exact theoretical centre of
# sampler call 2k (and 2k+1) is a rational n_k / g_k with
#     g_0 = q                       (first two calls)
#     g_{n-1} = m_1 = ||(g,-f)||^2  (last two calls)
# Nobody has checked this against a running implementation.  We do, by testing
# whether g * c is (very near) an integer for the predicted g.
#
# The consequence is what makes this worth measuring.  A centre n/g is
#   an INTEGER      iff g | n                     -- possible for any g
#   a HALF-INTEGER  iff 2n = g(2j+1)              -- possible ONLY if g is EVEN
#                                                    (g odd forces g | n, an integer)
# Falcon's SamplerZ splits at integers (floor), so it is sensitive at integer
# centres.  The countermeasure's NewSamplerZ splits at half-integers (round), so
# it is sensitive at half-integer centres.  Therefore:
#
#   * first two calls:  g_0 = q = 12289, ODD  -> no half-integer centres, EVER.
#     Rounding immunises these positions completely.
#   * last two calls:   g_{n-1} = ||(g,-f)||^2, which the reference key
#     generator always makes EVEN -> half-integer centres occur at the same
#     density 1/g as integer centres did.  Rounding immunises nothing here.
#
# And the last two calls are exactly the positions from which section 5 of that
# paper recovers the whole private key; the first two give only a short lattice
# vector it deems useless.  So deploying part 1 of the countermeasure without
# part 2 removes the harmless half of the exposure and leaves the dangerous half
# untouched -- which is what scripts/countermeasure_eval.jl measures directly.

using Falcon
using Printf

const F = Falcon

Falcon.eval(quote
    const MU = Float64[]
    const ON = Ref(false)
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

gap(x) = abs(x - round(x))
halfgap(x) = abs(x - 0.5 - round(x - 0.5))   # distance to nearest half-integer

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 25
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 60
    p = FALCON_512; n = p.n; L = 2n; q = p.q
    r = chacha20(collect(UInt8, 0x00:0x37))

    println("# denominator_check.jl  n=", n, "  ", nkeys, " keys x ", nsig, " messages")
    println()
    println("# Is g*c an integer for the predicted denominator g?")
    println("# (max residual over all samples; ~1e-10 or below means yes)")

    maxres_first = 0.0; maxres_last = 0.0
    nfirst = 0; nlast = 0
    halfnear_first = 0; halfnear_last = 0
    intnear_first = 0; intnear_last = 0
    tvals = Int[]

    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        t = Int(sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f))
        push!(tvals, t)
        for j in 1:nsig
            pt = hash_to_point(collect(codeunits("dc/$ki/$j")),
                               shake256(codeunits("dcsalt/$ki/$j"), SALT_LEN), n; q = q)
            rb = chacha20(shake256(codeunits("dcp/$ki/$j"), 56))
            empty!(F.MU); F.ON[] = true
            F.sample_preimage(sk, pt, x -> randombytes!(rb, x))
            F.ON[] = false
            length(F.MU) == L || continue
            for i in (1, 2)
                c = F.MU[i]; nfirst += 1
                maxres_first = max(maxres_first, gap(q * c))
                halfgap(c) < 1e-9 && (halfnear_first += 1)
                gap(c)     < 1e-9 && (intnear_first  += 1)
            end
            for i in (L-1, L)
                c = F.MU[i]; nlast += 1
                maxres_last = max(maxres_last, gap(t * c))
                halfgap(c) < 1e-9 && (halfnear_last += 1)
                gap(c)     < 1e-9 && (intnear_last  += 1)
            end
        end
    end

    @printf("  calls 1,2       : g = q = %d              max |g*c - round(g*c)| = %.3e   (%d samples)\n",
            q, maxres_first, nfirst)
    @printf("  calls 2n-1,2n   : g = ||(g,-f)||^2 (per key)  max |g*c - round(g*c)| = %.3e   (%d samples)\n",
            maxres_last, nlast)
    println()
    @printf("  ||(g,-f)||^2 over %d keys: all even? %s   (odd: %d)\n",
            nkeys, all(iseven, tvals) ? "YES" : "no", count(isodd, tvals))
    println()
    println("# Centres landing on the two kinds of discontinuity:")
    @printf("%-18s %12s %12s\n", "position", "integer", "half-integer")
    @printf("%-18s %12d %12d      <- floor(SamplerZ) / round(NewSamplerZ)\n",
            "calls 1,2", intnear_first, halfnear_first)
    @printf("%-18s %12d %12d\n", "calls 2n-1,2n", intnear_last, halfnear_last)
    println()
    println("# Arithmetic consequence (independent of how many samples we drew):")
    println("#   q = 12289 is ODD  -> n/q is never a half-integer -> rounding")
    println("#     immunises calls 1,2 COMPLETELY.")
    println("#   ||(g,-f)||^2 is EVEN on every reference key -> n/t can be a")
    println("#     half-integer, at the same density 1/t as integers -> rounding")
    println("#     immunises calls 2n-1,2n NOT AT ALL.")
    println("#   Calls 2n-1,2n are the full-key-recovery positions (2024/1709 sec 5).")
end

main()
