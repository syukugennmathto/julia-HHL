#!/usr/bin/env julia
#
# precision_law.jl -- how the signature-divergence rate depends on the size of
# the floating-point perturbation on the sampler centre.
#
#     julia --project=falcon falcon/scripts/precision_law.jl [keys] [sigs]
#
# ePrint 2024/1709 Lemma 1 predicts that, for a perturbation eps on the centre,
# the probability of an inconsistent execution is linear in eps (until eps ~ 1).
# We verify that law directly and calibrate its constant, by INJECTING a
# controlled perturbation of magnitude eps at every sampler call and measuring
# the resulting divergence rate, for a range of eps.  This connects to the
# specification's own precision analysis (sec 2.5.2), which bounds the effect
# of finite precision on the sampled *distribution* (delta_c + delta_sigma <=
# 2^-46) -- a different quantity from the reproducibility rate measured here.
#
# The naturally occurring respelling perturbations sit at |dmu| ~ 1e-14
# (scripts/position_profile.jl); this sweep places them on the rate(eps) curve
# and lets one read off the precision at which the rate reaches 1 (a guaranteed
# leak on every signature).

using Falcon
using Printf

const F = Falcon

# a copy of samplerz_isigma that adds EPS[] to the centre when ON[]
Falcon.eval(quote
    const EPS = Ref(0.0)
    const ON = Ref(false)
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        m = ON[] ? mu + EPS[] : mu
        s = Int(floor(m)); r = m - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

function rate_at(sk, pt, nsig, eps, tag)
    ndiff = 0
    for j in 1:nsig
        st = shake256(codeunits("$tag/$j"), 56)
        r1 = chacha20(st); r2 = chacha20(st)
        F.EPS[] = 0.0; F.ON[] = false
        _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
        F.EPS[] = eps; F.ON[] = true
        _, b = F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
        F.ON[] = false
        Int.(a) != Int.(b) && (ndiff += 1)
    end
    return ndiff
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
    p = FALCON_512; n = p.n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)
    r = chacha20(collect(UInt8, 0x00:0x37))
    keys = [falcon_keygen(n, k -> randombytes!(r, k))[1] for _ in 1:nkeys]

    println("# precision_law.jl  n=", n, "  ", nkeys, " keys x ", nsig,
            " sigs = ", nkeys*nsig, " per eps")
    println("# eps = 2^-t injected at every sampler centre; rate = P(signature differs)")
    @printf("%4s  %-11s  %10s  %-11s\n", "t", "eps", "events", "rate")
    for t in (8, 10, 12, 14, 16, 18, 20, 22)
        eps = 2.0^(-t)
        tot = 0
        for (ki, sk) in enumerate(keys)
            tot += rate_at(sk, pt, nsig, eps, "pl/$t/$ki")
        end
        @printf("%4d  %-11.3e  %10d  %-11.3e\n", t, eps, tot, tot/(nkeys*nsig))
    end
    println()
    println("# For reference, the naturally-occurring respelling perturbations")
    println("# (position_profile.jl) sit at |dmu| ~ 1e-14, i.e. t ~ 46.")
end

main()
