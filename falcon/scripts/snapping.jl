#!/usr/bin/env julia
#
# snapping.jl -- a countermeasure that needs no key-generation change.
#
#     julia --project=falcon falcon/scripts/snapping.jl
#
# ePrint 2024/1709 section 7.1 fixes the sensitivity in two parts: round instead
# of floor, AND force ||(g,-f)||^2 odd.  Part 2 is impossible on the reference
# key generator, and part 1 alone leaves every key-recovery position exposed --
# we recover keys straight through it (scripts/countermeasure_break.jl).
#
# There is a third option that needs neither part.  At the sensitive positions
# the exact centre is a rational n/g whose denominator g the SIGNER already
# knows: g = q at the first two calls, g = ||(g,-f)||^2 at the last two.  The
# floating-point centre mu_hat approximates n/g to about 1e-13, so
#
#     g * mu_hat  is within  g * 1e-13 ~ 1.6e-9  of the integer n
#
# -- seven orders of magnitude of margin -- and therefore
#
#     n = round(g * mu_hat)          recovers the numerator EXACTLY,
#     s = fld(n, g)                  is then exact integer arithmetic.
#
# Measured margin (scripts/denominator_check.jl, 100000 samples per end):
# max |g*mu - round(g*mu)| = 1.4e-9 at the first two calls and 2.4e-8 at the
# last two, against the 0.5 that would be needed to break it.
#
# The result is that the split of the centre stops depending on floating point
# at exactly the places where that dependence is dangerous, WITHOUT touching key
# generation, without changing the sampled distribution (the fractional part
# r = mu - s is still computed in floating point, and Lemma 2 of that paper says
# the sampler is insensitive to perturbations that do not cross the split), and
# without needing the odd-norm keys the reference cannot produce.
#
# This script demonstrates it on the divergent pairs already found: every one of
# them stops diverging.

using Falcon
using Printf

const F = Falcon

Falcon.eval(quote
    const CALLN   = Ref(0)      # which sampler call we are at, 1-based
    const NCALLS  = Ref(1024)   # 2n
    const GLAST   = Ref(0)      # ||(g,-f)||^2 for the current key
    const GFIRST  = Ref(12289)  # q
    const SNAP    = Ref(false)

    "The denominator of the exact centre at this position, or 0 if unknown."
    function snap_g(i)
        (i == 1 || i == 2) && return GFIRST[]
        (i >= NCALLS[] - 1) && return GLAST[]
        return 0
    end

    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        CALLN[] += 1
        s = if SNAP[]
                g = snap_g(CALLN[])
                g == 0 ? Int(floor(mu)) : Int(fld(round(Int128, big(g) * mu), big(g)))
            else
                Int(floor(mu))
            end
        r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        CALLN[] += 1
        s = if SNAP[]
                g = snap_g(CALLN[])
                g == 0 ? Int(floor(mu)) : Int(fld(round(Int128, big(g) * mu), big(g)))
            else
                Int(floor(mu))
            end
        r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

respelled_tree(sk) = begin
    t = with_spec_spelling(cdiv = true, ldl = true) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

"Replay one known-divergent (arm, key, sig) with snapping off and on."
function replay(arm, keyidx, sigidx)
    p = FALCON_512; n = p.n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)
    r = chacha20(collect(UInt8, 0x00:0x37))
    local sk
    for _ in 1:keyidx; sk = falcon_keygen(n, k -> randombytes!(r, k))[1]; end
    F.GLAST[] = Int(sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f))
    F.NCALLS[] = 2n
    sk2 = arm == "A1" ? sk : respelled_tree(sk)
    sign2 = arm == "A1" ? with_spec_ffsampling : (f -> f())
    st = shake256(codeunits("$arm/$keyidx/$sigidx"), 56)
    out = String[]
    for snap in (false, true)
        F.SNAP[] = snap
        r1 = chacha20(st); r2 = chacha20(st)
        F.CALLN[] = 0; _, a = F.sample_preimage(sk,  pt, x -> randombytes!(r1, x))
        F.CALLN[] = 0; _, b = sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
        d = count(Int.(a) .!= Int.(b))
        push!(out, @sprintf("snapping %-3s : %4d of %d coefficients differ  %s",
                            snap ? "ON" : "off", d, n, d == 0 ? "AGREE" : "DIVERGE"))
    end
    F.SNAP[] = false
    @printf("%s key %d sig %d\n", arm, keyidx, sigidx)
    for l in out; println("   ", l); end
    return out
end

"Snapping must not change signatures that were never near a boundary."
function sanity(nkeys, nsig)
    p = FALCON_512; n = p.n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)
    r = chacha20(collect(UInt8, 0x00:0x37))
    changed = 0; total = 0
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        F.GLAST[] = Int(sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f))
        F.NCALLS[] = 2n
        for j in 1:nsig
            st = shake256(codeunits("snap/$ki/$j"), 56)
            F.SNAP[] = false; r1 = chacha20(st); F.CALLN[] = 0
            s1a, s2a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
            F.SNAP[] = true;  r2 = chacha20(st); F.CALLN[] = 0
            s1b, s2b = F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
            total += 1
            (Int.(s2a) != Int.(s2b)) && (changed += 1)
            F.SNAP[] = false
        end
    end
    @printf("\nsanity: snapping changed %d of %d signatures (same spelling, same PRNG)\n",
            changed, total)
    println("   (should be 0: snapping only makes the split exact, it does not move it)")
end

println("# snapping.jl -- replaying every divergent pair this project has found")
println("# with the exact-rational split turned off, then on")
println()
for (arm, k, s) in (("A2", 70, 534), ("A2", 65, 1239), ("A2", 68, 1885), ("A1", 4, 105))
    replay(arm, k, s)
end
sanity(20, 60)
