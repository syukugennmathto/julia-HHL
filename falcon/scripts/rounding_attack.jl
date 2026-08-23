#!/usr/bin/env julia
#
# rounding_attack.jl -- on-demand key recovery from a brief perturbation of the
# floating-point ENVIRONMENT, with no source change and no waiting.
#
#     julia --project=falcon falcon/scripts/rounding_attack.jl [keys] [sigs]
#
# scripts/rounding_mode.jl shows that flipping the IEEE-754 rounding direction
# changes EVERY signature (1800 of 1800, against a nearest-vs-nearest control of
# 0 of 1800).  That is a determinism failure, but on its own it is not an attack:
# a whole-signature flip desynchronises the sampler, and the resulting difference
# is unstructured, which is exactly the case ePrint 2024/1709 section 5 says is
# useless for key recovery.
#
# The attack is to flip BRIEFLY.  The rounding mode is per-thread process state;
# an adversary who can set it for a short window -- a co-resident library, a
# signal handler, another thread in the same process, a hypervisor restoring an
# FPU context -- can arrange for only the TAIL of the tree traversal to be
# perturbed, and then restore it.  Then the two runs differ only in the last few
# sampled integers, Delta z0 is 2-sparse, and section 5.1's recovery applies:
#
#     ds0 = dz0 * g,  ds1 = -dz0 * f,  dz0 = a + b x^{n/2}
#
# which is a search over (a,b) in [-19,19]^2.  Unlike the natural events measured
# elsewhere in this project (7.5e-6 to 1.9e-5 per signature, so tens of thousands
# of signatures of waiting), this is available ON DEMAND on the next signature.
#
# Requires the same syndrome to be signed twice, i.e. the derandomized variants
# (IBE key extraction, SNARK-friendly signatures, aggregation) that 2024/1709
# identifies -- the adversary supplies the second signature by re-signing with
# the window active.

using Falcon
using Printf

const F = Falcon
const FE_TONEAREST = Cint(0)
const FE_UPWARD    = Cint(0x800)
const FE_DOWNWARD  = Cint(0x400)
setround(m) = ccall(:fesetround, Cint, (Cint,), m)

Falcon.eval(quote
    const CALLN    = Ref(0)
    const FLIP_AT  = Ref(0)     # flip the mode when the call counter reaches this
    const FLIP_END = Ref(0)     # restore it when the counter reaches this
    const FLIP_MODE = Ref(Cint(0))
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        CALLN[] += 1
        if FLIP_AT[] != 0
            CALLN[] == FLIP_AT[]  && ccall(:fesetround, Cint, (Cint,), FLIP_MODE[])
            CALLN[] == FLIP_END[] && ccall(:fesetround, Cint, (Cint,), Cint(0))
        end
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

mul_sparse(v, a, b) = begin
    n = length(v); h = n ÷ 2; out = a .* v
    for k in 0:(n-1)
        src = k - h
        out[k+1] -= b * (src >= 0 ? v[src+1] : -v[src+n+1])
    end
    out
end
function exact_div(v, d)
    d == 0 && return nothing
    out = similar(v)
    for i in eachindex(v)
        q_, rr = divrem(v[i], d); rr == 0 || return nothing; out[i] = q_
    end
    out
end
function recover(ds0, ds1, h, q)
    redq(v) = Int[mod(c, q) for c in v]
    for a in -19:19, b in -19:19
        (a == 0 && b == 0) && continue
        d = a*a + b*b
        g0 = exact_div(mul_sparse(ds0, a, b), d); g0 === nothing && continue
        f0 = exact_div(mul_sparse(ds1, a, b), d); f0 === nothing && continue
        for (gc, fc) in ((g0, -f0), (-g0, f0), (g0, f0), (-g0, -f0))
            (maximum(abs, gc) <= 256 && maximum(abs, fc) <= 256) || continue
            redq(gc) == redq(F.polymulq(redq(fc), h)) || continue
            return (a, b, gc, fc)
        end
    end
    nothing
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 25
    p = FALCON_512; n = p.n; q = p.q; L = 2n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = q)
    setround(FE_TONEAREST)
    r = chacha20(collect(UInt8, 0x00:0x37))
    keys = [falcon_keygen(n, k -> randombytes!(r, k))[1] for _ in 1:nkeys]

    println("# rounding_attack.jl -- a brief rounding-mode window, then key recovery")
    println("# window = [call k, call 2n]; outside it the mode is FE_TONEAREST")
    println()
    @printf("%-8s %8s %10s %10s %12s\n", "k", "signed", "diverged", "recovered", "rate")
    ks = length(ARGS) >= 3 ? [parse(Int,ARGS[3])] : [L-1, L-2, L-3, L-4, L-6]
    for k in ks
        tot = 0; div = 0; rec = 0
        for (ki, sk) in enumerate(keys)
            h = F.polydivq(Int[mod(c,q) for c in sk.g], Int[mod(c,q) for c in sk.f])
            tf = Int.(sk.f); tg = Int.(sk.g)
            for j in 1:nsig
                st = shake256(codeunits("ra/$k/$ki/$j"), 56)
                setround(FE_TONEAREST); F.FLIP_AT[] = 0; F.CALLN[] = 0
                r1 = chacha20(st)
                a1, a2 = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
                setround(FE_TONEAREST)
                F.FLIP_AT[] = k; F.FLIP_END[] = L; F.FLIP_MODE[] = FE_UPWARD; F.CALLN[] = 0
                r2 = chacha20(st)
                b1, b2 = F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
                setround(FE_TONEAREST); F.FLIP_AT[] = 0
                tot += 1
                ds0 = Int.(a1) .- Int.(b1); ds1 = Int.(a2) .- Int.(b2)
                (count(!=(0), ds0) + count(!=(0), ds1)) == 0 && continue
                div += 1
                res = recover(ds0, ds1, h, q)
                if res !== nothing
                    rot(v,i) = [m < i ? -v[m-i+n+1] : v[m-i+1] for m in 0:(n-1)]
                    _,_,gc,fc = res
                    ok = any(sg -> any(i -> gc == sg .* rot(tg,i) && fc == sg .* rot(tf,i), 0:(n-1)), (1,-1))
                    ok && (rec += 1)
                end
            end
        end
        @printf("%-8d %8d %10d %10d %11.3g\n", k, tot, div, rec, rec/tot)
    end
    setround(FE_TONEAREST)
    println()
    println("# 'recovered' = the (a,b) search returned exactly the stored key up to the")
    println("# NTRU lattice's rotation/negation symmetries.  Compare with the natural")
    println("# event rates measured elsewhere: 1.9e-5 (A1), 7.5e-6 (A2).")
end

main()
