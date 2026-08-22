#!/usr/bin/env julia
#
# countermeasure_break.jl -- with ePrint 2024/1709's own countermeasure
# (Algorithm 4, "NewSamplerZ") deployed exactly as an implementer can deploy it
# today, recover the private key anyway.
#
#     julia --project=falcon falcon/scripts/countermeasure_break.jl [keys] [sigs]
#
# "As an implementer can deploy it today" means: part 1 of section 7.1 (round
# instead of floor) but NOT part 2 (odd ||(g,-f)||^2), because the reference key
# generator cannot produce an odd-norm key -- it forces the coefficient sums of
# both f and g odd, so their sum is even (0 odd keys in 3200 generated here).
#
# Why that still breaks: a centre n/g is a half-integer iff 2n = g(2j+1), which
# needs g EVEN.  At the last two sampler calls g = ||(g,-f)||^2, always even on
# reference keys, so half-integer centres -- NewSamplerZ's discontinuity -- occur
# at the same density 1/g that integer centres did for the original sampler.  And
# the last two calls are exactly the positions from which section 5 of that paper
# recovers the whole key.
#
# The scan is targeted rather than blind: we sign once, look at the last two
# centres, and only pay for the comparison run when one of them is within 1e-12
# of a half-integer.  That is roughly a 2x saving over comparing every pair.

using Falcon
using Printf

const F = Falcon
setprecision(BigFloat, 256)

function build_new_rcdt()
    smax = BigFloat("1.8205"); den = 2 * smax * smax
    w = [i == 0 ? BigFloat("0.5") : exp(-(BigFloat(i)^2 - BigFloat(i)) / den) for i in 0:18]
    p = w ./ sum(w)
    UInt128[UInt128(round(BigInt, sum(p[(k+1):end]) * BigFloat(2)^72)) for k in 1:18]
end

Falcon.eval(quote
    const NEW_RCDT = $(build_new_rcdt())
    const MU = Float64[]
    const ON = Ref(false)
    function newbasesampler(rb)
        bytes = rb(RCDT_PREC >> 3); u = UInt128(0)
        for i in 1:(RCDT_PREC >> 3); u |= UInt128(bytes[i]) << (8*(i-1)); end
        y = 0; for elt in NEW_RCDT; y += Int(u < elt); end
        return y
    end
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = round(mu); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            yp = newbasesampler(rb); b = Int(rb(1)[1]) & 1; y = (2b-1)*yp
            x = ((y-r)^2)*dss - (yp*yp - yp)*INV_2SIGMA2
            berexp(x, ccs, rb) && return y + Int(s)
        end
    end
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = round(mu); r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            yp = newbasesampler(rb); b = Int(rb(1)[1]) & 1; y = (2b-1)*yp
            x = ((y-r)^2)*dss - (yp*yp - yp)*INV_2SIGMA2
            berexp(x, ccs, rb) && return y + Int(s)
        end
    end
end)

halfgap(x) = abs(x - 0.5 - round(x - 0.5))

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

"Recover (f,g) from one last-two divergence, per 2024/1709 section 5.1."
function recover(ds0, ds1, h, q, truef, trueg)
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
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1500
    p = FALCON_512; n = p.n; q = p.q; L = 2n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = q)
    r = chacha20(collect(UInt8, 0x00:0x37))

    println("# countermeasure_break.jl -- NewSamplerZ (2024/1709 Alg 4) deployed,")
    println("# on reference keys (||(g,-f)||^2 even, part 2 unavailable).")
    println("# scanning ", nkeys*nsig, " signatures for a half-integer straddle")
    println()
    flagged = 0; diverged = 0; recovered = 0
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        h = F.polydivq(Int[mod(c,q) for c in sk.g], Int[mod(c,q) for c in sk.f])
        for j in 1:nsig
            st = shake256(codeunits("cb/$ki/$j"), 56)
            rb = chacha20(st)
            empty!(F.MU); F.ON[] = true
            F.sample_preimage(sk, pt, x -> randombytes!(rb, x))
            F.ON[] = false
            length(F.MU) == L || continue
            (halfgap(F.MU[L-1]) < 1e-12 || halfgap(F.MU[L]) < 1e-12) || continue
            flagged += 1

            r1 = chacha20(st); r2 = chacha20(st)
            a1, a2 = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
            b1, b2 = with_spec_ffsampling() do
                F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
            end
            ds0 = Int.(a1) .- Int.(b1); ds1 = Int.(a2) .- Int.(b2)
            nd = count(!=(0), ds0) + count(!=(0), ds1)
            nd == 0 && continue
            diverged += 1
            @printf("DIVERGENCE key %d sig %d : %d of %d coefficients differ\n", ki, j, nd, 2n)
            res = recover(ds0, ds1, h, q, Int.(sk.f), Int.(sk.g))
            if res === nothing
                println("   key recovery: FAILED on this pair")
            else
                a, b, gc, fc = res
                # the NTRU lattice's shortest vectors are the rotations
                # +/- x^i (f,g); any of them is an equally usable signing key.
                rot(v, i) = [k < i ? -v[k - i + n + 1] : v[k - i + 1] for k in 0:(n-1)]
                tf = Int.(sk.f); tg = Int.(sk.g)
                which = ""
                for i in 0:(n-1), sgn in (1, -1)
                    if gc == sgn .* rot(tg, i) && fc == sgn .* rot(tf, i)
                        which = sgn == 1 ? "= x^$i (f,g)" : "= -x^$i (f,g)"
                        break
                    end
                end
                nrm(a_, b_) = sum(x -> Int128(x)^2, a_) + sum(x -> Int128(x)^2, b_)
                ok = which != ""
                recovered += ok
                @printf("   *** KEY RECOVERED *** (a,b)=(%d,%d)  reproduces h; %s\n",
                        a, b, ok ? "exactly the stored key up to lattice symmetry " * which :
                             "short (f,g), norm " * string(nrm(fc,gc)) * " vs true " * string(nrm(tf,tg)))
            end
        end
    end
    println()
    @printf("# %d near-half-integer centres flagged, %d diverged, %d keys recovered, in %d signatures\n",
            flagged, diverged, recovered, nkeys*nsig)
    println("# A recovery here means the deployable half of the section 7.1")
    println("# countermeasure does not prevent full key recovery.")
end

main()
