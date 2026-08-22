#!/usr/bin/env julia
#
# linear_form.jl -- the exact centre of the first two sampler calls is a PUBLIC
# LINEAR FORM in the secret f.
#
#     julia --project=falcon falcon/scripts/linear_form.jl
#
# In `sample_preimage`, t1 = (-point_fft * b)/q with b = B0[1,2] = -f, so as a
# ring element t1 = (c*f)/q exactly, where c is the public hashed message.  The
# first descent samples z1 from t1, so the first two sampler centres are
# coefficients of (c*f)/q.  Verified here against exact Int128 negacyclic
# arithmetic: over every key and message tried, the centres match
#
#     call 1 : mu = (c*f)[n/2 - 1] / q          (index 256, 1-based)
#     call 2 : mu = (c*f)[n   - 1] / q          (index 512, 1-based)
#
# to within 4.6e-14, i.e. floating-point error alone, at FIXED indices.
#
# Two consequences.
#
# (1) Heuristic 1 becomes a THEOREM at the first two positions.  ePrint
#     2024/1709 Remark 1 explains why its Heuristic 1 cannot be made a theorem
#     in general; at these two positions it can, because the denominator is
#     exactly q and the numerator is an explicit integer.  "Integer centre" is
#     exactly the condition (c*f)[n/2-1] = 0 (mod q).
#
# (2) That condition is LINEAR IN THE SECRET, with public coefficients.  Each
#     detected integer-centre event gives one equation <c, x^{n/2} f> = 0 mod q;
#     n independent events determine f mod q, and since ||f|| << q that
#     determines f.  This is a different attack shape from ePrint 2024/1709
#     section 5.1, which reads the secret out of the difference VECTOR of two
#     signatures: here the leak is the EVENT, and the difference vector is not
#     needed at all.  It applies at the first two calls -- the positions that
#     paper dismisses as yielding only a short lattice vector.
#
#     What it needs is an oracle for the event.  The respelling perturbation
#     does NOT provide one: scripts/first_two_probe.jl manufactures integer
#     centres at these positions and none of them straddles, because |dmu| there
#     sits below the shared rounding offset (docs/debug_log.md #062).  A larger
#     perturbation source (FMA), a conformance-mismatch report, or a side
#     channel would.  We state the shape and its requirement; we do not claim a
#     working attack on plain Falcon.

using Falcon, Printf
const F = Falcon
# In sample_preimage:  t1 = (-point_fft * b)/q  with b = B0[1,2] = -f,
# so t1 = (c * f)/q  exactly, as a ring element.  The first sampler calls come
# from the t1 branch.  If so, the exact centre at calls 1,2 is a coefficient of
# (c*f)/q -- a PUBLIC LINEAR FORM in the secret f -- and "integer centre" is the
# linear condition (c*f)[j] = 0 mod q.  Test it in exact integer arithmetic.
Falcon.eval(quote
    const MU = Float64[]; const ON = Ref(false)
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
"negacyclic product in Z[x]/(x^n+1), exact"
function negamul(a::Vector{Int}, b::Vector{Int})
    n = length(a); w = zeros(Int128, n)
    for i in 0:n-1, j in 0:n-1
        k = i+j; c = Int128(a[i+1])*b[j+1]
        w[mod(k,n)+1] += k < n ? c : -c
    end
    w
end
function main(nkeys, nmsg)
    p = FALCON_512; n = p.n; q = p.q
    r = chacha20(collect(UInt8, 0x00:0x37))
    best = fill(Inf, 4); idx = zeros(Int, 4)
    checked = 0
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        f = Int.(sk.f)
        for j in 1:nmsg
            pt = hash_to_point(collect(codeunits("lin/$ki/$j")),
                               shake256(codeunits("lins/$ki/$j"), SALT_LEN), n; q=q)
            rb = chacha20(shake256(codeunits("linp/$ki/$j"), 56))
            empty!(F.MU); F.ON[]=true
            F.sample_preimage(sk, pt, x -> randombytes!(rb, x)); F.ON[]=false
            length(F.MU) == 2n || continue
            cf = negamul(Int.(pt), f)          # exact  c*f  in Z[x]/(x^n+1)
            checked += 1
            # find which coefficient of (c*f)/q each of the first two centres is
            for (slot, call) in enumerate((1,2))
                mu = F.MU[call]
                # search all coefficients for the matching one
                bestd = Inf; bi = 0
                for k in 1:n
                    d = abs(mu - Float64(cf[k])/q)
                    d < bestd && (bestd = d; bi = k)
                end
                if bestd < best[slot]; best[slot] = bestd; idx[slot] = bi; end
                if slot <= 2 && checked <= 3
                    @printf("  key %d msg %d  call %d: mu = %.10f   (c*f)[%d]/q = %.10f   |diff| = %.3g\n",
                            ki, j, call, mu, bi, Float64(cf[bi])/q, bestd)
                end
            end
        end
    end
    println()
    @printf("over %d signatures: best match |mu - (c*f)[k]/q|  call 1: %.3g (k=%d)   call 2: %.3g (k=%d)\n",
            checked, best[1], idx[1], best[2], idx[2])
end
main(3, 4)
