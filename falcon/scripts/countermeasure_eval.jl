#!/usr/bin/env julia
#
# countermeasure_eval.jl -- an independent implementation and evaluation of the
# countermeasure proposed in ePrint 2024/1709 section 7.1.
#
#     julia --project=falcon falcon/scripts/countermeasure_eval.jl <mode> [args]
#
#       chisq                      -- validate NewSamplerZ's output distribution
#       rate <arm> <keys> <sigs> <chunk> <out>   -- divergence rate with NewSamplerZ
#       parity <keys>              -- check the m1/m2/m3 parity argument on real keys
#
# ===========================================================================
# WHAT THE COUNTERMEASURE IS, AND WHY IT HAS TWO PARTS
# ===========================================================================
#
# The sensitivity exists because SamplerZ splits its centre as
# `s = floor(mu); r = mu - s`, and `floor` is discontinuous at the integers,
# where Falcon's centres land with probability 1/q (first two calls) or
# 1/||(g,-f)||^2 (last two).  Algorithm 4 of that paper ("NewSamplerZ") rounds
# to NEAREST instead:
#
#     r  <- c - round(c)                   in [-1/2, 1/2)
#     y+ <- NewBaseSampler()               one-sided, centred at 1/2
#     b  <- {0,1};  y <- (2b-1) y+         note: NOT b + (2b-1)y+
#     x  <- (y-r)^2/(2 sigma^2) - (y+^2 - y+)/(2 sigma_max^2)
#     return y + round(c) w.p. (sigma_min/sigma) exp(-x)
#
# That moves the discontinuity from the integers to the HALF-integers.  On its
# own this does nothing, because a centre n_k/g_k can be a half-integer just as
# easily as an integer.  So the countermeasure has a SECOND part: force the six
# in-precision denominators g_0, g_1, g_2, g_{n-3}, g_{n-2}, g_{n-1} to be odd,
# which (via m_1 = t, m_2 = t^2 - 2u^2, m_3 = t^3 - 2t(u^2+v^2+w^2) + 2u(v-w)^2)
# reduces to requiring t = ||(g,-f)||^2 ODD.  With g_k odd, n_k/g_k = j + 1/2
# would need 2 n_k = g_k (2j+1) with g_k odd, forcing g_k | n_k and hence an
# integer -- so half-integer centres cannot occur at all.
#
# THE CATCH, which the paper states in one sentence and does not measure: the C
# reference implementation only ever generates keys with ||(g,-f)||^2 EVEN (it
# forces the coefficient sums of BOTH f and g odd, so their sum is even).  So on
# the reference key generator, part 2 is unavailable and part 1 alone is what an
# implementer would deploy.  This script measures what that costs.
#
# The prediction, before running: with an EVEN t = ||(g,-f)||^2, write t = 2m.
# A centre n/t is a half-integer iff 2n = t(2j+1) iff n = m(2j+1), i.e. exactly
# one residue class of n modulo t -- the same density 1/t as the integer-centre
# condition it replaced.  So the rate should be UNCHANGED, not reduced.
#
# CONSTANT TIME: out of scope, as everywhere in this project.

using Falcon
using Printf
using Statistics

const F = Falcon
setprecision(BigFloat, 256)

# ---------------------------------------------------------------------------
# NewBaseSampler's reverse CDT
# ---------------------------------------------------------------------------
# The paper's D: proportional to rho_{sigma_max, 1/2}(i) for 1 <= i <= 18, and
# to (1/2) rho_{sigma_max, 1/2}(0) for i = 0.  Since
# rho_{s,1/2}(i) = exp(-(i-1/2)^2/(2 s^2)) and the constant exp(-1/(8 s^2))
# cancels in the normalisation, the weights are w(i) = exp(-(i^2-i)/(2 s^2))
# for i >= 1 and w(0) = 1/2.  Halving i = 0 is what makes y = (2b-1) y+ come out
# symmetric: y = 0 is reachable from both values of b.
#
# The table is built at 256 bits and stored the way Falcon stores RCDT: entry k
# is 2^72 * P(X >= k), and the sampler counts how many entries a 72-bit uniform
# falls below.
function build_new_rcdt()
    smax = BigFloat("1.8205")
    den = 2 * smax * smax
    w = [i == 0 ? BigFloat("0.5") : exp(-(BigFloat(i)^2 - BigFloat(i)) / den)
         for i in 0:18]
    tot = sum(w)
    p = w ./ tot
    tail = [sum(p[(k + 1):end]) for k in 1:18]      # P(X >= k), k = 1..18
    scale = BigFloat(2)^72
    return UInt128[UInt128(round(BigInt, t * scale)) for t in tail], p
end

const RCDT_NEW, PMF_NEW = build_new_rcdt()

# ---------------------------------------------------------------------------
# NewSamplerZ, installed over the library's samplers
# ---------------------------------------------------------------------------
Falcon.eval(quote
    const NEW_RCDT = $(RCDT_NEW)
    const USE_NEW = Ref(true)

    function newbasesampler(randombytes)
        bytes = randombytes(RCDT_PREC >> 3)
        u = UInt128(0)
        for i in 1:(RCDT_PREC >> 3)
            u |= UInt128(bytes[i]) << (8 * (i - 1))
        end
        y = 0
        for elt in NEW_RCDT
            y += Int(u < elt)
        end
        return y
    end

    "Algorithm 4 of ePrint 2024/1709, in the sigma form."
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        USE_NEW[] || return _samplerz_floor(mu, sigma, sigmin, rb)
        s = round(mu); r = mu - s
        dss = 1 / (2 * sigma * sigma); ccs = sigmin / sigma
        while true
            yp = newbasesampler(rb)
            b = Int(rb(1)[1]) & 1
            y = (2b - 1) * yp
            x = ((y - r)^2) * dss - (yp * yp - yp) * INV_2SIGMA2
            berexp(x, ccs, rb) && return y + Int(s)
        end
    end

    "Algorithm 4, in the reciprocal-sigma form the C reference's tree uses."
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        USE_NEW[] || return _samplerz_isigma_floor(mu, isigma, sigmin, rb)
        s = round(mu); r = mu - s
        dss = 0.5 * (isigma * isigma); ccs = isigma * sigmin
        while true
            yp = newbasesampler(rb)
            b = Int(rb(1)[1]) & 1
            y = (2b - 1) * yp
            x = ((y - r)^2) * dss - (yp * yp - yp) * INV_2SIGMA2
            berexp(x, ccs, rb) && return y + Int(s)
        end
    end

    # the originals, kept reachable for controls
    function _samplerz_floor(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        s = Int(floor(mu)); r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
    function _samplerz_isigma_floor(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

# ---------------------------------------------------------------------------
# mode: chisq -- does NewSamplerZ actually sample the right distribution?
# ---------------------------------------------------------------------------
function mode_chisq()
    println("# chi-square: NewSamplerZ (2024/1709 Alg 4) against D_{Z,sigma,mu}")
    println("# an independent check that the countermeasure preserves correctness")
    println()
    r = chacha20(shake256(codeunits("cm/chisq"), 56))
    rb = k -> randombytes!(r, k)
    N = 200_000
    for (mu, sigma) in ((0.0, 1.5), (0.3, 1.5), (-0.4, 1.7), (7.25, 1.29), (0.5, 1.8))
        cnt = Dict{Int,Int}()
        for _ in 1:N
            z = F.samplerz(mu, sigma, 1.2778336969128337, rb)
            cnt[z] = get(cnt, z, 0) + 1
        end
        lo = minimum(keys(cnt)); hi = maximum(keys(cnt))
        wts = [exp(-(BigFloat(k) - BigFloat(mu))^2 / (2 * BigFloat(sigma)^2)) for k in lo:hi]
        pk = wts ./ sum(wts)
        chi = 0.0; df = 0
        for (i, k) in enumerate(lo:hi)
            e = Float64(pk[i]) * N
            e < 20 && continue
            o = get(cnt, k, 0)
            chi += (o - e)^2 / e; df += 1
        end
        df -= 1
        @printf("  mu=%6.2f sigma=%.3f : chi2 = %8.2f, df = %2d, chi2/df = %.3f  %s\n",
                mu, sigma, chi, df, chi / df,
                chi / df < 2.0 ? "OK" : "*** SUSPECT ***")
    end
    println()
    println("# chi2/df near 1 means the rounding-based sampler is distributionally sound.")
end

# ---------------------------------------------------------------------------
# mode: parity -- the m1/m2/m3 argument, on real keys
# ---------------------------------------------------------------------------
function mode_parity(nkeys)
    println("# the parity argument of 2024/1709 section 7.1, on real keys")
    println("# m1 = t, m2 = t^2 - 2u^2, m3 = t^3 - 2t(u^2+v^2+w^2) + 2u(v-w)^2")
    println("# with t = <b0,b0> = ||(g,-f)||^2, u = <b0,b2>, v = <b0,b4>, w = <b0,b5>")
    println()
    r = chacha20(collect(UInt8, 0x00:0x37))
    nodd = 0
    for ki in 1:nkeys
        sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
        t = sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f)
        isodd(t) && (nodd += 1)
    end
    @printf("  %d keys from the reference key generator: %d have ODD ||(g,-f)||^2\n",
            nkeys, nodd)
    @printf("  -> part 2 of the countermeasure is %s on this key generator\n",
            nodd == 0 ? "UNAVAILABLE (every key is even)" : "available")
end

# ---------------------------------------------------------------------------
# mode: rate -- divergence rate with NewSamplerZ, on reference (even-t) keys
# ---------------------------------------------------------------------------
function respelled_tree(sk)
    t = with_spec_spelling(cdiv = true, ldl = true) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

function mode_rate(arm, nkeys, nsig, chunk, out)
    p = FALCON_512; n = p.n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)
    build, sign2 = arm == "A1" ? (identity, with_spec_ffsampling) :
                                 (respelled_tree, f -> f())
    r = chacha20(shake256(codeunits("cm-keys/$chunk"), 56))
    ev = 0; halfnear = 0
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        sk2 = build(sk)
        for j in 1:nsig
            st = shake256(codeunits("cm/$arm/$chunk/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            _, a = F.sample_preimage(sk,  pt, x -> randombytes!(r1, x))
            _, b = sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
            Int.(a) != Int.(b) && (ev += 1)
        end
    end
    open(out, "a") do io
        @printf(io, "%s %d %d %d %d\n", arm, chunk, nkeys, nsig, ev)
    end
    tot_ev = 0; tot_n = 0
    for line in eachline(out)
        f = split(line)
        length(f) >= 5 && f[1] == arm || continue
        tot_n += parse(Int, f[3]) * parse(Int, f[4]); tot_ev += parse(Int, f[5])
    end
    @printf("NewSamplerZ  arm %s chunk %d: %d events in %d sigs\n", arm, chunk, ev, nkeys*nsig)
    @printf("RUNNING TOTAL %s : %d events in %d signatures (rate %.3g)\n",
            arm, tot_ev, tot_n, tot_ev / max(tot_n, 1))
end

const MODE = length(ARGS) >= 1 ? ARGS[1] : "chisq"
if MODE == "chisq"
    mode_chisq()
elseif MODE == "parity"
    mode_parity(length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200)
elseif MODE == "rate"
    mode_rate(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]),
              parse(Int, ARGS[5]), ARGS[6])
else
    error("unknown mode $MODE")
end
