#!/usr/bin/env julia
#
# window.jl -- the two-sided condition under which the first-two event channel
# of sections 6.4-6.5 has an oracle, measured rather than assumed.
#
#     julia --project=falcon falcon/scripts/window.jl [mode] [events] [sigs]
#     mode = eta | lower | upper | fma | noise | all   (default: all)
#
# Section 6.5 recovers the key from integer-centre EVENTS at sampler call 1,
# by Gaussian elimination over F_q.  It needs an oracle: something that tells
# the adversary which messages produced an event.  A difference between two
# implementations is such an oracle exactly when its perturbation |dmu| on the
# centre is neither too small nor too large:
#
#     |eta|  <<  |dmu|  <<  the false-positive ceiling
#
# LOWER EDGE.  At an event the true centre is an exact integer m, but both
# implementations compute m + eta with the SAME accumulated rounding error eta
# (measured here, mode `eta`).  floor() differs between them only if the two
# perturbed values land on opposite sides of m, so if |dmu| << |eta| the signs
# agree and nothing straddles.  This is why the respellings of section 5 --
# |dmu| = 2.7e-15 (A1) and exactly 0 (A2) at the first two calls, section 6.3 --
# produce no first-two straddles at all (10 manufactured events, 0 straddles,
# scripts/first_two_probe.jl).
#
# UPPER EDGE.  A perturbation of size dmu also straddles the ~2n INTERIOR
# centres, whose fractional parts are essentially uniform, at rate ~2n|dmu| per
# signature.  Those are false positives: the signatures differ but no event
# occurred.  The solve of section 6.5 needs n-1 rows that are all true, so it
# tolerates a false fraction of about 1/n; with a true-event rate of
# (1/q) * P(straddle | event) that gives
#
#     kappa * dmu  <=  (1/n) * (1/q) * P(straddle | event),
#
# with kappa the constant of the linear law rate = kappa*dmu (mode `upper`;
# scripts/precision_law.jl measures the same law over a coarser range).
#
# THE QUESTION.  Is the window non-empty, and does any REALISTIC difference land
# in it?  Mode `fma` answers the second half with the one difference that needs
# no disagreement between implementers at all: C99 6.5p8 lets the compiler
# contract a*b+c into a single fused multiply-add, GCC does so by default
# (-ffp-contract=fast) and clang does not.  Same source, same machine, two
# mainstream compilers at their defaults.  ePrint 2024/1709 section 6.2 raises
# FMA as a plausible divergence source; here it is an arm with a number on it.
#
# As in scripts/event_solve.jl the syndromes are drawn uniformly rather than
# hashed: hash_to_point's output is uniform and the measured event rate matches
# (section 6.5).  Using real hashes would only add SHAKE time.

using Falcon
using Printf
using Random
using Statistics

const F = Falcon

Falcon.eval(quote
    const MU    = Float64[]
    const REC   = Ref(false)
    const PERT  = Float64[]
    const PIDX  = Ref(0)
    @inline function _wpert()
        PIDX[] += 1
        (isempty(PERT) || PIDX[] > length(PERT)) ? 0.0 : PERT[PIDX[]]
    end
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        m = mu + _wpert()
        REC[] && push!(MU, m)
        s = Int(floor(m)); r = m - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        m = mu + _wpert()
        REC[] && push!(MU, m)
        s = Int(floor(m)); r = m - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

"Coefficient of f_j in (c*f)[k] over Z[x]/(x^n+1), reduced mod q."
function form_row(c::Vector{Int}, k::Int, n::Int, q::Int)
    row = Vector{Int}(undef, n)
    @inbounds for j in 0:(n-1)
        idx = k - j
        row[j+1] = idx >= 0 ? c[idx+1] : mod(-c[idx+n+1], q)
    end
    return row
end

"One signature under a given per-call perturbation vector; returns (s1, s2, centres)."
function run_arm(sk, pt, seed::Vector{UInt8}, pert::Vector{Float64}; fma::Bool = false)
    rb = chacha20(copy(seed))
    empty!(F.MU); empty!(F.PERT); append!(F.PERT, pert)
    F.PIDX[] = 0; F.REC[] = true
    s = fma ? with_fma(() -> F.sample_preimage(sk, pt, x -> randombytes!(rb, x))) :
              F.sample_preimage(sk, pt, x -> randombytes!(rb, x))
    F.REC[] = false; empty!(F.PERT)
    return (Int.(s[1]), Int.(s[2]), copy(F.MU))
end

"Draw uniform syndromes until `want` of them satisfy (c*f)[n/2-1] = 0 (mod q)."
function collect_events(rng, f::Vector{Int}, n::Int, q::Int, want::Int)
    k = n ÷ 2 - 1
    ev = Vector{Vector{Int}}(); tried = 0
    while length(ev) < want
        tried += 1
        c = rand(rng, 0:(q-1), n)
        row = form_row(c, k, n, q)
        v = 0
        @inbounds for j in 1:n; v = (v + row[j] * f[j]) % q; end
        mod(v, q) == 0 && push!(ev, c)
    end
    return ev, tried
end

pertvec(rng, L, d) = d == 0.0 ? Float64[] : Float64[d * (2rand(rng) - 1) for _ in 1:L]

# ---------------------------------------------------------------------------

function mode_eta(sk, f, n, q, nev)
    rng = MersenneTwister(11)
    ev, tried = collect_events(rng, f, n, q, nev)
    @printf("\n## eta -- the common rounding error at an exact-integer centre\n")
    @printf("# %d events from %d uniform syndromes (rate %.3g, 1/q = %.3g)\n",
            length(ev), tried, length(ev)/tried, 1/q)
    etas = Float64[]
    for (i, c) in enumerate(ev)
        _, _, mu = run_arm(sk, c, shake256(codeunits("w/eta/$i"), 56), Float64[])
        push!(etas, mu[1] - round(mu[1]))
    end
    a = abs.(etas)
    @printf("|eta| at call 1 over %d events:  min %.3g  median %.3g  max %.3g\n",
            length(a), minimum(a), median(a), maximum(a))
    @printf("sign(eta) > 0 in %d of %d  (a straddle needs the two arms to differ in sign)\n",
            count(>(0), etas), length(etas))
    return median(a), ev
end

function mode_lower(sk, ev, n; deltas = [0.0, 1e-16, 1e-15, 1e-14, 1e-13,
                                          1e-12, 1e-11, 1e-10, 1e-8])
    L = 2n
    @printf("\n## lower -- P(straddle at call 1 | exact integer centre) vs perturbation\n")
    @printf("# %d events; arm A unperturbed, arm B perturbed by delta*U[-1,1] at every call\n",
            length(ev))
    @printf("%-10s  %8s  %8s  %11s  %s\n", "delta", "straddle", "sigdiff", "P(straddle)", "95% CI")
    for d in deltas
        rng = MersenneTwister(2026)
        ns = 0; nd = 0
        for (i, c) in enumerate(ev)
            seed = shake256(codeunits("w/low/$i"), 56)
            p = pertvec(rng, L, d)
            a1, a2, mua = run_arm(sk, c, seed, Float64[])
            b1, b2, mub = run_arm(sk, c, seed, p)
            floor(mua[1]) != floor(mub[1]) && (ns += 1)
            (a1 != b1 || a2 != b2) && (nd += 1)
        end
        m = length(ev)
        ph = ns/m; se = sqrt(max(ph*(1-ph), 1e-12)/m)
        @printf("%-10.3g  %8d  %8d  %11.3f  %.3f-%.3f\n", d, ns, nd, ph,
                max(0.0, ph-1.96se), min(1.0, ph+1.96se))
    end
end

function mode_upper(sk, n, q, nsig; deltas = [1e-4, 3e-5, 1e-5, 3e-6])
    L = 2n
    rng = MersenneTwister(77)
    @printf("\n## upper -- the interior false-positive law, rate = kappa * delta\n")
    @printf("# %d random syndromes per delta (events are ~1/q of them, negligible here)\n", nsig)
    @printf("%-10s  %8s  %-11s  %-11s\n", "delta", "diffs", "rate", "kappa=rate/delta")
    ks = Float64[]
    for d in deltas
        nd = 0
        for j in 1:nsig
            c = rand(rng, 0:(q-1), n)
            seed = shake256(codeunits("w/up/$j"), 56)
            p = pertvec(rng, L, d)
            a1, a2, _ = run_arm(sk, c, seed, Float64[])
            b1, b2, _ = run_arm(sk, c, seed, p)
            (a1 != b1 || a2 != b2) && (nd += 1)
        end
        r = nd/nsig; push!(ks, r/d)
        @printf("%-10.3g  %8d  %-11.3g  %-11.3g\n", d, nd, r, r/d)
    end
    return median(ks)
end

function mode_fma(sk, f, ev, n, q, nsig)
    L = 2n
    rng = MersenneTwister(5)
    @printf("\n## fma -- GCC's -ffp-contract=fast against clang's default, as an arm\n")
    d1 = Float64[]; dmid = Float64[]; dlast = Float64[]
    nd = 0
    for j in 1:nsig
        c = rand(rng, 0:(q-1), n)
        seed = shake256(codeunits("w/fma/$j"), 56)
        a1, a2, mua = run_arm(sk, c, seed, Float64[])
        b1, b2, mub = run_arm(sk, c, seed, Float64[]; fma = true)
        (a1 != b1 || a2 != b2) && (nd += 1)
        length(mua) == L && length(mub) == L || continue
        push!(d1, abs(mua[1] - mub[1])); push!(d1, abs(mua[2] - mub[2]))
        push!(dlast, abs(mua[L-1] - mub[L-1])); push!(dlast, abs(mua[L] - mub[L]))
        for t in 3:(L-2); push!(dmid, abs(mua[t] - mub[t])); end
    end
    @printf("# %d random syndromes, both arms on the same PRNG tape\n", nsig)
    @printf("|dmu| first two : median %.3g  max %.3g  (exactly zero in %d of %d)\n",
            median(d1), maximum(d1), count(==(0.0), d1), length(d1))
    @printf("|dmu| interior  : median %.3g  max %.3g\n", median(dmid), maximum(dmid))
    @printf("|dmu| last two  : median %.3g  max %.3g\n", median(dlast), maximum(dlast))
    @printf("signature divergence rate: %d/%d = %.3g\n", nd, nsig, nd/nsig)

    ns = 0; nsd = 0
    for (i, c) in enumerate(ev)
        seed = shake256(codeunits("w/fmae/$i"), 56)
        a1, a2, mua = run_arm(sk, c, seed, Float64[])
        b1, b2, mub = run_arm(sk, c, seed, Float64[]; fma = true)
        floor(mua[1]) != floor(mub[1]) && (ns += 1)
        (a1 != b1 || a2 != b2) && (nsd += 1)
    end
    @printf("on %d manufactured events: %d straddled call 1, %d signatures differed\n",
            length(ev), ns, nsd)
    return ns/max(length(ev),1), median(d1)
end

"""
How informative the one-bit oracle is: of the signatures on which the two arms
disagree, what fraction carries a first-two EVENT (an F_q equation) and what
fraction is a last-two divergence (which carries no equation, and which ePrint
2024/1709 section 5.1 already exploits by a different route)?

Scanning costs one arm only: a near-integer centre is detected from arm A, and
arm B is run just on the hits.
"""
function mode_noise(sk, n, q, npair)
    L = 2n
    rng = MersenneTwister(31)
    @printf("\n## noise -- what fraction of FMA disagreements carries an equation\n")
    h1 = 0; s1 = 0; h2 = 0; s2 = 0
    for j in 1:npair
        c = rand(rng, 0:(q-1), n)
        seed = shake256(codeunits("w/nz/$j"), 56)
        a1, a2, mua = run_arm(sk, c, seed, Float64[])
        length(mua) == L || continue
        first = abs(mua[1] - round(mua[1])) < 1e-9 || abs(mua[2] - round(mua[2])) < 1e-9
        last  = abs(mua[L-1] - round(mua[L-1])) < 1e-9 || abs(mua[L] - round(mua[L])) < 1e-9
        (first || last) || continue
        b1, b2, mub = run_arm(sk, c, seed, Float64[]; fma = true)
        if first
            h1 += 1
            (floor(mua[1]) != floor(mub[1]) || floor(mua[2]) != floor(mub[2])) && (s1 += 1)
        end
        if last
            h2 += 1
            length(mub) == L &&
                (floor(mua[L-1]) != floor(mub[L-1]) || floor(mua[L]) != floor(mub[L])) && (s2 += 1)
        end
    end
    @printf("# %d syndromes scanned with arm A; arm B run only on the hits\n", npair)
    @printf("first two : %d integer centres (rate %.3g, 2/q = %.3g), %d straddled -> informative rate %.3g\n",
            h1, h1/npair, 2/q, s1, s1/npair)
    @printf("last two  : %d integer centres (rate %.3g), %d straddled -> uninformative rate %.3g\n",
            h2, h2/npair, s2, s2/npair)
    tot = s1 + s2
    @printf("of the disagreements the one-bit oracle reports, %d of %d = %.0f%% carry an F_q equation\n",
            s1, tot, tot == 0 ? NaN : 100*s1/tot)
    return s1, s2
end

function main()
    mode = length(ARGS) >= 1 ? ARGS[1] : "all"
    nev  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    nsig = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1500
    p = FALCON_512; n = p.n; q = p.q
    r = chacha20(collect(UInt8, 0x00:0x37))
    sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
    f = Int.(sk.f)
    @printf("# window.jl  n=%d  q=%d  mode=%s  events=%d  sigs=%d\n", n, q, mode, nev, nsig)

    etam = NaN; ev = Vector{Vector{Int}}(); kappa = NaN; pfma = NaN; dfma = NaN
    if mode in ("all", "eta", "lower", "fma")
        etam, ev = mode_eta(sk, f, n, q, nev)
    end
    mode in ("all", "lower") && mode_lower(sk, ev, n)
    if mode in ("all", "upper")
        kappa = mode_upper(sk, n, q, nsig)
    end
    if mode in ("all", "fma")
        pfma, dfma = mode_fma(sk, f, ev, n, q, nsig)
    end
    mode == "noise" && mode_noise(sk, n, q, nsig)

    if mode == "all"
        @printf("\n## the window\n")
        fpmax = (1/n) * (1/q) * 0.5
        dmax = fpmax / kappa
        @printf("lower edge  ~ |eta|                 = %.3g\n", etam)
        @printf("upper edge  = (1/n)(1/q)(1/2)/kappa = %.3g   (kappa = %.3g)\n", dmax, kappa)
        @printf("width                               = %.2f decades\n", log10(dmax/etam))
        @printf("respellings (sec 6.3): |dmu| = 2.7e-15 (A1), 0 (A2)  -> below the lower edge\n")
        @printf("FMA contraction      : |dmu| = %.3g, straddle rate %.3f  -> %s\n",
                dfma, pfma, (dfma > etam/10 && dfma < dmax) ? "INSIDE the window" :
                            (dfma >= dmax ? "above the window" : "below the window"))
    end
end

main()
