#!/usr/bin/env julia
#
# position_profile.jl -- per sampler-call-position profile of (i) how close the
# centre gets to an integer and (ii) how large a perturbation the A2 respelling
# (complex division + D11) puts on that centre.  One instrument, two questions:
#
#   #2  Is the sensitive set really only positions {0,1,2n-2,2n-1}?  We measure
#       the near-integer margin at EVERY position, so a fifth sensitive position
#       would show up as an interior position with a small margin.
#   #4  ePrint 2024/1709 sec 6.1 leaves open why, conditional on an integer
#       centre, the LAST two calls diverge more readily than the first two, for
#       the dyn/tree difference.  If the respelling perturbs the last-two
#       centres by a larger |dmu| than the first-two, that is the explanation.
#
#     julia --project=falcon falcon/scripts/position_profile.jl [keys] [sigs]
#
# For each position i in 0..2n-1 we accumulate, over many (key, signature):
#   * min |mu - round(mu)|                 -- the near-integer margin
#   * count of |mu - round(mu)| < 1e-6     -- near-integer events
#   * mean and max |mu_C - mu_spec|        -- the respelling perturbation
# where mu_C is the C spelling and mu_spec is with cdiv+ldl in spec form, on
# identical key, message and PRNG state.

using Falcon
using Printf
using Statistics

const F = Falcon

# log (mu) at every call, for whichever run is active
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
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        ON[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

function respelled_tree(sk)
    t = with_spec_spelling(cdiv = true, ldl = true) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 40
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    arm   = length(ARGS) >= 3 ? ARGS[3] : "A2"     # A2 (tree respelling) or A1 (bottom-level)
    p = FALCON_512; n = p.n; L = 2n

    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)
    r = chacha20(collect(UInt8, 0x00:0x37))

    minfrac = fill(Inf, L)
    nearcount = zeros(Int, L)
    sumdmu = zeros(Float64, L)
    maxdmu = zeros(Float64, L)
    ndraw = 0

    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        sk2 = arm == "A2" ? respelled_tree(sk) : sk   # A1 uses same tree, spec ffsampling
        for j in 1:nsig
            st = shake256(codeunits("pp/$arm/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            empty!(F.MU); F.ON[] = true
            F.sample_preimage(sk, pt, x -> randombytes!(r1, x)); A = copy(F.MU)
            empty!(F.MU)
            if arm == "A1"
                with_spec_ffsampling(() -> F.sample_preimage(sk, pt, x -> randombytes!(r2, x)))
            else
                F.sample_preimage(sk2, pt, x -> randombytes!(r2, x))
            end
            B = copy(F.MU)
            F.ON[] = false
            (length(A) == L && length(B) == L) || continue
            ndraw += 1
            for i in 1:L
                fr = abs(A[i] - round(A[i]))
                fr < minfrac[i] && (minfrac[i] = fr)
                fr < 1e-6 && (nearcount[i] += 1)
                d = abs(A[i] - B[i])
                sumdmu[i] += d
                d > maxdmu[i] && (maxdmu[i] = d)
            end
        end
    end

    @printf("# position_profile.jl  arm=%s  n=%d  %d keys x %d sigs = %d draws/position\n",
            arm, n, nkeys, nsig, ndraw)
    println()
    show_positions = vcat(1:6, (L-5):L)
    println("# per-position: min |frac|, near(<1e-6) count, mean|dmu|, max|dmu|")
    @printf("%6s  %-12s  %6s  %-11s  %-11s\n", "pos", "min|frac|", "near", "mean|dmu|", "max|dmu|")
    for i in show_positions
        tag = (i <= 2) ? " <- first two" : (i >= L-1) ? " <- last two" : ""
        @printf("%6d  %-12.3e  %6d  %-11.3e  %-11.3e%s\n",
                i-1, minfrac[i], nearcount[i], sumdmu[i]/max(ndraw,1), maxdmu[i], tag)
    end
    println()
    # interior summary (positions 7 .. 2n-6)
    interior = 7:(L-6)
    @printf("# interior positions %d..%d: min|frac| over all = %.3e, total near-events = %d\n",
            first(interior)-1, last(interior)-1,
            minimum(minfrac[interior]), sum(nearcount[interior]))
    println()
    # accumulation shape: mean |dmu| over 16 blocks across all 2n positions,
    # to see whether the perturbation grows along the traversal (issue #4 mechanism)
    println("# accumulation of |dmu| across the traversal (16 blocks of ", L÷16, " positions):")
    for blk in 0:15
        lo = blk*(L÷16) + 1; hi = (blk+1)*(L÷16)
        m = mean(sumdmu[lo:hi]) / max(ndraw,1)
        @printf("  pos %4d..%4d : mean|dmu| %.3e  %s\n", lo-1, hi-1, m, "#"^clamp(round(Int, log10(m+1e-300)+18), 0, 40))
    end
    println()
    # the #4 question, quantified:
    firsttwo_dmu = (sumdmu[1] + sumdmu[2]) / (2 * max(ndraw,1))
    lasttwo_dmu  = (sumdmu[L-1] + sumdmu[L]) / (2 * max(ndraw,1))
    @printf("# ASYMMETRY (issue #4): mean |dmu|  first-two %.3e   last-two %.3e   ratio %.2f\n",
            firsttwo_dmu, lasttwo_dmu, lasttwo_dmu / max(firsttwo_dmu, 1e-300))
    @printf("#   max |dmu|          first-two %.3e   last-two %.3e\n",
            max(maxdmu[1],maxdmu[2]), max(maxdmu[L-1],maxdmu[L]))
end

main()
