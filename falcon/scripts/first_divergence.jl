#!/usr/bin/env julia
#
# first_divergence.jl -- for one divergent (arm, key, signature) triple
# reported by scripts/divergence_rate.jl, find the first sampler call at which
# the two spellings disagree and print what differed there.
#
#     julia --project=falcon falcon/scripts/first_divergence.jl A1 4 105
#     julia --project=falcon falcon/scripts/first_divergence.jl A2 70 534
#
# The arm names, the key ordering and the PRNG seeds all match
# scripts/divergence_rate.jl exactly, so a case printed by that script can be
# handed straight to this one.
#
# WHAT TO LOOK FOR.  The mechanism (ePrint 2024/1709, Lemma 1) is that
# `SamplerZ` opens with `s = floor(mu)`, so a discrepancy of 1e-13 in the
# centre becomes a difference of 1 in `s` exactly when the two evaluations
# straddle an integer.  A divergence that is really this mechanism shows
#
#   * `mu` within about 1e-13 of an integer in both runs,
#   * `floor(mu)` differing between them,
#   * at call 1, 2, 2n-1 or 2n -- Heuristic 1 makes integer centres negligible
#     anywhere else, and scripts/heuristic1_check.jl confirms that directly
#     (0 in 102 million draws away from the ends).
#
# Position matters for more than tidiness.  Section 5 of that paper recovers
# the whole private key from a single pair differing in the LAST two calls; a
# difference in the first two yields only a short lattice vector, which is not
# enough.
#
# The two `samplerz` definitions below are copies of the library's, with a log
# added.  They are copies rather than wrappers because a wrapper that captured
# the original function and then redefined the name recursed into itself.

using Falcon, Printf
const F = Falcon

Falcon.eval(quote
    const LOG = NamedTuple{(:mu,:w,:z,:used),Tuple{Float64,Float64,Int,Int}}[]
    const ON = Ref(false)
    const USED = Ref(0)
    function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, rb)
        s = Int(floor(mu)); r = mu - s
        dss = 1/(2*sigma*sigma); ccs = sigmin/sigma
        u0 = USED[]
        while true
            z0 = basesampler(rb); USED[] += 9
            b = Int(rb(1)[1]) & 1; USED[] += 1
            z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            if berexp(x, ccs, rb)
                ON[] && push!(LOG, (mu=mu, w=sigma, z=z+s, used=USED[]-u0))
                return z + s
            end
        end
    end
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        u0 = USED[]
        while true
            z0 = basesampler(rb); USED[] += 9
            b = Int(rb(1)[1]) & 1; USED[] += 1
            z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            if berexp(x, ccs, rb)
                ON[] && push!(LOG, (mu=mu, w=isigma, z=z+s, used=USED[]-u0))
                return z + s
            end
        end
    end
end)

# --- the arms, kept identical to scripts/divergence_rate.jl ----------------

function respelled_tree(sk; cdiv::Bool, ldl::Bool)
    t = with_spec_spelling(cdiv = cdiv, ldl = ldl) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

function recip_tree(sk)
    recip = 1 / sk.params.sigma
    t = F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    patch!(x) = x isa F.FFLDLNode ? (patch!(x.left); patch!(x.right)) :
                (x.isigma = sqrt(real(x.value[1])) * recip)
    patch!(t)
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

const ARMS = Dict(
    "A1"  => (build = identity, sign2 = with_spec_ffsampling),
    "A2"  => (build = sk -> respelled_tree(sk, cdiv = true,  ldl = true),  sign2 = f -> f()),
    "A2a" => (build = sk -> respelled_tree(sk, cdiv = true,  ldl = false), sign2 = f -> f()),
    "A2b" => (build = sk -> respelled_tree(sk, cdiv = false, ldl = true),  sign2 = f -> f()),
    "B"   => (build = recip_tree, sign2 = f -> f()),
)

function main(arm::String, keyidx::Int, sigidx::Int)
    haskey(ARMS, arm) ||
        error("unknown arm $arm; expected one of " * join(sort(collect(keys(ARMS))), ", "))
    a = ARMS[arm]
    p = FALCON_512
    n = p.n

    r = chacha20(collect(UInt8, 0x00:0x37))
    local sk
    for _ in 1:keyidx
        sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
    end
    sk2 = a.build(sk)
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), 512; q = p.q)
    st = shake256(codeunits("$arm/$keyidx/$sigidx"), 56)

    F.ON[] = true
    empty!(F.LOG); F.USED[] = 0; r1 = chacha20(st)
    F.sample_preimage(sk, pt, x -> randombytes!(r1, x)); A = copy(F.LOG)
    empty!(F.LOG); F.USED[] = 0; r2 = chacha20(st)
    a.sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
    B = copy(F.LOG); F.ON[] = false

    @printf("arm %s  key %d  sig %d: %d vs %d sampler calls (2n = %d)\n",
            arm, keyidx, sigidx, length(A), length(B), 2n)
    for i in 1:min(length(A), length(B))
        A[i].z == B[i].z && continue
        pos = i == 1 || i == 2 ? "one of the FIRST two calls (short vector only)" :
              i >= 2n - 1     ? "one of the LAST two calls (full key recovery)" :
                                "the interior -- Heuristic 1 says this should not happen"
        @printf("first differing draw: call %d of %d -- %s\n", i, 2n, pos)
        @printf("  C-spelling  mu %.17g  w %.17g  z %d  bytes %d\n",
                A[i].mu, A[i].w, A[i].z, A[i].used)
        @printf("  spec        mu %.17g  w %.17g  z %d  bytes %d\n",
                B[i].mu, B[i].w, B[i].z, B[i].used)
        @printf("  |dmu| %.3g   distance to the nearest integer: C %.3g  spec %.3g\n",
                abs(A[i].mu - B[i].mu),
                abs(A[i].mu - round(A[i].mu)), abs(B[i].mu - round(B[i].mu)))
        @printf("  floor(mu): C %.17g  spec %.17g  -> %s\n",
                floor(A[i].mu), floor(B[i].mu),
                floor(A[i].mu) == floor(B[i].mu) ?
                    "EQUAL, so this is NOT the floor mechanism" : "differ (the floor straddle)")
        @printf("  widths equal: %s\n", A[i].w == B[i].w)
        k = count(j -> A[j].mu != B[j].mu, 1:i-1)
        @printf("  of the %d earlier calls, %d had a differing centre but the same z\n", i - 1, k)
        return
    end
    println("no differing draw found (?)")
end

length(ARGS) == 3 || error("usage: first_divergence.jl <arm> <key> <sig>")
main(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]))
