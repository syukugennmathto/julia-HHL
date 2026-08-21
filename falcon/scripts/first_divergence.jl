#!/usr/bin/env julia
#
# first_divergence.jl -- for one divergent (key, signature) pair reported by
# scripts/divergence_rate.jl, find the first sampler call at which the two
# spellings disagree and print what differed there.
#
#     julia --project=falcon falcon/scripts/first_divergence.jl 6 764
#
# Every case measured so far is a nearly-integer centre at call 1023 or 1024
# of 1024 -- positions 2n-2 and 2n-1, exactly where ePrint 2024/1709 predicts
# the sensitivity.  `SamplerZ` starts with `s = floor(mu)`, and `floor` is
# discontinuous, so a 1e-13 discrepancy in `mu` becomes a difference of 1 in
# `s`.  See scripts/divergence_rate.jl for the full argument.
#
# The two `samplerz` definitions below are copies of the library's, with a log
# added.  They are copies rather than wrappers because a wrapper that captured
# the original function and then redefined the name recursed into itself.

# Find the first sampler call at which the two spellings diverge, and say what
# differed there.
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
function main(keyidx, sigidx)
    p = FALCON_512
    r = chacha20(collect(UInt8, 0x00:0x37))
    local sk
    for _ in 1:keyidx; sk = falcon_keygen(512, k -> randombytes!(r,k))[1]; end
    msg = collect(codeunits("power analysis message"))
    pt = hash_to_point(msg, shake256(codeunits("power salt"), SALT_LEN), 512; q=p.q)
    st = shake256(codeunits("prng/$keyidx/$sigidx"), 56)
    F.ON[] = true
    empty!(F.LOG); F.USED[] = 0; r1 = chacha20(st)
    F.sample_preimage(sk, pt, x -> randombytes!(r1, x)); A = copy(F.LOG)
    empty!(F.LOG); F.USED[] = 0; r2 = chacha20(st)
    with_spec_ffsampling(() -> F.sample_preimage(sk, pt, x -> randombytes!(r2, x)))
    B = copy(F.LOG); F.ON[] = false
    @printf("key %d sig %d: %d vs %d sampler calls\n", keyidx, sigidx, length(A), length(B))
    for i in 1:min(length(A), length(B))
        if A[i].z != B[i].z
            @printf("first differing draw: call %d\n", i)
            @printf("  C-spelling  mu %.17g  w %.17g  z %d  bytes %d\n", A[i].mu, A[i].w, A[i].z, A[i].used)
            @printf("  spec        mu %.17g  w %.17g  z %d  bytes %d\n", B[i].mu, B[i].w, B[i].z, B[i].used)
            @printf("  |dmu| %.3g   mu - floor(mu): C %.17g  spec %.17g\n",
                    abs(A[i].mu-B[i].mu), A[i].mu-floor(A[i].mu), B[i].mu-floor(B[i].mu))
            @printf("  floor(mu) equal: %s   1/w(C) %.17g vs w(spec) %.17g\n",
                    floor(A[i].mu)==floor(B[i].mu), 1/A[i].w, B[i].w)
            # how many draws before this one already had differing mu?
            k = count(j -> A[j].mu != B[j].mu, 1:i-1)
            @printf("  of the %d earlier calls, %d had a differing centre but the same z\n", i-1, k)
            return
        end
    end
    println("no differing draw found (?)")
end
main(parse(Int,ARGS[1]), parse(Int,ARGS[2]))
