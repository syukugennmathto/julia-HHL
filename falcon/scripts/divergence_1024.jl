#!/usr/bin/env julia
#
# divergence_1024.jl -- the A2 (complex division + D11) and A1 (hand-unrolled
# bottom levels) arms of scripts/divergence_rate.jl, at FALCON-1024.
#
#     julia --project=falcon falcon/scripts/divergence_1024.jl [keys] [sigs]
#
# The mechanism (ePrint 2024/1709, Lemma 1) is degree-independent: an integer
# centre passing through `floor`.  What changes with degree is the RATE, since
# Heuristic 1 puts the last-two-call probability at 1/||(g,-f)||^2, and that
# norm grows with n.  So this is a rate measurement, not a new phenomenon.
#
# The only n=1024-specific care: samplerz's INV_2SIGMA2 constant in the
# library is sigma-max-dependent and already correct for both degrees, and
# with_spec_spelling / respelled_tree are degree-generic, so the arms are the
# same code as divergence_rate.jl with FALCON_1024 substituted for FALCON_512.

using Falcon
using Printf

const F = Falcon

point() = hash_to_point(collect(codeunits("power analysis message")),
                        shake256(codeunits("power salt"), SALT_LEN), 1024;
                        q = FALCON_1024.q)

function respelled_tree(sk; cdiv::Bool, ldl::Bool)
    t = with_spec_spelling(cdiv = cdiv, ldl = ldl) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

function rate(nkeys, nsig, tag; build = identity, sign2 = f -> f())
    pt = point()
    r = chacha20(collect(UInt8, 0x00:0x37))
    n = 0; ndiff = 0; ncoef = 0; ncoefdiff = 0
    divergent = Tuple{Int,Int,Int}[]
    for ki in 1:nkeys
        sk = falcon_keygen(1024, k -> randombytes!(r, k))[1]
        sk2 = build(sk)
        for j in 1:nsig
            st = shake256(codeunits("$tag/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
            _, b = sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
            n += 1
            d = count(Int.(a) .!= Int.(b))
            ncoef += length(a); ncoefdiff += d
            d > 0 && (ndiff += 1; push!(divergent, (ki, j, d)))
        end
    end
    return (n, ndiff, ncoef, ncoefdiff, divergent)
end

function report(title, res)
    n, nd, nc, ncd, div = res
    println("## ", title)
    @printf("  signatures differing   : %d of %d   (rate %.3g)\n", nd, n, nd / n)
    @printf("  coefficients differing : %d of %d\n", ncd, nc)
    if nd == 0
        @printf("  zero events -> 95%% upper bound on the rate: %.3g\n", 3 / n)
    else
        @printf("  when one differs, %.0f of 1024 coefficients do\n", ncd / nd)
        for d in div
            @printf("    %s\n", string(d))
        end
    end
    println()
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 50
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
    println("# divergence_1024.jl  n=1024  ", nkeys, " keys x ", nsig, " signatures per arm")

    # ||(g,-f)||^2 for context: Heuristic 1's last-two-call denominator
    let r0 = chacha20(collect(UInt8, 0x00:0x37)),
        sk0 = falcon_keygen(1024, k -> randombytes!(r0, k))[1]
        nsq = sum(x -> Int128(x)^2, sk0.g) + sum(x -> Int128(x)^2, sk0.f)
        @printf("# ||(g,-f)||^2 (one key) = %d   (1/that = %.3g)\n", nsq, 1 / nsq)
    end
    println()

    report("A2 -- complex division + D11 (not in 2024/1709)",
           rate(nkeys, nsig, "A2_1024"; build = sk -> respelled_tree(sk, cdiv = true, ldl = true)))
    report("A1 -- hand-unrolled bottom levels (= 2024/1709 6.1, positive control)",
           rate(nkeys, nsig, "A1_1024"; sign2 = with_spec_ffsampling))
end

main()
