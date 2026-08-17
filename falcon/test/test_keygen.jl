# test_keygen.jl -- key generation (module 6 completed with module 7's sampler).
#
# The random stream is derived from a short label via SHAKE256 rather than
# embedded, so these vectors stay small; our SHAKE256 is itself checked against
# hashlib in test_shake.jl, so the derivation is not a weak link.
#
# The byte counts are the interesting part.  gen_poly draws 4096 samples
# whatever n is, and folds them; drawing n directly would give the same
# distribution and a different stream.  Recording how many bytes were consumed
# pins that -- a distribution test never could.

# Golden vectors are included by runtests.jl.

@testset "keygen" begin

    @testset "SIGMA_FG_MIN" begin
        @test SIGMA_FG_MIN == SIGMA_FG_BASE - 0.001
        # samplerz's precondition, with room to spare
        @test 1 < SIGMA_FG_MIN < SIGMA_FG_BASE < 1.8205
    end

    @testset "gen_poly against the reference, byte for byte" begin
        for (n, label, used, want) in GEN_POLY_KAT
            rb = ReplayBytes(shake256(codeunits(label), used))
            f = gen_poly(n, rb)
            @test f == want
            @test Falcon.consumed(rb) == used     # exact byte accounting
            @test length(f) == n
        end
    end

    @testset "gen_poly's fold is 4096 samples, not n" begin
        # Every degree consumes about the same number of bytes, because the
        # number of samplerz calls is 4096 regardless of n.  If someone
        # "optimised" this to draw n samples, the n = 8 case would consume
        # ~512x fewer bytes than n = 1024 -- so this is the assertion that
        # catches that change.
        used = [u for (_, _, u, _) in GEN_POLY_KAT]
        @test maximum(used) / minimum(used) < 1.05
        # and roughly 4096 samples at ~17 bytes each
        @test all(u -> 40_000 < u < 120_000, used)
    end

    @testset "gen_poly's output has the right width" begin
        # Effective sigma should be sigma_fg = 1.17*sqrt(q/(2n)), which is what
        # folding 4096/n samples of width SIGMA_FG_BASE gives.
        for (n, _, _, f) in GEN_POLY_KAT
            n < 64 && continue                   # too few coefficients to measure
            p = params(n == 512 ? 512 : 1024)
            want = n == 512 ? p.sigma_fg : (n == 1024 ? p.sigma_fg :
                                            1.17 * sqrt(Q / (2n)))
            v = sum(abs2, Float64.(f)) / length(f)
            @test isapprox(sqrt(v), want; rtol = 0.25)
        end
    end

    @testset "ntru_gen against the reference" begin
        for (n, label, used, wf, wg, wF, wG) in NTRU_GEN_KAT
            rb = ReplayBytes(shake256(codeunits(label), used))
            f, g, F, G = ntru_gen(n, rb)
            @test f == wf
            @test g == wg
            @test F == wF
            @test G == wG
            @test Falcon.consumed(rb) == used
        end
    end

    @testset "generated keys satisfy every condition keygen tested for" begin
        for (n, _, _, f, g, F, G) in NTRU_GEN_KAT
            # 1. the NTRU equation, exactly
            @test ntru_equation_holds(f, g, F, G)
            # 2. the Gram-Schmidt norm is within the acceptance bound
            @test gs_norm_ok(Float64.(f), Float64.(g))
            # 3. f is invertible mod q, so the public key h = g/f exists
            @test is_invertible_zq(Int.(f))
            h = polydivq(Int.(mod.(g, Q)), Int.(mod.(f, Q)))
            @test polymulq(h, Int.(mod.(f, Q))) == Int.(mod.(g, Q))
        end
    end

    @testset "keygen at the target dimension, self-validated" begin
        # n = 512 with no oracle: generate, then check the exact identity and
        # every acceptance condition.  Slow (the solve dominates) but this is
        # the dimension that matters.
        rng = chacha20(collect(UInt8, 0:55))
        src = bytesource(rng)
        f, g, F, G = ntru_gen(512, src)
        @test length(f) == length(g) == length(F) == length(G) == 512
        @test ntru_equation_holds(f, g, F, G)
        @test gs_norm_ok(Float64.(f), Float64.(g))
        @test is_invertible_zq(Int.(f))
        # (F, G) must be reduced, or the key is valid but useless
        fg_bits = max(maximum(bitsize, f), maximum(bitsize, g))
        @test maximum(bitsize, F) <= fg_bits + 16
        @test maximum(bitsize, G) <= fg_bits + 16
        # the public key round-trips
        fq = Int.(mod.(f, Q)); gq = Int.(mod.(g, Q))
        h = polydivq(gq, fq)
        @test polymulq(h, fq) == gq
    end

    @testset "gen_poly argument checking" begin
        rb = ReplayBytes(shake256(codeunits("genpoly-8"), 70001))
        @test_throws ArgumentError gen_poly(4096, rb)
        @test_throws ArgumentError gen_poly(3, rb)
    end
end
