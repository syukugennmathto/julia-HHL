using Statistics, Random
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
            # `sampler = :spec` selects the specification's gen_poly.  This
            # branch defaults to the C reference's CDT sampler, which produces
            # the same distribution by a different route and so cannot
            # reproduce a recorded vector (README.md, docs/debug_log.md #043).
            f, g, F, G = ntru_gen(n, rb; sampler = :spec)
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

    @testset "the CDT sampler: same distribution, different realisation" begin
        # THIS BRANCH's default sampler (docs/debug_log.md #043).  It cannot be
        # checked against a recorded vector -- that is the whole point of the
        # divergence -- so it is checked against the distribution it is supposed
        # to have, and against the constraints it is supposed to enforce.

        rng = MersenneTwister(20260827)

        # 1. The table is the one the C reference carries, and it says what its
        #    header says: entry 0 is P(x = 0), entry k is P(x >= k+1 | x > 0),
        #    scaled by 2^63, for sigma = 1.17*sqrt(q/(2*1024)).
        @test length(GAUSS_1024_12289) == 27
        @test GAUSS_1024_12289[end] == 0            # the table terminates
        @test issorted(GAUSS_1024_12289[2:end]; rev = true)
        let top = big(2)^63, sigma = 1.17 * sqrt(Q / 2048)
            p0 = Float64(GAUSS_1024_12289[1] / top)
            # P(x=0) for a discrete Gaussian of this width
            norm = sum(exp(-k^2 / (2 * sigma^2)) for k in -60:60)
            @test isapprox(p0, 1 / norm; rtol = 1e-3)
            # P(x=1 | x>0) follows from the same distribution
            p1given = 1 - Float64(GAUSS_1024_12289[2] / top)
            want = exp(-1 / (2 * sigma^2)) /
                   sum(exp(-k^2 / (2 * sigma^2)) for k in 1:60)
            @test isapprox(p1given, want; rtol = 1e-3)
        end

        # 2. Summing 2^(10-logn) draws gives the specification's sigma_fg, which
        #    is what makes this a *different route to the same distribution*
        #    rather than a different distribution.
        for n in (512, 1024)
            src = bytesource(chacha20(collect(UInt8, 0x20:0x57)))
            v = Int[]
            for _ in 1:30
                append!(v, Int.(gen_poly_cdt(n, src)))
            end
            want = SIGMA_FG_BASE * sqrt(4096 / n)
            @test isapprox(std(v), want; rtol = 0.05)
            @test abs(mean(v)) < 0.2
            @test all(c -> -127 <= c <= 127, v)      # fits the key format's byte
        end

        # 3. The constraint the specification's sampler does *not* impose: the
        #    coefficient sum is odd, so Res(f, x^n+1) is odd and the binary GCD
        #    at the bottom of the descent cannot fail on a factor of two.
        let src = bytesource(chacha20(collect(UInt8, 0x60:0x97)))
            for _ in 1:50
                @test isodd(sum(Int.(gen_poly_cdt(512, src))))
            end
        end

        # 3b. The buffered reader the sampler draws through hands out the same
        #     64-bit words a naive eight-bytes-at-a-time reader would.  It has
        #     to: 512 is the ChaCha20 buffer size and is divisible by 8, so
        #     neither request pattern triggers the reference's end-of-buffer
        #     discard, and the streams coincide byte for byte.
        #
        #     This is the regression test for #043 -- not for the boxing (that
        #     was a speed bug, and speed is not a test) but for the rewrite
        #     that fixed it, which touched how bytes are read.
        let a = chacha20(collect(UInt8, 0x01:0x38)),
            b = chacha20(collect(UInt8, 0x01:0x38))
            buffered = Falcon.BufferedU64(k -> randombytes!(a, k))
            naive = () -> begin
                v = randombytes!(b, 8)
                r = UInt64(0)
                for i in 1:8
                    r |= UInt64(v[i]) << (8 * (i - 1))
                end
                r
            end
            # 1000 words spans several refills of the 64-word block
            @test all(buffered() == naive() for _ in 1:1000)
        end

        # 4. `sampler = :spec` still selects the specification's route, which is
        #    what keeps the recorded vectors above meaningful.
        @test_throws ArgumentError ntru_gen(8, ReplayBytes(zeros(UInt8, 8)); sampler = :nope)
        let n = 8
            rb1 = ReplayBytes(shake256(codeunits("sampler/spec"), 1 << 18))
            rb2 = ReplayBytes(shake256(codeunits("sampler/spec"), 1 << 18))
            a = ntru_gen(n, rb1; sampler = :spec)
            b = ntru_gen(n, rb2; sampler = :spec)
            @test a == b                              # deterministic
            @test ntru_equation_holds(a[1], a[2], a[3], a[4])
        end

        # 5. Keys from the fast path are real keys.
        let src = bytesource(chacha20(collect(UInt8, 0xa0:0xd7)))
            f, g, F, G = ntru_gen(512, src)
            @test ntru_equation_holds(f, g, F, G)
            @test gs_norm_ok(Float64.(f), Float64.(g))
            @test is_invertible_zq(Int.(f))
            @test isodd(sum(Int.(f))) && isodd(sum(Int.(g)))
        end
    end
end
