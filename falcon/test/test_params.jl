# test_params.jl -- module 1.
#
# There is nothing to "run" in params.jl, so these tests check the only things
# that can go wrong with a table of constants: a mistyped digit, and a value
# that is inconsistent with the formula it is supposed to come from.

@testset "params" begin

    @testset "modulus" begin
        @test Q == 12289
        @test Q == 12 * 1024 + 1
        @test Q == 3 * 2^12 + 1
        # q - 1 must be divisible by 2n for the NTT to exist at n = 1024.
        @test (Q - 1) % (2 * 1024) == 0
        @test isprime_trial(Q)
    end

    @testset "degree-independent lengths" begin
        @test HEAD_LEN == 1
        @test SALT_LEN == 40
        @test SEED_LEN == 56
    end

    @testset "byte lengths agree with the C reference macros" begin
        for p in (FALCON_512, FALCON_1024)
            @test p.pubkey_bytes == Falcon.c_pubkey_size(p.logn)
            @test p.privkey_bytes == Falcon.c_privkey_size(p.logn)
            @test p.sig_bytes == Falcon.c_sig_padded_size(p.logn)
        end
        # Spot values, independently computed by hand from falcon.h:
        @test FALCON_512.pubkey_bytes == 897
        @test FALCON_512.privkey_bytes == 1281
        @test FALCON_512.sig_bytes == 666
        @test FALCON_1024.pubkey_bytes == 1793
        @test FALCON_1024.privkey_bytes == 2305
        @test FALCON_1024.sig_bytes == 1280

        # The macros must be monotone in logn -- falcon.h states this
        # explicitly ("increasing logn cannot result in a shorter length"),
        # and it is a cheap check that the shifts were transcribed correctly
        # at *every* degree, not just the two we use.
        for logn in 2:9
            @test Falcon.c_pubkey_size(logn) <= Falcon.c_pubkey_size(logn + 1)
            @test Falcon.c_privkey_size(logn) <= Falcon.c_privkey_size(logn + 1)
            @test Falcon.c_sig_padded_size(logn) <= Falcon.c_sig_padded_size(logn + 1)
        end

        # The PADDED size must leave room for header + salt + at least one
        # byte of compressed signature.
        for p in (FALCON_512, FALCON_1024)
            @test p.sig_bytes > HEAD_LEN + SALT_LEN
        end
    end

    @testset "sigma_min is the smoothing parameter of Z" begin
        # [derived]: sigma_min = eta_eps(Z) at eps = 1/sqrt(2^64 n^3).
        # The reference literals and this closed form are evaluated by
        # different means, so we only require agreement to ~1e-11 relative.
        for p in (FALCON_512, FALCON_1024)
            eta = smoothing_eta(falcon_eps(p.n))
            @test isapprox(eta, p.sigma_min; rtol = 1e-11)
        end
        # The exponent pattern that identified the formula in the first place:
        # log2(1/eps) = 32 + 1.5*logn, i.e. 45.5 and 47.0.
        @test isapprox(log2(1 / falcon_eps(512)), 45.5; atol = 1e-9)
        @test isapprox(log2(1 / falcon_eps(1024)), 47.0; atol = 1e-9)
    end

    @testset "sigma = 1.17 * sqrt(q) * sigma_min" begin
        for p in (FALCON_512, FALCON_1024)
            @test isapprox(gram_schmidt_quality() * sqrt(p.q) * p.sigma_min,
                           p.sigma; rtol = 1e-11)
        end
    end

    @testset "sigma ordering" begin
        for p in (FALCON_512, FALCON_1024)
            # The per-node width handed to samplerz must land in
            # [sigma_min, sigma_max]; that only makes sense if the interval is
            # non-empty, and sigma_max must be the CDT base sampler's width.
            @test p.sigma_min < p.sigma_max
            @test p.sigma_max == 1.8205
            # The lattice width is far larger than the per-coefficient widths.
            # (The ratio is about 91 at both degrees.  The first version of
            # this test asserted > 100 * sigma_max, which is simply false --
            # a wrong test, not a wrong constant.  debug_log #018.)
            @test p.sigma > 50 * p.sigma_max
            @test 80 < p.sigma / p.sigma_max < 100
            # f, g are much narrower than the signature Gaussian.
            @test p.sigma_fg < p.sigma
        end
        # sigma and sigma_min both grow with n; sigma_fg shrinks (it is
        # 1.17*sqrt(q/2n)).
        @test FALCON_512.sigma < FALCON_1024.sigma
        @test FALCON_512.sigma_min < FALCON_1024.sigma_min
        @test FALCON_512.sigma_fg > FALCON_1024.sigma_fg
    end

    @testset "sigma_fg and its base" begin
        # The reference draws 4096 samples at SIGMA_FG_BASE and folds them in
        # blocks of k = 4096/n; the effective width is sqrt(k)*SIGMA_FG_BASE.
        for p in (FALCON_512, FALCON_1024)
            k = 4096 ÷ p.n
            @test isapprox(sqrt(k) * SIGMA_FG_BASE, p.sigma_fg; rtol = 1e-12)
        end
        @test isapprox(SIGMA_FG_BASE, 1.17 * sqrt(Q / 8192); rtol = 1e-12)
    end

    @testset "signature bound" begin
        # Both reference implementations agree on these, and they are
        # *inclusive* bounds equal to floor(beta^2).
        @test FALCON_512.sig_bound == 34034726
        @test FALCON_1024.sig_bound == 70265242

        # beta^2 should be within an order of magnitude of the expected
        # squared norm of an honest signature, 2n * sigma^2.  This is the
        # check that would catch a bound copied from the wrong row of the
        # table: the ratio is a small constant, not a factor of 1000.
        for p in (FALCON_512, FALCON_1024)
            expected = 2 * p.n * p.sigma^2
            @test 1.0 < p.sig_bound / expected < 2.0
        end
    end

    @testset "lookup" begin
        @test params(512) === FALCON_512
        @test params(1024) === FALCON_1024
        @test_throws ArgumentError params(256)
        @test_throws ArgumentError params(0)
    end

    @testset "provenance placeholders are honest" begin
        # These used to say TODO rather than carry an invented table number,
        # and this test asserted the TODO -- as a reminder that the PDF was
        # still missing, and to make removing the placeholder a deliberate act
        # rather than something that could drift.
        #
        # The PDF arrived on 2026-08-21 (docs/debug_log.md #046).  So the test
        # is inverted: the fields must now name the document and the table,
        # and test/test_spec.jl checks that the *values* are what that table
        # says.  A `spec_ref` naming a table nothing verifies would be back to
        # an invented citation, which is the thing being guarded against.
        for p in (FALCON_512, FALCON_1024)
            @test !occursin("TODO", p.spec_ref)
            @test occursin("Falcon spec v1.2", p.spec_ref)
            @test occursin("Table 3.3", p.spec_ref)
        end
        @test occursin("Falcon-512", FALCON_512.spec_ref)
        @test occursin("Falcon-1024", FALCON_1024.spec_ref)
    end
end
