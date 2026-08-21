# test_samplerz.jl -- module 7.
#
# A sampler is harder to test than anything so far, because its output is
# *supposed* to look random.  A subtly wrong sampler returns plausible integers
# forever and only shows up as a slightly-too-long signature, or a slow leak of
# the key.  So three independent kinds of check:
#
#   1. **The official KAT.**  3072 vectors of (mu, sigma, sigmin, random bytes,
#      expected z).  Byte-exact: the sampler must consume the same randomness in
#      the same order and reach the same answer.  This is the strong one -- it
#      pins the fixed-point arithmetic, the constants, and the byte accounting
#      all at once.
#   2. **The constants against what they claim to be.**  RCDT is supposed to be
#      the reverse CDF of a half-Gaussian at width 1.8205, and EXP_COEFFS is
#      supposed to approximate exp(-x).  Both are checkable against the
#      continuous functions they came from, which is what catches a transposed
#      digit that the KAT would also catch but less legibly.
#   3. **The distribution itself.**  Mean, variance and a chi-squared test over
#      a large sample, driven by the real PRNG.  This is the weakest check
#      statistically but the only one that would notice if the KAT and our code
#      were wrong in the *same* way.

# Golden vectors are included by runtests.jl.

@testset "samplerz" begin

    @testset "constants match the reference" begin
        @test RCDT == RCDT_REF
        @test EXP_COEFFS == EXP_COEFFS_REF
        @test length(RCDT) == 18
        @test length(EXP_COEFFS) == 13
        @test Falcon.RCDT_PREC == 72
        # RCDT must be strictly decreasing -- it is a reverse CDF
        @test all(i -> RCDT[i] > RCDT[i + 1], 1:(length(RCDT) - 1))
        @test RCDT[end] == 1
        # and must fit in 72 bits but not 64: the reason for UInt128
        @test RCDT[1] < UInt128(1) << 72
        @test RCDT[1] > typemax(UInt64)
        # the top exp coefficient is exactly 2^63, which is why UInt64 and not Int64
        @test EXP_COEFFS[end] == UInt64(1) << 63
        @test EXP_COEFFS[end] > typemax(Int64)
    end

    @testset "LN2 and ILN2 are the reference's truncated decimals" begin
        # NOT the correctly-rounded values.  Getting "helpful" here changes the
        # byte stream and breaks the KAT; see the docstring in samplerz.jl.
        @test Falcon.LN2 == 0.69314718056
        @test Falcon.ILN2 == 1.44269504089
        @test Falcon.LN2 != log(2)
        @test Falcon.ILN2 != 1 / log(2)
        # ... but they are the right numbers to about 1e-12
        @test isapprox(Falcon.LN2, log(2); rtol = 1e-11)
        @test isapprox(Falcon.ILN2, 1 / log(2); rtol = 1e-11)
    end

    @testset "RCDT really is a half-Gaussian reverse CDF at sigma = 1.8205" begin
        # RCDT[i] / 2^72 should be P(z0 >= i) for z0 ~ D_{Z+, 0, 1.8205}: the
        # Gaussian restricted to the NON-NEGATIVE integers, each weighted once.
        #
        # Not the "folded" distribution that counts z > 0 twice -- that model
        # gives P(z0 = 0) = 0.219 against the table's 0.359, and was this
        # test's first, wrong, guess (docs/debug_log.md #023).  No folding
        # happens because samplerz supplies the sign separately.
        #
        # The comparison is ABSOLUTE, not relative.  The table holds integers
        # spanning 3e21 down to 1, so for the deep tail the quantisation to an
        # integer dominates any relative measure -- the second, wrong, version
        # of this test used rtol and failed on the last five entries against a
        # perfectly good table.  Measured agreement is within 8 units
        # everywhere, so a bound of 16 is sharp: it would catch any transposed
        # digit above the least significant couple of places.
        setprecision(BigFloat, 300) do
            s = BigFloat("1.8205")
            w(z) = exp(-BigFloat(z)^2 / (2 * s^2))
            total = sum(w(z) for z in 0:400)
            scale = BigFloat(2)^72
            for i in 1:length(RCDT)
                model = (sum(w(z) for z in i:400) / total) * scale
                @test abs(BigInt(RCDT[i]) - round(BigInt, model)) <= 16
            end
            # and the top entry is the one that pins sigma itself: P(z0 = 0)
            @test isapprox(1 - Float64(BigFloat(RCDT[1]) / scale), 0.3594977747;
                           rtol = 1e-9)
        end
    end

    @testset "approxexp approximates 2^64 * ccs * exp(-x)" begin
        for (x, ccs, want) in APPROXEXP_KAT
            @test approxexp(x, ccs) == want           # exact, fixed point
        end
        # and the fixed-point value really is the function it claims to be.
        # NOTE the scale: 2^64, not the 2^63 the reference's docstring claims.
        # See docs/debug_log.md #023.
        for x in 0.0:0.05:0.69, ccs in (0.5, 0.75, 0.9)
            got = Float64(approxexp(x, ccs)) / 2.0^64
            @test isapprox(got, ccs * exp(-x); rtol = 1e-10)
        end
        # the doubling is real and deliberate: check the ratio explicitly
        @test approxexp(0.0, 0.5) == UInt64(1) << 63
        # and it never overflows UInt64, which needs ccs < 1 strictly
        @test approxexp(0.0, 0.9999999) < typemax(UInt64)
        # never zero on the domain berexp uses -- berexp subtracts 1 from it,
        # and on UInt64 a zero would wrap to all-ones and change the comparison
        for x in 0.0:0.01:0.6931, ccs in (0.5, 0.7499908532676649, 0.999999)
            @test approxexp(x, ccs) > 0
        end
    end

    @testset "basesampler against the reference" begin
        for (bytes, want) in BASESAMPLER_KAT
            @test basesampler(ReplayBytes(bytes)) == want
        end
        # extremes: u = 0 falls below every entry (z0 = 18); u = all ones falls
        # below none (z0 = 0).
        @test basesampler(ReplayBytes(zeros(UInt8, 9))) == 18
        @test basesampler(ReplayBytes(fill(0xff, 9))) == 0
        # always consumes exactly 9 bytes, whatever it returns
        for (bytes, _) in BASESAMPLER_KAT
            rb = ReplayBytes(bytes)
            basesampler(rb)
            @test Falcon.consumed(rb) == 9
        end
        @test_throws ArgumentError basesampler(ReplayBytes(zeros(UInt8, 8)))
    end

    @testset "the official samplerz KAT" begin
        # 3072 vectors from the reference's samplerz_KAT512 / samplerz_KAT1024.
        # Replayed with reversed_chunks = true, which is the harness's
        # convention (docs/debug_log.md #022) -- feeding them in natural order
        # gives different, equally plausible answers.
        nfail = 0
        for (mu, sigma, sigmin, octets, want) in SAMPLERZ_KAT
            rb = ReplayBytes(octets; reversed_chunks = true)
            got = samplerz(mu, sigma, sigmin, rb)
            got == want || (nfail += 1)
        end
        @test nfail == 0
        @test length(SAMPLERZ_KAT) == 3072
    end

    @testset "the KAT vectors respect the documented parameter range" begin
        # samplerz requires 1 < sigmin < sigma < MAX_SIGMA; if the vectors did
        # not satisfy it, passing them would prove less than it appears to.
        for (mu, sigma, sigmin, _, _) in SAMPLERZ_KAT
            @test 1 < sigmin < sigma < 1.8205
        end
    end

    @testset "the replay convention is load-bearing" begin
        # Feeding the same bytes in natural order must give a *different*
        # answer on a decent fraction of the vectors -- otherwise the
        # reversed_chunks flag would be decoration and the KAT above would not
        # actually be pinning the byte order.
        differing = 0
        checked = 0
        for (mu, sigma, sigmin, octets, want) in SAMPLERZ_KAT[1:200]
            rb = ReplayBytes(octets; reversed_chunks = false)
            got = try
                samplerz(mu, sigma, sigmin, rb)
            catch e
                e isa ArgumentError && (differing += 1; checked += 1; continue)
                rethrow()
            end
            checked += 1
            got == want || (differing += 1)
        end
        @test checked == 200
        @test differing > 100
    end

    @testset "byte accounting" begin
        # The number of bytes consumed is part of the specification: 9 for the
        # base sampler, 1 for the sign, 1..8 in berexp, per rejection round.
        # So consumption must be 10 + k for some 1 <= k <= 8 per round.
        for (mu, sigma, sigmin, octets, _) in SAMPLERZ_KAT[1:500]
            rb = ReplayBytes(octets; reversed_chunks = true)
            samplerz(mu, sigma, sigmin, rb)
            used = Falcon.consumed(rb)
            rounds, rem = divrem(used, 11)      # not exact; just bound it
            @test used >= 11                    # at least one full round
            @test used <= length(octets)
        end
    end

    @testset "running off the end of the replay is an error, not a wrong answer" begin
        mu, sigma, sigmin, octets, _ = SAMPLERZ_KAT[1]
        @test_throws ArgumentError samplerz(mu, sigma, sigmin,
                                            ReplayBytes(octets[1:5]; reversed_chunks = true))
    end

    @testset "the distribution, driven by the real PRNG" begin
        # The statistical check.  Weak on its own, but it is the only test that
        # would notice if our code and the KAT were wrong in the same way --
        # e.g. if the constants were self-consistent but described the wrong
        # width.
        rng = chacha20(collect(UInt8, 0:55))
        src = bytesource(rng)
        sigma = 1.7
        sigmin = FALCON_512.sigma_min
        N = 200_000

        for mu in (0.0, 0.5, -3.25)
            samples = Vector{Int}(undef, N)
            for i in 1:N
                samples[i] = samplerz(mu, sigma, sigmin, src)
            end
            m = sum(samples) / N
            v = sum((samples .- m) .^ 2) / (N - 1)

            # The discrete Gaussian's mean is mu and its variance is sigma^2,
            # up to O(exp(-2 pi^2 sigma^2)) corrections that are far below the
            # sampling error here.  Standard error of the mean is sigma/sqrt(N).
            se = sigma / sqrt(N)
            @test abs(m - mu) < 5 * se
            # Standard error of the variance is roughly sigma^2*sqrt(2/N).
            sev = sigma^2 * sqrt(2 / N)
            @test abs(v - sigma^2) < 5 * sev

            # Chi-squared goodness of fit against the discrete Gaussian.
            lo, hi = floor(Int, mu) - 8, floor(Int, mu) + 8
            obs = zeros(Int, hi - lo + 1)
            outside = 0
            for z in samples
                if lo <= z <= hi
                    obs[z - lo + 1] += 1
                else
                    outside += 1
                end
            end
            wts = [exp(-(z - mu)^2 / (2 * sigma^2)) for z in lo:hi]
            tot = sum(exp(-(z - mu)^2 / (2 * sigma^2)) for z in (lo - 60):(hi + 60))
            expc = N .* (wts ./ tot)
            chi2 = sum((obs[i] - expc[i])^2 / expc[i] for i in eachindex(obs)
                       if expc[i] >= 10)
            dof = count(>=(10), expc) - 1
            # 17 cells, so dof is about 16; the 1e-6 upper tail is around 60.
            # A wrong width would blow past this by orders of magnitude.
            @test chi2 < 3 * dof + 30
            @test outside < N * 0.001
        end
    end

    @testset "samplerz is deterministic given its randomness" begin
        for (mu, sigma, sigmin, octets, want) in SAMPLERZ_KAT[1:50]
            a = samplerz(mu, sigma, sigmin, ReplayBytes(octets; reversed_chunks = true))
            b = samplerz(mu, sigma, sigmin, ReplayBytes(octets; reversed_chunks = true))
            @test a == b == want
        end
    end

    @testset "shifting the centre by an integer shifts the sample" begin
        # samplerz splits mu into floor(mu) + r and only the fractional part
        # reaches the arithmetic, so the same bytes at mu and mu+k must give
        # answers differing by exactly k.  This pins that split.
        for (mu, sigma, sigmin, octets, want) in SAMPLERZ_KAT[1:100]
            for k in (-7, 1, 13)
                got = samplerz(mu + k, sigma, sigmin,
                               ReplayBytes(octets; reversed_chunks = true))
                @test got == want + k
            end
        end
    end
end
