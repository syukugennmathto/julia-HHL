# test_spec.jl -- the specification's own tables, transcribed and checked.
#
# WHY THIS FILE EXISTS
# --------------------
# From #002 until 2026-08-21 this project had no route to the specification
# PDF.  Every constant was pinned instead by two mutually independent
# reference implementations (C and Python) that agree exactly, and every
# `spec_ref` field said "TODO", on the principle that an invented table number
# is worse than an absent one (docs/debug_log.md #002, #046).
#
# The PDF is now in hand: *Falcon: Fast-Fourier Lattice-based Compact
# Signatures over NTRU*, Specification **v1.2 -- 01/10/2020**, 67 pages.
#
# This file transcribes the tables and equations that the implementation
# depends on, and checks them.  It is deliberately a *transcription*, not a
# re-derivation: the point is that the literals in `src/` can be compared,
# digit by digit, against a source that is neither of the two implementations
# they were originally taken from.  Where the specification gives a formula,
# the formula is evaluated too, so that a mistyped digit in the table and a
# mistyped digit in `src/` would have to coincide to pass.
#
# Section, equation and table numbers below refer to that document and are the
# values now written into `FalconParams.spec_ref`.

@testset "the specification's tables" begin

    # -----------------------------------------------------------------------
    # Table 3.3 -- Falcon parameter sets (spec p.51)
    # -----------------------------------------------------------------------
    @testset "Table 3.3: parameter sets" begin
        # Transcribed from the PDF.  The specification prints the standard
        # deviations to 12 significant figures with digit grouping; the
        # implementations carry full Float64 literals, so the comparison is at
        # the precision the specification actually states, not tighter.
        #
        #                              Falcon-512        Falcon-1024
        #   Ring degree n              512               1024
        #   Modulus q                  12289
        #   Standard deviation sigma   165.736 617 183   168.388 571 447
        #   sigma_min                  1.277 833 697     1.298 280 334
        #   sigma_max                  1.8205
        #   Max. sig. square norm      34 034 726        70 265 242
        #   Public key bytelength      897               1 793
        #   Signature bytelength       666               1 280
        for (p, n, sigma, sigmin, beta2, pk, sig) in (
                (FALCON_512,  512,  165.736617183, 1.277833697, 34034726, 897,  666),
                (FALCON_1024, 1024, 168.388571447, 1.298280334, 70265242, 1793, 1280))
            @test p.n == n
            @test p.q == 12289
            # 12 significant figures is what the table prints
            @test p.sigma     ≈ sigma  rtol = 5e-12
            @test p.sigma_min ≈ sigmin rtol = 5e-10
            @test p.sigma_max == 1.8205
            @test p.sig_bound == beta2          # exact: an integer in the table
            @test p.pubkey_bytes == pk          # exact
            @test p.sig_bytes == sig            # exact
        end

        # The table gives one sigma_max for both parameter sets.
        @test FALCON_512.sigma_max == FALCON_1024.sigma_max

        # Table 3.3 does NOT list the private key bytelength ("Private key
        # size (not listed above) is about three times that of a signature").
        # So those two literals remain [C-ref]-only, and that is recorded
        # rather than papered over.
        @test FALCON_512.privkey_bytes == 1281
        @test FALCON_1024.privkey_bytes == 2305
    end

    # -----------------------------------------------------------------------
    # Section 2.6 -- how the parameters are derived (spec pp.17-19)
    # -----------------------------------------------------------------------
    @testset "Section 2.6: the derivations behind Table 3.3" begin
        q = Q

        # (2.10)  q = 12*1024 + 1 = 12289, "the smallest prime of the form
        #         k*2n + 1", chosen to make the NTT possible.
        @test q == 12 * 1024 + 1
        @test q == 3 * 2^12 + 1

        # (2.12)  each coefficient of f and g is sampled from D_{Z, sigma_fg}
        #         with sigma_fg = 1.17 * sqrt(q / 2n).
        for p in (FALCON_512, FALCON_1024)
            @test p.sigma_fg ≈ 1.17 * sqrt(q / (2 * p.n)) rtol = 1e-15
        end

        # (2.11)  ||B||_GS <= 1.17 * sqrt(q) -- the rejection threshold that
        #         `gram_schmidt_quality()` implements.
        @test gram_schmidt_quality() == 1.17

        # (2.13)  sigma = (1/pi) * sqrt(log(4n(1+1/eps))/2) * 1.17 * sqrt(q),
        #         with eps <= 1/sqrt(Qs * lambda), Qs = 2^64 (NIST), and
        #         lambda = 128 for Level I, 256 for Level V.
        #
        # This is the equation the *signature* width comes from, and it
        # reproduces the Table 3.3 literals to the last bits.
        spec_sigma(n, eps) = (1 / pi) * sqrt(log(4 * n * (1 + 1 / eps)) / 2) *
                             1.17 * sqrt(q)
        Qs = 2.0^64
        @test spec_sigma(512, 1 / sqrt(Qs * 128)) ≈ FALCON_512.sigma  rtol = 1e-15
        @test spec_sigma(1024, 1 / sqrt(Qs * 256)) ≈ FALCON_1024.sigma rtol = 1e-15

        # sigma_min uses a DIFFERENT epsilon, and the specification does not
        # print the one it used -- Table 3.3 gives the value only.  Inverting
        # the smoothing parameter of Z recovers eps = 1/sqrt(2^64 * n^3), and
        # that reproduces both literals to 3e-13 (docs/debug_log.md #046).
        #
        # Recording the discrepancy is the point: `sigma` and `sigma_min` are
        # NOT the same Gaussian's parameter at the same epsilon, and reading
        # (2.13) as if they were is a plausible and wrong simplification.
        for p in (FALCON_512, FALCON_1024)
            @test smoothing_eta(falcon_eps(p.n)) ≈ p.sigma_min rtol = 1e-12
            # ...and the (2.13) epsilon does NOT give sigma_min
            lambda = p.n == 512 ? 128.0 : 256.0
            @test !isapprox(smoothing_eta(1 / sqrt(Qs * lambda)), p.sigma_min;
                            rtol = 1e-3)
        end

        # (2.14)  beta = tau_sig * sigma * sqrt(2n) with tau_sig = 1.1, and the
        #         acceptance bound is floor(beta^2).  Exact integer agreement.
        for p in (FALCON_512, FALCON_1024)
            @test Int(floor((1.1 * p.sigma * sqrt(2 * p.n))^2)) == p.sig_bound
        end
    end

    # -----------------------------------------------------------------------
    # Table 3.1 -- the base sampler's distribution (spec p.41)
    # -----------------------------------------------------------------------
    @testset "Table 3.1: the BaseSampler distribution" begin
        # pdt[i], scaled by 2^72, transcribed from the PDF.  RCDT is then
        # *derived* here rather than transcribed a second time, so that the
        # check against `src/samplerz.jl` goes through the relation the
        # specification states, (3.33) and the definitions on p.41:
        #
        #     cdt[i]  = sum_{j <= i} pdt[j]
        #     RCDT[i] = sum_{j > i} pdt[j] = 2^72 - cdt[i]
        pdt = BigInt[
            1697680241746640300030, 1459943456642912959616,
             928488355018011056515,  436693944817054414619,
             151893140790369201013,   39071441848292237840,
               7432604049020375675,    1045641569992574730,
                108788995549429682,       8370422445201343,
                   476288472308334,         20042553305308,
                      623729532807,            14354889437,
                         244322621,                3075302,
                             28626,                    197,
                                 1,
        ]
        @test length(pdt) == 19                     # support {0, ..., 18}

        # The distribution must be a distribution: the scaled masses sum to
        # exactly 2^72.  (Table 3.1's last cdt entry is 4722366482869645213696,
        # which is 2^72 -- a nice check that the transcription of nineteen
        # 22-digit numbers has no typo anywhere.)
        @test sum(pdt) == big(2)^72
        @test sum(pdt) == big"4722366482869645213696"

        # cdt[0] and RCDT[0], spelled out in the table, as a second anchor.
        @test pdt[1] == big"1697680241746640300030"
        @test big(2)^72 - pdt[1] == big"3024686241123004913666"

        # Now the actual claim: src/samplerz.jl's RCDT is Table 3.1's RCDT.
        cum = big(0)
        derived = BigInt[]
        for i in 1:18                                # RCDT[0..17]
            cum += pdt[i]
            push!(derived, big(2)^72 - cum)
        end
        @test length(RCDT) == 18
        @test BigInt.(RCDT) == derived

        # RCDT[18] = 0 and is not stored: BaseSampler's loop runs i = 0..17
        # (algorithm 12, line 3), so an entry that can never be less than u
        # would only be dead weight.
        @test big(2)^72 - sum(pdt) == 0

        # (3.33): chi(i) = 2^-72 * pdt[i], and chi is "extremely close to the
        # half-Gaussian D_{Z+, sigma_max}".  Check that closeness directly --
        # this is what would catch a transposition that preserved the sum.
        let smax = 1.8205
            norm = sum(exp(-i^2 / (2 * smax^2)) for i in 0:200)
            for i in 0:18
                want = exp(-i^2 / (2 * smax^2)) / norm
                got = Float64(pdt[i + 1] / big(2)^72)
                # the tail entries are tiny; compare relatively where that is
                # meaningful and absolutely where it is not
                @test isapprox(got, want; rtol = 0.02, atol = 1e-18)
            end
        end

        # RCDT_PREC is the 72 of "scaled by a factor 2^72".
        @test Falcon.RCDT_PREC == 72
        # sigma_max is the width of that half-Gaussian, and Table 3.3's value.
        # `src/samplerz.jl` spells it inline in INV_2SIGMA2.
        @test Falcon.INV_2SIGMA2 == 1 / (2 * FALCON_512.sigma_max^2)
    end

    # -----------------------------------------------------------------------
    # Table 3.2 -- SamplerZ test vectors (spec pp.44-45)
    # -----------------------------------------------------------------------
    @testset "Table 3.2: SamplerZ test vectors" begin
        # Sixteen vectors, transcribed from the PDF.  These matter because
        # they are an *independent* source from the reference implementation's
        # own samplerz_KAT files (which test_samplerz.jl already replays, all
        # 3072 of them).  Two sources agreeing is the whole method of this
        # project; here the second source is the specification itself.
        #
        # Byte order: these replay with reversed_chunks = TRUE, the same
        # convention as the reference's KAT files.
        #
        # That answers a question #022 left open.  #022 found that the
        # official KAT vectors only match if each chunk is fed reversed, and
        # concluded "this convention is written down nowhere but in the
        # harness code".  It is written down here.  The specification's own
        # p.43 prose -- "at each iteration, the first 9 random bytes are used
        # by BaseSampler, the next one by line 5 and the last one(s) by
        # BerExp" -- fixes *which* bytes go to which step, not the order in
        # which the nine are assembled into the 72-bit word, and Table 3.2 is
        # written in the same order the KAT files are.  Checked both ways:
        # reversed_chunks = false matches 4 of 16 (docs/debug_log.md #046).
        #
        # Table 3.2 line 1 is also, byte for byte, the KAT vector that #022
        # failed on: mu = -91.90471153063714, expected -92, got -95.
        spec_vectors = [
            (-91.90471153063714,  1.7037990414754918, "0fc5442ff043d66e91d1eacac64ea5450a22941edc6c", -92),
            ( -8.322564895434937, 1.7037990414754918, "f4da0f8d8444d1a77265c2ef6f98bbbb4bee7db8d9b3",  -8),
            (-19.096516109216804, 1.7035823083824078, "db47f6d7fb9b19f25c36d6b9334d477a8bc0be68145d", -20),
            (-11.335543982423326, 1.7035823083824078,
             "ae41b4f5209665c74d00dcc1a8168a7bb516b3190cb42c1ded26cd52aed770eca7dd334e0547bcc3c163ce0b", -12),
            (  7.9386734193997555, 1.6984647769450156,
             "31054166c1012780c603ae9b833cec73f2f41ca5807cc89c92158834632f9b1555",   8),
            (-28.990850086867255, 1.6984647769450156, "737e9d68a50a06dbbc6477", -30),
            ( -9.071257914091655, 1.6980782114808988, "a98ddd14bf0bf22061d632", -10),
            (-43.88754568839566,  1.6980782114808988, "3cbf6818a68f7ab9991514", -41),
            (-58.17435547946095,  1.7010983419195522, "6f8633f5bfa5d26848668e3d5ddd46958e97630410587c", -61),
            (-43.58664906684732,  1.7010983419195522, "272bc6c25f5c5ee53f83c43a361fbc7cc91dc783e20a", -46),
            (-34.70565203313315,  1.7009387219711465, "45443c59574c2c3b07e2e1d9071e6d133dbe32754b0a", -34),
            (-44.36009577368896,  1.7009387219711465,
             "6ac116ed60c258e2cbaeab728c4823e6da36e18d08da5d0cc104e21cc7fd1f5ca8d9dbb675266c928448059e", -44),
            (-21.783037079346236, 1.6958406126012802, "68163bc1e2cbf3e18e7426", -23),
            (-39.68827784633828,  1.6958406126012802, "d6a1b51d76222a705a0259", -40),
            (-18.488607061056847, 1.6955259305261838, "f0523bfaa8a394bf4ea5c10f842366fde286d6a30803", -22),
            (-48.39610939101591,  1.6955259305261838, "87bd87e63374cee62127fc6931104aab64f136a0485b", -50),
        ]
        @test length(spec_vectors) == 16

        # "Table 3.2: Test vectors for SamplerZ (sigma_min = 1.277 833 697)",
        # i.e. the Falcon-512 sigma_min.
        sigmin = FALCON_512.sigma_min

        hexbytes(s) = UInt8[parse(UInt8, s[i:i+1]; base = 16) for i in 1:2:length(s)]

        for (mu, sigma, hex, want) in spec_vectors
            octets = hexbytes(hex)
            rb = ReplayBytes(octets; reversed_chunks = true)   # see above
            @test samplerz(mu, sigma, sigmin, rb) == want
            # Every byte in the vector is consumed: the specification split
            # the strings by iteration, so a sampler that stopped early or ran
            # long would be taking a different path even if the answer matched.
            @test Falcon.consumed(rb) == length(octets)
        end

        # The documented parameter range holds for all of them (1 < sigma_min
        # < sigma' < sigma_max), which is what makes passing them meaningful.
        for (_, sigma, _, _) in spec_vectors
            @test 1 < sigmin < sigma < FALCON_512.sigma_max
        end

        # Lines 4 and 12 are the long ones -- the specification calls out that
        # BerExp occasionally needs more than one byte (p.43), and line 9 is
        # its worked example of that.  Assert the vectors still exercise the
        # multi-iteration path, so that a future edit cannot quietly replace
        # them with sixteen single-shot cases.  One iteration is 11 bytes
        # (9 for BaseSampler, 1 for the sign bit, 1 for BerExp), so the
        # lengths present -- 11, 22, 23, 33, 44 -- are 1, 2, 3 and 4
        # iterations, with the 23 being line 9's second BerExp byte.
        @test sort(unique(length(hexbytes(v[3])) for v in spec_vectors)) ==
              [11, 22, 23, 33, 44]

        # And the convention is load-bearing here too, not decoration: fed in
        # natural order the same vectors mostly give different, equally
        # plausible integers.  Without this, the `reversed_chunks = true`
        # above could be silently wrong and nothing would say so.
        let wrong = count(spec_vectors) do (mu, sigma, hex, want)
                rb = ReplayBytes(hexbytes(hex))
                got = try
                    samplerz(mu, sigma, sigmin, rb)
                catch e
                    e isa ArgumentError || rethrow()
                    nothing                     # ran off the end of the string
                end
                got != want
            end
            @test wrong == 12
        end
    end
end
