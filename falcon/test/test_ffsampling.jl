# test_ffsampling.jl -- module 8.
#
# The tree is the hardest object in the project to inspect: at n = 512 it has
# 512 leaves and 511 internal nodes of complex vectors, so "diff it against the
# C implementation by eye" is not a plan.  Hence the approach settled on when
# the project was planned: get the toy degrees right first (n = 8, 16, 32,
# where the whole tree fits on a screen), and use *properties* at n = 512.
#
# The properties available are unusually good here:
#
#   * the Gram matrix must be Hermitian and positive definite -- if it is not,
#     the adjoint's sign is wrong (module 3) and LDL* is meaningless;
#   * L * D * adj(L)^T must reconstruct G, which is the factorisation's
#     defining equation and needs no reference at all;
#   * the leaf sigmas must land in [sigma_min, sigma_max], which is exactly
#     samplerz's precondition and exactly what key generation's Gram-Schmidt
#     bound was chosen to guarantee.  This one ties three modules together.
#
# And the byte-exact check: ffSampling driven by a fixed stream must consume
# the same bytes and return the same integers as the reference.

# Golden vectors are included by runtests.jl.

# Rebuild a key and its normalised tree from a recorded label.
function _tree_from_label(n, label, used, sigma)
    rb = ReplayBytes(shake256(codeunits(label), used))
    f, g, F, G = ntru_gen(n, rb)
    B = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    B[1, 1] = fft(Float64.(g)); B[1, 2] = fft(Float64.(-f))
    B[2, 1] = fft(Float64.(G)); B[2, 2] = fft(Float64.(-F))
    return (falcon_tree(B, sigma), f, g, F, G, B)
end

@testset "ffsampling" begin

    @testset "the Gram matrix is Hermitian and positive definite" begin
        for (n, label, used, sigma, _, _) in FFLDL_KAT
            _, _, _, _, _, B = _tree_from_label(n, label, used, sigma)
            G = gram_fft(B)
            # diagonal entries real and strictly positive
            for i in 1:2
                @test maximum(abs, imag.(G[i, i])) < 1e-8 * maximum(abs, real.(G[i, i]))
                @test all(>(0), real.(G[i, i]))
            end
            # G[2,1] = conj(G[1,2])
            @test _maxabsdiff(G[2, 1], conj.(G[1, 2])) <
                  1e-9 * max(1.0, maximum(abs, G[1, 2]))
        end
    end

    @testset "L * D * adj(L)^T reconstructs G" begin
        # The defining equation of the factorisation.  No oracle needed.
        for (n, label, used, sigma, _, _) in FFLDL_KAT
            _, _, _, _, _, B = _tree_from_label(n, label, used, sigma)
            G = gram_fft(B)
            L10, D00, D11 = ldl_fft(G)
            one_ = ones(ComplexF64, n)
            # (L D L*)[1,1] = D00
            @test _maxabsdiff(D00, G[1, 1]) < 1e-9 * maximum(abs, G[1, 1])
            # (L D L*)[2,1] = L10 * D00
            @test _maxabsdiff(mul_fft(L10, D00), G[2, 1]) <
                  1e-8 * max(1.0, maximum(abs, G[2, 1]))
            # (L D L*)[2,2] = L10 D00 adj(L10) + D11
            rec = add_fft(mul_fft(mul_fft(L10, adj_fft(L10)), D00), D11)
            @test _maxabsdiff(rec, G[2, 2]) < 1e-8 * maximum(abs, G[2, 2])
            # D11 must be real and positive: it is a squared Gram-Schmidt norm
            @test maximum(abs, imag.(D11)) < 1e-8 * maximum(abs, real.(D11))
            @test all(>(0), real.(D11))
        end
    end

    @testset "tree shape" begin
        for (n, label, used, sigma, want_sigmas, want_l10) in FFLDL_KAT
            T, _, _, _, _, _ = _tree_from_label(n, label, used, sigma)
            @test nleaves(T) == n
            @test treedepth(T) == Int(log2(n))
            @test length(leaf_sigmas(T)) == n
            @test length(node_l10s(T)) == n - 1        # a binary tree with n leaves
            @test length(want_sigmas) == n
            @test length(want_l10) == n - 1
        end
    end

    @testset "ffLDL tree against the reference" begin
        for (n, label, used, sigma, want_sigmas, want_l10) in FFLDL_KAT
            T, _, _, _, _, _ = _tree_from_label(n, label, used, sigma)
            got = leaf_sigmas(T)
            @test _maxabsdiff(got, want_sigmas) < 1e-9 * maximum(abs, want_sigmas)
            gl = node_l10s(T)
            @test length(gl) == length(want_l10)
            for (a, b) in zip(gl, want_l10)
                @test _maxabsdiff(a, b) < 1e-9 * max(1.0, maximum(abs, b))
            end
        end
    end

    @testset "leaf sigmas satisfy samplerz's precondition" begin
        # This is the property that ties key generation, the parameters and the
        # sampler together: keygen rejects (f,g) unless the Gram-Schmidt norm is
        # at most 1.17^2 q, and *that* is what makes every leaf sigma land in
        # [sigma_min, sigma_max] so samplerz is called in range.
        p = FALCON_512
        for (n, label, used, sigma, _, _) in FFLDL_KAT
            T, _, _, _, _, _ = _tree_from_label(n, label, used, sigma)
            for s in leaf_sigmas(T)
                @test isfinite(s)
                @test p.sigma_min < s < p.sigma_max
            end
        end
    end

    @testset "an un-normalised tree is refused, not silently sampled" begin
        n, label, used, sigma, _, _ = FFLDL_KAT[1]
        rb = ReplayBytes(shake256(codeunits(label), used))
        f, g, F, G = ntru_gen(n, rb)
        B = Matrix{Vector{ComplexF64}}(undef, 2, 2)
        B[1, 1] = fft(Float64.(g)); B[1, 2] = fft(Float64.(-f))
        B[2, 1] = fft(Float64.(G)); B[2, 2] = fft(Float64.(-F))
        raw = ffldl_fft(gram_fft(B))            # NOT normalised
        t = (fft(zeros(n)), fft(zeros(n)))
        src = ReplayBytes(shake256(codeunits("x"), 10_000))
        @test_throws ArgumentError ffsampling_fft(t, raw, FALCON_512.sigma_min, src)
    end

    @testset "our roots and the reference's differ, and it shows here" begin
        # THE finding of this module.  Our FFT roots are correctly rounded and
        # agree with the C reference; the Python reference's table carries only
        # ~15 significant digits (docs/debug_log.md #010).  The difference is
        # 1e-16 relative and invisible everywhere until ffSampling, where the
        # FFT result is the *centre* of a discrete Gaussian and a last-ulp
        # change flips a rejection comparison.
        #
        # Measured: with our roots, 1 of 6 recorded vectors reproduces; with
        # the reference table installed, 6 of 6.  Both numbers are asserted,
        # because the point is not that one is right -- it is that a 1-ulp
        # table difference changes the signature.  docs/debug_log.md #025.
        function nmatching()
            k = 0
            for (n, label, kused, t0, t1, slabel, sused, wz0, wz1) in FFSAMPLING_KAT
                T, _, _, _, _, _ = _tree_from_label(n, label, kused, 165.7366171829776)
                src = ReplayBytes(shake256(codeunits(slabel), sused + 20_000))
                z0, z1 = ffsampling_fft((fft(t0), fft(t1)), T,
                                        FALCON_512.sigma_min, src)
                (round.(Int, ifft(z0)) == wz0 && round.(Int, ifft(z1)) == wz1) && (k += 1)
            end
            return k
        end
        @test nmatching() < length(FFSAMPLING_KAT)          # our roots: not all
        @test with_fft_roots(nmatching, ROOTS_C) == length(FFSAMPLING_KAT)
        # and the tables really are different, but only just
        reset_fft_roots!()
        for n in (8, 512)
            @test Falcon.fft_roots(n) != ROOTS_C[n]
            @test _maxabsdiff(Falcon.fft_roots(n), ROOTS_C[n]) < 1e-13
        end
    end

    @testset "ffnp against the reference (deterministic)" begin
        for (n, label, used, t0, t1, wz0, wz1) in FFNP_KAT
            T, _, _, _, _, _ = _tree_from_label(n, label, used, 165.7366171829776)
            z0, z1 = ffnp_fft((fft(t0), fft(t1)), T)
            @test round.(Int, ifft(z0)) == wz0
            @test round.(Int, ifft(z1)) == wz1
        end
    end

    @testset "ffSampling against the reference, byte for byte" begin
        # Replayed with the *reference's* root table installed, because that is
        # what the recorded vectors were produced with.  See the testset above:
        # with our own (more accurate) roots these do not reproduce, and that
        # is the finding rather than a defect.
        with_fft_roots(ROOTS_C) do
        for (n, label, kused, t0, t1, slabel, sused, wz0, wz1) in FFSAMPLING_KAT
            T, _, _, _, _, _ = _tree_from_label(n, label, kused, 165.7366171829776)
            src = ReplayBytes(shake256(codeunits(slabel), sused))
            z0, z1 = ffsampling_fft((fft(t0), fft(t1)), T,
                                    FALCON_512.sigma_min, src)
            @test round.(Int, ifft(z0)) == wz0
            @test round.(Int, ifft(z1)) == wz1
            @test Falcon.consumed(src) == sused     # exact byte accounting
        end
        end
    end

    @testset "the sampled point gives a short signature" begin
        # The property that actually matters, and the one GPV guarantees:
        # (t - z) * B must be SHORT.  The residual (t - z) itself is not small
        # coefficient-wise -- measured around 32 at n = 512 -- because the
        # sampling happens in the Gram-Schmidt basis, not coordinate-wise.  An
        # earlier version of this test asserted |t - z| < 8 and was simply
        # wrong about what ffSampling promises (docs/debug_log.md #026).
        #
        # The expectation is ||s||^2 ~ 2n * sigma^2, since there are 2n
        # coordinates each of variance sigma^2.  A factor of two of headroom
        # catches a sampler that is systematically too wide without being
        # flaky.
        sigma = 165.7366171829776
        for (n, label, kused, t0, t1, slabel, _, _, _) in FFSAMPLING_KAT
            T, _, _, _, _, B = _tree_from_label(n, label, kused, sigma)
            # generous byte budget: with our roots the rejection path differs
            # from the reference's, so the recorded byte count does not apply
            src = ReplayBytes(shake256(codeunits(slabel), 200_000))
            tf = (fft(t0), fft(t1))
            z0, z1 = ffsampling_fft(tf, T, FALCON_512.sigma_min, src)

            # the sampled values are integers
            @test _maxabsdiff(ifft(z0), round.(ifft(z0))) < 1e-6
            @test _maxabsdiff(ifft(z1), round.(ifft(z1))) < 1e-6

            d0 = sub_fft(tf[1], z0)
            d1 = sub_fft(tf[2], z1)
            s0 = ifft(add_fft(mul_fft(d0, B[1, 1]), mul_fft(d1, B[2, 1])))
            s1 = ifft(add_fft(mul_fft(d0, B[1, 2]), mul_fft(d1, B[2, 2])))
            nrm = sum(abs2, s0) + sum(abs2, s1)
            @test 0 < nrm < 2 * (2n) * sigma^2
        end
    end

    @testset "the target dimension, by property" begin
        # n = 512 with a real key: no reference, only the properties above.
        # This is the dimension the tree is actually used at, and the one where
        # eyeballing is impossible.
        rng = chacha20(collect(UInt8, 0:55))
        src = bytesource(rng)
        f, g, F, G = ntru_gen(512, src)
        B = Matrix{Vector{ComplexF64}}(undef, 2, 2)
        B[1, 1] = fft(Float64.(g)); B[1, 2] = fft(Float64.(-f))
        B[2, 1] = fft(Float64.(G)); B[2, 2] = fft(Float64.(-F))

        Gm = gram_fft(B)
        @test all(>(0), real.(Gm[1, 1]))
        @test _maxabsdiff(Gm[2, 1], conj.(Gm[1, 2])) < 1e-6 * maximum(abs, Gm[1, 2])

        T = falcon_tree(B, FALCON_512.sigma)
        @test nleaves(T) == 512
        @test treedepth(T) == 9
        p = FALCON_512
        for s in leaf_sigmas(T)
            @test p.sigma_min < s < p.sigma_max
        end

        # And the property that ties this module to the signature bound: the
        # sampled lattice point, expressed through the basis, must be short
        # enough that verification would accept it.  Measured ||s||^2 is about
        # 2.8e7 against beta^2 = 3.4e7, which is exactly the expected
        # 2n*sigma^2 = 2.81e7 -- the sampler is running at the width the
        # parameters were designed for.
        for _ in 1:3
            t = (fft(fill(0.5, 512)), fft(fill(-0.5, 512)))
            z0, z1 = ffsampling_fft(t, T, p.sigma_min, src)
            @test _maxabsdiff(ifft(z0), round.(ifft(z0))) < 1e-6
            d0 = sub_fft(t[1], z0); d1 = sub_fft(t[2], z1)
            s0 = ifft(add_fft(mul_fft(d0, B[1, 1]), mul_fft(d1, B[2, 1])))
            s1 = ifft(add_fft(mul_fft(d0, B[1, 2]), mul_fft(d1, B[2, 2])))
            nrm = sum(abs2, s0) + sum(abs2, s1)
            @test nrm <= p.sig_bound            # verification would accept
            @test nrm > 0.5 * (2 * 512) * p.sigma^2   # and it is not degenerate
        end
    end
end
