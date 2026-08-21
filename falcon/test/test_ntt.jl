# test_ntt.jl -- module 4.
#
# Three independent lines of attack, because an NTT that is subtly wrong tends
# to be self-consistent (it round-trips fine while computing the wrong
# transform, and then only the *product* is wrong):
#
#   1. against the Python reference's vectors (test/vectors/ntt_kat.jl);
#   2. against the definition -- ntt(f)[j] must equal f(w_j) by Horner;
#   3. against poly.jl -- the NTT product must equal the schoolbook product.
#
# (2) is the one that would catch a transform that round-trips but evaluates at
# the wrong points; (3) is the one that would catch a wrong root *ordering*,
# since a permuted-but-consistent table still round-trips.

# Golden vectors are included by runtests.jl.

@testset "ntt" begin

    @testset "the root of unity" begin
        @test Falcon.NTT_ZETA == 7
        @test Falcon.find_ntt_zeta() == Falcon.NTT_ZETA
        # order exactly 2048
        @test powermod(Falcon.NTT_ZETA, 2048, Q) == 1
        @test powermod(Falcon.NTT_ZETA, 1024, Q) != 1
        # zeta^1024 = -1, which is what makes the odd powers roots of x^n+1
        @test powermod(Falcon.NTT_ZETA, 1024, Q) == Q - 1
        # 2n must divide q-1 for every degree we support
        @test (Q - 1) % 2048 == 0
        @test_throws ArgumentError Falcon.find_ntt_zeta(Q, 4096 * 4)
    end

    @testset "root tables match the reference" begin
        for (n, want) in ROOTS_ZQ
            @test Falcon.ntt_roots(n) == want
        end
        # ... and are genuinely the roots of x^n + 1, in +- pairs, ordered
        # with the smaller residue first.
        for n in sort(collect(keys(ROOTS_ZQ)))
            w = Falcon.ntt_roots(n)
            @test length(w) == n
            @test length(unique(w)) == n
            @test all(x -> powermod(x, n, Q) == Q - 1, w)      # w^n = -1
            for i in 1:2:n
                @test w[i] + w[i + 1] == Q                     # the pair is (r, q-r)
                @test w[i] < w[i + 1]                          # smaller first
            end
            # the parent relation that defines the tower
            if n > 2
                parent = Falcon.ntt_roots(n ÷ 2)
                for i in 1:(n ÷ 2)
                    @test powermod(w[2i - 1], 2, Q) == parent[i]
                end
            end
        end
    end

    @testset "degree checking" begin
        @test_throws ArgumentError Falcon.ntt_roots(1)
        @test_throws ArgumentError Falcon.ntt_roots(3)
        @test_throws ArgumentError Falcon.ntt_roots(2048)
        @test_throws ArgumentError ntt([1])
        @test_throws ArgumentError ntt([1, 2, 3])
    end

    @testset "INV2_Q" begin
        @test mod(2 * Falcon.INV2_Q, Q) == 1
        @test Falcon.INV2_Q == 6145
    end

    @testset "invmod_q" begin
        for x in (1, 2, 3, 7, 1479, Q - 1, 12288)
            @test mod(x * Falcon.invmod_q(x), Q) == 1
        end
        @test_throws DivideError Falcon.invmod_q(0)
        @test_throws DivideError Falcon.invmod_q(Q)      # 0 mod q
    end

    @testset "forward NTT against the reference" begin
        for (f, want) in NTT_KAT
            @test ntt(f) == want
        end
    end

    @testset "NTT is evaluation at the roots" begin
        # The definition, checked by Horner.  If this passes, the transform is
        # correct regardless of how the butterflies are arranged.
        for (f, want) in NTT_KAT
            n = length(f)
            w = Falcon.ntt_roots(n)
            @test want == [_evalpoly_mod(f, w[j], Q) for j in 1:n]
        end
    end

    @testset "round trip" begin
        for (f, F) in NTT_KAT
            @test intt(F) == f
            @test intt(ntt(f)) == f
            @test ntt(intt(F)) == F
        end
    end

    @testset "split_ntt / merge_ntt are mutually inverse" begin
        for (_, F) in NTT_KAT
            length(F) < 4 && continue
            f0, f1 = Falcon.split_ntt(F)
            @test Falcon.merge_ntt(f0, f1) == F
        end
        # and they correspond to polysplit/polymerge in the coefficient domain
        for (f, F) in NTT_KAT
            length(f) < 4 && continue
            f0, f1 = polysplit(f)
            g0, g1 = Falcon.split_ntt(F)
            @test ntt(f0) == g0
            @test ntt(f1) == g1
        end
    end

    @testset "pointwise operations" begin
        for (fn, gn, wmul, wdiv) in NTT_POINTWISE
            @test ntt_mul(fn, gn) == wmul
            @test ntt_div(fn, gn) == wdiv
            @test ntt_mul(wdiv, gn) == fn         # (f/g)*g == f
            @test ntt_add(fn, gn) == Int[mod(fn[i] + gn[i], Q) for i in eachindex(fn)]
            @test ntt_sub(ntt_add(fn, gn), gn) == fn
        end
    end

    @testset "NTT product agrees with the schoolbook product" begin
        # The cross-check that matters: poly.jl computes the negacyclic
        # convolution by definition, this module computes it by evaluation and
        # interpolation.  Agreement means both got x^n = -1 right, and it also
        # pins the root *ordering* -- a permuted table would still round-trip
        # but would give a different product here.
        for (f, g, _, _, wmul) in NEGACYCLIC_ZQ
            @test polymulq_ntt(f, g) == wmul
            @test polymulq_ntt(f, g) == polymulq(f, g)
        end
    end

    @testset "division in R_q" begin
        for (f, g, want) in DIV_ZQ
            @test polydivq(f, g) == want
            @test polymulq(want, g) == f          # closes the loop via poly.jl
        end
    end

    @testset "invertibility" begin
        for (f, zero_at) in NOT_INVERTIBLE
            F = ntt(f)
            @test F[zero_at + 1] == 0             # vectors record a 0-based index
            @test !is_invertible_zq(f)
            @test_throws DivideError polydivq(f, f)
            g = [1; zeros(Int, length(f) - 1)]
            @test_throws DivideError polydivq(g, f)
        end
        # 1 is invertible; so is a random polynomial, essentially always
        for n in (8, 512)
            one_ = [1; zeros(Int, n - 1)]
            @test is_invertible_zq(one_)
            @test polydivq(one_, one_) == one_
        end
        for (f, _, _, _, _) in NEGACYCLIC_ZQ
            if is_invertible_zq(f)
                inv_f = polydivq([1; zeros(Int, length(f) - 1)], f)
                @test polymulq(inv_f, f) == [1; zeros(Int, length(f) - 1)]
            end
        end
    end

    @testset "linearity" begin
        for (f, F) in NTT_KAT
            n = length(f)
            g = circshift(f, 1)
            G = ntt(g)
            @test ntt(polyaddq(f, g)) == ntt_add(F, G)
            @test ntt(polysubq(f, g)) == ntt_sub(F, G)
        end
    end

    @testset "the iterative in-place NTT agrees with everything else" begin
        # Regression test for docs/debug_log.md #038.  ntt.jl now carries two
        # transforms: the recursive one above, which follows the specification's
        # tree ordering, and an iterative Cooley-Tukey/Gentleman-Sande pair used
        # only to multiply.  The second leaves its intermediate values in
        # bit-reversed order, which is *not* the specification's order -- that is
        # sound precisely because multiplication never looks at an individual
        # coordinate, and this testset is what holds that claim in place.
        rng = MersenneTwister(20260821)
        for n in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
            for _ in 1:3
                f = rand(rng, 0:(Q - 1), n)
                g = rand(rng, 0:(Q - 1), n)
                # against the schoolbook convolution, which shares no code
                @test polymulq_fast(f, g) == polymulq(f, g)
                # and against the recursive NTT, which shares no code either
                @test polymulq_fast(f, g) == polymulq_ntt(f, g)
            end
            # inputs that are not already reduced: signatures carry centred
            # coefficients, so negative values reach this path in verification
            f = rand(rng, -30000:30000, n)
            g = rand(rng, 0:(Q - 1), n)
            @test polymulq_fast(f, g) == polymulq(f, g)

            # the transforms are mutually inverse in their own ordering
            a = UInt32.(rand(rng, 0:(Q - 1), n))
            b = copy(a)
            Falcon.ntt_ip!(b, ntt_zetas(n))
            @test b != a || n == 1                    # it did something
            Falcon.intt_ip!(b, ntt_zetas(n))
            @test b == a
        end

        # The root table really is negacyclic: psi^n = -1, not +1.  Getting
        # this wrong gives a *cyclic* convolution, which is a different ring
        # and would show up only as wrong signatures.
        for n in (2, 8, 512, 1024)
            psi = Falcon._psi(n)
            @test powermod(psi, n, Q) == Q - 1
            @test powermod(psi, 2n, Q) == 1
            @test length(ntt_zetas(n)) == n
        end

        # in-place really is in place: no allocation once the tables are warm
        let n = 512, a = UInt32.(rand(rng, 0:(Q - 1), n)), z = ntt_zetas(n)
            Falcon.ntt_ip!(copy(a), z)                # warm up
            b = copy(a)
            @test (@allocated Falcon.ntt_ip!(b, z)) == 0
        end
    end

    @testset "the inverse twiddle table folds what it claims to fold" begin
        # Regression test for docs/debug_log.md #039.  intt_zetas precomputes
        # q - zeta *and* multiplies the single entry the last level uses by 1/n,
        # so a mistake here scales the whole transform by a constant and would
        # show up only as a wrong signature.  The relationship is asserted
        # directly rather than via a round trip, which would hide a
        # compensating error.
        for n in (8, 16, 64, 512, 1024)
            z = ntt_zetas(n)
            iz = intt_zetas(n)
            @test length(iz) == n
            ninv = Falcon.ntt_ninv(n)
            @test mod(Int(ninv) * n, Q) == 1
            for k in 1:n
                plain = mod(Q - Int(z[k]), Q)
                if k == 2 && n > 8
                    # the folded entry: the last level's only group
                    @test Int(iz[k]) == mod(plain * Int(ninv), Q)
                else
                    @test Int(iz[k]) == plain
                end
            end
        end
        # n = 8: the fused block is the whole transform, so nothing is folded
        # and every coefficient is scaled in the final pass instead.
        let z = ntt_zetas(8), iz = intt_zetas(8)
            @test all(k -> Int(iz[k]) == mod(Q - Int(z[k]), Q), 1:8)
        end
    end
end
