# test_poly.jl -- module 3.
#
# Two kinds of check:
#
#   * against the Python reference (test/vectors/poly_kat.jl).  For the mod-q
#     product this is meaningful in a strong sense: the reference computes it
#     with an NTT, we compute it by schoolbook convolution, so the two agree
#     only if both are right about the negacyclic sign.
#   * against mathematical properties that hold for *any* correct
#     implementation (ring axioms, the adjoint being an involution, split and
#     merge being mutually inverse).  These catch the cases the fixed vectors
#     happen to miss.

# Golden vectors are included by runtests.jl.

@testset "poly" begin

    @test POLY_Q == Q       # the vectors were generated at the same modulus

    @testset "negacyclic product over Z (reference vectors)" begin
        for (f, g, want) in NEGACYCLIC_Z
            @test polymul(f, g) == want
        end
    end

    @testset "negacyclic product over Z (small cases by hand)" begin
        # n = 2, R = Z[x]/(x^2+1) = the Gaussian integers.
        # (a + b x)(c + d x) = (ac - bd) + (ad + bc) x
        for (a, b, c, d) in ((1, 0, 0, 1), (0, 1, 0, 1), (2, 3, 4, 5), (-1, 7, 3, -2))
            @test polymul([a, b], [c, d]) == [a * c - b * d, a * d + b * c]
        end
        # x^n = -1 in the most direct possible form: x * x^{n-1} = -1.
        for n in (2, 4, 8, 16)
            x = zeros(Int, n); x[2] = 1
            xtop = zeros(Int, n); xtop[n] = 1
            want = zeros(Int, n); want[1] = -1
            @test polymul(x, xtop) == want
        end
    end

    @testset "ring axioms over Z" begin
        for (f, g, _) in NEGACYCLIC_Z
            n = length(f)
            one_ = zeros(BigInt, n); one_[1] = 1
            zero_ = zeros(BigInt, n)
            @test polymul(f, one_) == f
            @test polymul(f, zero_) == zero_
            @test polymul(f, g) == polymul(g, f)                    # commutative
            @test polyadd(f, g) == polyadd(g, f)
            @test polysub(f, g) == polyadd(f, polyneg(g))
            @test polyadd(f, polyneg(f)) == zero_
            # distributivity
            @test polymul(f, polyadd(g, one_)) == polyadd(polymul(f, g), f)
        end
        # associativity, on the small cases only (it is O(n^2) three times)
        for (f, g, _) in NEGACYCLIC_Z[1:3]
            h = circshift(g, 1)
            @test polymul(polymul(f, g), h) == polymul(f, polymul(g, h))
        end
    end

    @testset "overflow is detected, not silently wrapped" begin
        # The last reference vector was generated with coefficients of size
        # ~10^9 at n = 64 precisely so that the exact product does not fit in
        # Int64.  BigInt must give the right answer; Int64 must complain
        # loudly rather than return a wrapped one.
        f, g, want = NEGACYCLIC_Z[end]
        @test length(f) == 64
        @test maximum(abs, want) > typemax(Int64)
        @test polymul(f, g) == want                       # BigInt path

        fi = Int64.(f); gi = Int64.(g)
        @test_throws OverflowError polymul(fi, gi)        # Int64 path

        # ... and a case that does fit must go through unharmed.
        small = Int64[1, 2, 3, 4]
        @test polymul(small, small) == Int64.(polymul(BigInt.(small), BigInt.(small)))
    end

    @testset "adjoint" begin
        for (f, want) in ADJ_Z
            @test polyadj(f) == want
        end
        for (f, g, _) in NEGACYCLIC_Z
            # adj is an involution ...
            @test polyadj(polyadj(f)) == f
            # ... a ring anti-homomorphism (here: homomorphism, R is commutative) ...
            @test polyadj(polymul(f, g)) == polymul(polyadj(f), polyadj(g))
            @test polyadj(polyadd(f, g)) == polyadd(polyadj(f), polyadj(g))
            # ... and f*adj(f) is self-adjoint, i.e. its own adjoint.
            ffa = polymul(f, polyadj(f))
            @test polyadj(ffa) == ffa
            # The constant coefficient of f*adj(f) is ||f||^2: this is the
            # identity that makes the Gram matrix work.
            @test ffa[1] == sqnorm(f)
        end
    end

    @testset "sqnorm" begin
        @test sqnorm(Int[]) == 0
        @test sqnorm([3, 4]) == 25
        @test sqnorm([3], [4]) == 25
        @test sqnorm([1, 2], [3, 4]) == 30
        # exact for values that overflow Float64's integer range
        big = BigInt(2)^40
        @test sqnorm([big, big]) == 2 * big^2
        @test sqnorm([big, big]) isa BigInt
    end

    @testset "arithmetic mod q (reference vectors)" begin
        for (f, g, wadd, wsub, wmul) in NEGACYCLIC_ZQ
            @test polyaddq(f, g) == wadd
            @test polysubq(f, g) == wsub
            @test polymulq(f, g) == wmul
            # representatives must be reduced
            @test all(x -> 0 <= x < Q, wmul)
        end
    end

    @testset "mod q agrees with the exact product reduced" begin
        for (f, g, _) in NEGACYCLIC_Z
            exact = polymul(f, g)
            fq = Int[mod(c, Q) for c in f]
            gq = Int[mod(c, Q) for c in g]
            @test polymulq(fq, gq) == Int[mod(c, Q) for c in exact]
        end
    end

    @testset "centred representatives" begin
        @test centered([0, 1, Q - 1, Q ÷ 2, Q ÷ 2 + 1]) ==
              [0, 1, -1, Q ÷ 2, Q ÷ 2 + 1 - Q]
        # every centred value lies in (-q/2, q/2]
        for (f, _, _, _, _) in NEGACYCLIC_ZQ
            c = centered(f)
            @test all(x -> -(Q ÷ 2) <= x <= Q ÷ 2, c)
            @test Int[mod(x, Q) for x in c] == f
        end
    end

    @testset "split / merge" begin
        for (f, _, _) in NEGACYCLIC_Z
            f0, f1 = polysplit(f)
            @test length(f0) == length(f1) == length(f) ÷ 2
            @test polymerge(f0, f1) == f
            # f0 holds the even-index coefficients, f1 the odd ones
            @test f0 == f[1:2:end]
            @test f1 == f[2:2:end]
        end
        @test_throws ArgumentError polysplit([1, 2, 3])
        @test_throws DimensionMismatch polymerge([1, 2], [1])

        # The defining identity: f(x) = f0(x^2) + x * f1(x^2).  Check it by
        # evaluating both sides at a random point in C, which needs no
        # machinery beyond Horner.
        for (f, _, _) in NEGACYCLIC_Z
            f0, f1 = polysplit(f)
            z = 0.3 + 0.7im
            lhs = _evalpoly_c(f, z)
            rhs = _evalpoly_c(f0, z^2) + z * _evalpoly_c(f1, z^2)
            @test isapprox(lhs, rhs; rtol = 1e-9, atol = 1e-9 * max(1, abs(lhs)))
        end
    end

    @testset "dimension checking" begin
        @test_throws DimensionMismatch polyadd([1, 2], [1, 2, 3])
        @test_throws DimensionMismatch polysub([1, 2], [1, 2, 3])
        @test_throws DimensionMismatch polymul([1, 2], [1, 2, 3])
        @test_throws DimensionMismatch polymulq([1, 2], [1, 2, 3])
    end
end
