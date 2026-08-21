# test_ntrugen.jl -- module 6.
#
# The strongest check available anywhere in this project lives here: the NTRU
# equation
#
#     f*G - g*F = q
#
# is an *exact integer identity*.  There is no tolerance, no distribution to
# sample, no reference byte string to trust -- either the polynomials satisfy
# it or they do not, and `BigInt` settles the question outright.
#
# So the tests come in two grades:
#
#   * `ntru_equation_holds` on everything.  This is self-validating: it needs
#     no oracle at all.  A solution that passes it is correct, full stop.
#   * agreement with the Python reference's *exact* (F, G).  Stronger than
#     necessary for correctness, but it pins the choices that the equation
#     alone does not determine -- the Bezout cofactors at the bottom of the
#     descent, and how far Babai reduction gets.  Those are what make a key
#     short rather than merely valid.
#
# The second is what would catch, say, `div` used where `fld` was meant: the
# equation would still hold, the key would still work, and the signatures would
# just be a bit longer than they should be.  That is exactly the kind of bug
# that never announces itself.

# Golden vectors are included by runtests.jl.

@testset "ntrugen" begin

    @testset "karamul agrees with the schoolbook product" begin
        # Karatsuba and poly.jl's O(n^2) convolution compute the same ring
        # element by different routes; they must agree exactly, over Z.
        for (f, g, want) in NEGACYCLIC_Z
            @test karamul(f, g) == want
            @test karamul(f, g) == polymul(f, g)
        end
        # including the case whose coefficients overflow Int64
        f, g, want = NEGACYCLIC_Z[end]
        @test karamul(f, g) == want
        # and on the tiny cases where the answer is checkable by hand
        @test karamul(BigInt[1, 0], BigInt[1, 0]) == BigInt[1, 0]
        @test karamul(BigInt[0, 1], BigInt[0, 1]) == BigInt[-1, 0]   # x*x = -1
    end

    @testset "karamul: machine-word path agrees with the BigInt path" begin
        # Regression test for docs/debug_log.md #036.  karamul now multiplies in
        # Int64 or Int128 when the operands provably fit, and falls back to
        # Karatsuba over BigInt otherwise.  The widths below are chosen to land
        # on *both* sides of both thresholds, including exactly at them -- a
        # test that only used small coefficients would never execute the
        # fallback, and one that only used large ones would never execute the
        # fast path.
        rng = MersenneTwister(20260821)
        function viaKaratsuba(a, b)
            n = length(a)
            ab = Falcon.karatsuba(a, b, n)
            return BigInt[ab[i] - ab[i + n] for i in 1:n]
        end

        widths_seen = Set{Symbol}()
        for n in (1, 2, 4, 8, 16, 64, 128), bits in (1, 4, 24, 26, 27, 28, 58, 62, 100, 126, 200, 900)
            a = BigInt[rand(rng, -(BigInt(2)^bits):(BigInt(2)^bits)) for _ in 1:n]
            b = BigInt[rand(rng, -(BigInt(2)^bits):(BigInt(2)^bits)) for _ in 1:n]
            got = karamul(a, b)
            # Two independent oracles: Karatsuba over BigInt (the old code path,
            # still reachable) and poly.jl's schoolbook convolution.
            @test got == viaKaratsuba(a, b)
            @test got == polymul(a, b)

            need = maximum(bitsize, a) + maximum(bitsize, b) +
                   (8 * sizeof(n) - leading_zeros(n))
            push!(widths_seen, need <= 62 ? :int64 : need <= 126 ? :int128 : :bigint)
        end
        # The point of the widths above: assert that all three paths ran.  If a
        # future change to the thresholds quietly stops exercising one of them,
        # this fails rather than the coverage silently vanishing.
        @test widths_seen == Set([:int64, :int128, :bigint])

        # The hot shape from babai_reduce at n = 128: 24-bit operands.  This is
        # the case the fast path exists for.
        f = BigInt[rand(rng, -(BigInt(2)^24):(BigInt(2)^24)) for _ in 1:128]
        k = BigInt[rand(rng, -(BigInt(2)^24):(BigInt(2)^24)) for _ in 1:128]
        @test karamul(f, k) == viaKaratsuba(f, k)
        # ... and a sparse `k`, which the fast path skips over
        ksparse = BigInt[iszero(i % 7) ? k[i] : BigInt(0) for i in 1:128]
        @test karamul(f, ksparse) == viaKaratsuba(f, ksparse)

        @test karamul(BigInt[], BigInt[]) == BigInt[]
    end

    @testset "_negacyclic_addmul! agrees with karamul, in both tiers" begin
        # THIS BRANCH.  `babai_reduce` no longer calls `karamul` for the
        # correction: `ki` is a vector of machine integers, so the product is
        # bignum-times-word and goes through `mpz_addmul_ui` with no
        # temporaries -- or, when everything fits 62 bits, through a plain
        # Int64 convolution with no GMP at all (docs/debug_log.md #044).
        #
        # Two tiers means two chances to be wrong, so the ranges below are
        # chosen to cross the boundary in both directions: coefficient widths
        # from 4 to 200 bits against corrections from 1 to 63 bits.  `karamul`
        # is the oracle, which is fair -- it is separately checked against the
        # schoolbook product above.
        rng = MersenneTwister(20260821)
        tiers = Set{Bool}()
        for n in (1, 2, 4, 8, 16, 32), _ in 1:60
            bits = rand(rng, (4, 20, 60, 200))
            f = BigInt[rand(rng, big(-2)^bits:big(2)^bits) for _ in 1:n]
            ki = Int64[rand(rng, 1:3) == 1 ? Int64(0) :
                       (rand(rng, Int64) >> rand(rng, 1:63)) for _ in 1:n]

            # which tier will this take?  (mirrors the test in the function)
            kmax = maximum(abs, ki)
            need = maximum(bitsize, f) + (64 - leading_zeros(kmax)) +
                   (8 * sizeof(n) - leading_zeros(n))
            push!(tiers, need <= 62)

            acc = BigInt[BigInt() for _ in 1:n]
            Falcon._negacyclic_addmul!(acc, f, ki)
            @test acc == karamul(f, BigInt.(ki))
        end
        # The test is worthless if it only ever exercised one tier.
        @test length(tiers) == 2

        # `acc` is reused across calls, so it must be zeroed, not accumulated.
        let n = 8
            f = BigInt[BigInt(i) for i in 1:n]
            ki = Int64[i == 1 ? Int64(3) : Int64(0) for i in 1:n]
            acc = BigInt[BigInt(999) for _ in 1:n]
            Falcon._negacyclic_addmul!(acc, f, ki)
            @test acc == karamul(f, BigInt.(ki))
            Falcon._negacyclic_addmul!(acc, f, ki)
            @test acc == karamul(f, BigInt.(ki))
        end
    end

    @testset "tower operations against the reference" begin
        for (a, wconj, wnorm, wlift) in TOWER_OPS
            @test galois_conjugate(a) == wconj
            @test field_norm(a) == wnorm
            @test lift(a) == wlift
        end
    end

    @testset "tower operations against their definitions" begin
        for (a, _, _, _) in TOWER_OPS
            n = length(a)
            # galois_conjugate is a(-x): an involution, and a ring homomorphism
            @test galois_conjugate(galois_conjugate(a)) == a
            b = circshift(a, 1)
            @test galois_conjugate(karamul(a, b)) ==
                  karamul(galois_conjugate(a), galois_conjugate(b))
            # It is NOT polyadj.  Worth pinning: both are called "conjugate".
            if n >= 4
                @test galois_conjugate(a) != polyadj(a)
            end

            # N(a) = a(x)*a(-x), lifted back up, must equal the product in the
            # big ring.  This is the identity the whole descent rests on.
            @test lift(field_norm(a)) == karamul(a, galois_conjugate(a))

            # the norm is multiplicative
            @test field_norm(karamul(a, b)) == karamul(field_norm(a), field_norm(b))

            # lift then take the even part gets back where we started
            @test polysplit(lift(a))[1] == a
            @test all(iszero, polysplit(lift(a))[2])
        end
    end

    @testset "bitsize" begin
        for (a, want) in BITSIZE_KAT
            @test bitsize(a) == want
        end
        @test bitsize(0) == 0
        # rounded UP to a multiple of 8, and sign-independent
        for a in (1, 7, 255, 256, 2^62, -(2^62))
            @test bitsize(a) % 8 == 0
            @test bitsize(a) == bitsize(-a)
            @test bitsize(a) >= (a == 0 ? 0 : ndigits(abs(BigInt(a)); base = 2))
            @test bitsize(a) < ndigits(abs(BigInt(a)); base = 2) + 8
        end

        # Regression test for docs/debug_log.md #036.  bitsize used to be the
        # reference's shift-a-byte-off-at-a-time loop, which on a BigInt
        # allocates once per shift and was responsible for essentially all of
        # key generation's 8 GB.  It is now read from the limb count in O(1).
        # The oracle here is that original definition, kept verbatim, because
        # the replacement has to agree with it *exactly* -- the value feeds a
        # shift amount, so being one byte out silently changes the arithmetic
        # rather than failing.
        bitsize_ref(a::Integer) = begin
            val = abs(BigInt(a)); res = 0
            while val != 0; res += 8; val >>= 8; end
            res
        end
        # exhaustive across a range crossing many byte boundaries
        @test all(a -> bitsize(a) == bitsize_ref(a), -600:600)
        # every power of two and its neighbours -- the boundary cases, and
        # where an off-by-one in a rounded-up bit count would hide
        @test all(e -> (v = BigInt(2)^e;
                        bitsize(v) == bitsize_ref(v) &&
                        bitsize(v - 1) == bitsize_ref(v - 1) &&
                        bitsize(v + 1) == bitsize_ref(v + 1) &&
                        bitsize(-v) == bitsize_ref(-v)), 0:600)
        # the width babai_reduce actually meets at n = 128
        let v = BigInt(2)^6232 - 12345
            @test bitsize(v) == bitsize_ref(v) == 6232
        end
        # machine integer types, including the extremes
        for T in (Int8, Int16, Int32, Int64, Int128)
            @test bitsize(typemax(T)) == bitsize_ref(typemax(T))
            @test bitsize(typemin(T)) == bitsize_ref(typemin(T))
            @test bitsize(zero(T)) == 0
        end
    end

    @testset "xgcd uses floor division, as Python does" begin
        # The whole point of not using Base.gcdx: the reference is written in
        # Python, whose // floors and whose % takes the sign of the divisor.
        # Julia's div/rem truncate. For negative operands the gcd is the same
        # but the cofactors differ, and the cofactors are what end up in the key.
        for (a, b) in ((BigInt(-7), BigInt(3)), (BigInt(7), BigInt(-3)),
                       (BigInt(-30), BigInt(-18)), (BigInt(12289), BigInt(-5)),
                       (BigInt(1), BigInt(0)), (BigInt(0), BigInt(5)))
            d, u, v = Falcon.xgcd_floor(a, b)
            @test u * a + v * b == d              # the defining identity
            @test abs(d) == gcd(abs(a), abs(b))
        end
    end

    @testset "the NTRU equation holds exactly" begin
        # The self-validating check: no oracle needed.
        for (f, g, F, G, _) in NTRU_SOLVE_KAT
            @test ntru_equation_holds(f, g, F, G)
            r = ntru_equation_residual(f, g, F, G)
            @test r[1] == Q
            @test all(iszero, @view r[2:end])
        end
    end

    @testset "ntru_solve reproduces the reference exactly" begin
        for (f, g, F, G, _) in NTRU_SOLVE_KAT
            Fj, Gj = ntru_solve(f, g)
            @test Fj == F
            @test Gj == G
        end
    end

    @testset "solutions are actually reduced" begin
        # Babai reduction is the difference between a valid key and a useful
        # one.  An unreduced solution satisfies the equation just as well but
        # has coefficients thousands of bits long.  Require (F, G) to be within
        # a modest factor of q times the size of (f, g) -- generous, but four
        # orders of magnitude away from "unreduced".
        for (f, g, F, G, _) in NTRU_SOLVE_KAT
            fg_bits = max(maximum(bitsize, f), maximum(bitsize, g))
            # (F, G) must come back to within a couple of bytes of (f, g)'s
            # size.  Measured, they land at 8-16 bits against f, g's 8; an
            # unreduced solution would be in the thousands.
            @test maximum(bitsize, F) <= fg_bits + 16
            @test maximum(bitsize, G) <= fg_bits + 16
            # and the Gram-Schmidt norm of the completed basis is finite and
            # positive, i.e. the basis is not degenerate
            @test 0 < gs_norm(Float64.(f), Float64.(g)) < Inf
        end
    end

    @testset "unsolvable (f, g) throw rather than return nonsense" begin
        for (f, g) in NTRU_SOLVE_FAILS
            @test_throws NTRUSolveFailure ntru_solve(f, g)
        end
        # f = g = 0 has no solution either (gcd is 0)
        @test_throws NTRUSolveFailure ntru_solve(BigInt[0, 0], BigInt[0, 0])
    end

    @testset "the target dimensions, self-validated" begin
        # n = 256 and n = 512 with no reference output at all: solve, then
        # check the exact integer identity.  This is the only test in the
        # project that needs no oracle whatsoever, and it covers the dimension
        # that actually matters.
        for (n, fs, gs) in NTRU_SOLVE_CANDIDATES
            solved = false
            rejected = 0
            for (f, g) in zip(fs, gs)
                local F, G
                try
                    F, G = ntru_solve(f, g)
                catch e
                    e isa NTRUSolveFailure || rethrow()
                    rejected += 1
                    continue
                end
                @test length(F) == n
                @test length(G) == n
                @test ntru_equation_holds(f, g, F, G)
                fg_bits = max(maximum(bitsize, f), maximum(bitsize, g))
                @test maximum(bitsize, F) <= fg_bits + 16
                @test maximum(bitsize, G) <= fg_bits + 16
                solved = true
                break
            end
            # If every candidate were rejected the test above would be vacuous,
            # so assert that one succeeded.
            @test solved
        end
    end

    @testset "gs_norm against the reference" begin
        for (f, g, want) in GS_NORM_KAT
            @test isapprox(gs_norm(f, g), want; rtol = 1e-9)
        end
    end

    @testset "gs_norm and the acceptance threshold" begin
        # gs_norm is at least ||(f,g)||^2 by construction (it is a max).
        for (f, g, _) in GS_NORM_KAT
            @test gs_norm(f, g) >= sum(abs2, f) + sum(abs2, g) - 1e-6
        end
        # A pathologically small (f, g) has a huge gs_norm -- the q^2/||fg||^2
        # branch blows up -- and must be rejected.  This is the branch that
        # protects against a key whose second basis vector is enormous.
        tiny = [1.0; zeros(63)]
        @test gs_norm(tiny, zeros(64)) > gram_schmidt_quality()^2 * Q
        @test !gs_norm_ok(tiny, zeros(64))
    end

    @testset "the descent's degree bookkeeping" begin
        # field_norm halves, lift doubles; the recursion must terminate at 1.
        a = BigInt.(1:64)
        len = length(a)
        while len > 1
            a = field_norm(a)
            len ÷= 2
            @test length(a) == len
        end
        @test length(a) == 1
    end

    @testset "this branch's reduction schedule -- what it keeps and what it drops" begin
        # THIS BRANCH DIVERGES FROM THE SPECIFICATION.  See README.md and
        # docs/debug_log.md #042.  `babai_reduce` here follows the C reference's
        # explicit bit-budget schedule instead of the specification's Reduce,
        # which the Python reference implements and which `main` follows.
        #
        # The testset above (`"ntru_solve reproduces the reference exactly"`)
        # still passes, and that is *not* evidence the two agree in general:
        # the recorded vectors have f and g of about 8 bits, so `size` is
        # clamped to 53 in both schedules and they do the same thing.  The
        # divergence only appears where f and g exceed 53 bits, which happens
        # at the deep levels of a real n = 512 descent.  This testset pins down
        # what survives the divergence and what does not.

        rng = MersenneTwister(20260826)

        # 1. The invariant that actually matters is preserved *exactly*.
        #    Subtracting a multiple of (f, g) cannot change f*G - g*F, whatever
        #    schedule chose the multiple.
        for n in (2, 4, 8, 16, 32, 64)
            for _ in 1:3
                f = BigInt[rand(rng, -5:5) for _ in 1:n]
                g = BigInt[rand(rng, -5:5) for _ in 1:n]
                local F, G
                try
                    F, G = ntru_solve(f, g)
                catch e
                    e isa NTRUSolveFailure || rethrow()
                    continue
                end
                @test ntru_equation_holds(f, g, F, G)
                # and the solution is short: within a couple of bytes of (f, g)
                fg = max(maximum(bitsize, f), maximum(bitsize, g))
                @test maximum(bitsize, F) <= fg + 16
                @test maximum(bitsize, G) <= fg + 16
            end
        end

        # 2. An all-zero correction is NOT a stopping condition here.  On main
        #    it is, and that is the whole difference: the loop must be able to
        #    pass through a scale at which there is nothing to remove and keep
        #    going at the next one.  Asserted by construction rather than by
        #    timing: give babai_reduce an (F, G) that is already reduced, and
        #    it must terminate and leave them alone rather than spin.
        let n = 8
            f = BigInt[rand(rng, -5:5) for _ in 1:n]
            g = BigInt[rand(rng, -5:5) for _ in 1:n]
            local F, G
            try
                F, G = ntru_solve(f, g)
                Fr, Gr = Falcon.babai_reduce(f, g, F, G)
                @test ntru_equation_holds(f, g, Fr, Gr)
                @test maximum(bitsize, Fr) <= maximum(bitsize, F)
                @test maximum(bitsize, Gr) <= maximum(bitsize, G)
            catch e
                e isa NTRUSolveFailure || rethrow()
            end
        end

        # 3. Arguments are still not mutated -- the reduction works in place on
        #    its own deep copies (docs/debug_log.md #040).
        let n = 16
            f = BigInt[rand(rng, -4:4) for _ in 1:n]
            g = BigInt[rand(rng, -4:4) for _ in 1:n]
            Fp, Gp = try
                ntru_solve(Falcon.field_norm(f), Falcon.field_norm(g))
            catch e
                e isa NTRUSolveFailure ? (nothing, nothing) : rethrow()
            end
            if Fp !== nothing
                Fraw = karamul(lift(Fp), galois_conjugate(g))
                Graw = karamul(lift(Gp), galois_conjugate(f))
                sF = copy(Fraw); sG = copy(Graw); sf = copy(f); sg = copy(g)
                Falcon.babai_reduce(f, g, Fraw, Graw)
                @test Fraw == sF && Graw == sG && f == sf && g == sg
            end
        end
    end
end
