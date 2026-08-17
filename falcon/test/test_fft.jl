# test_fft.jl -- module 5.
#
# Everything here is approximate, so every comparison needs a tolerance, and
# every tolerance needs a justification.  The ones used below come from
# measurements, not from taste:
#
#   * against the *Python* reference: loose (1e-9 relative).  Not because our
#     FFT is inaccurate but because theirs is -- their root table is written
#     out to ~15 decimal digits, which is up to ~370 ulp of error at n = 512.
#     Measured disagreement at n <= 1024, |c| <= 100: 3.9e-12.
#   * against the *C* reference: tight (1e-10 relative).  Its table is given to
#     27 digits and is correctly rounded, so it agrees with us much better than
#     Python does.  Measured disagreement: 2.3e-13.
#   * round trip: 1e-9 relative.  Measured: 6.4e-14 at |c| <= 100.
#
# If a tolerance here ever needs loosening, that is a finding, not a chore:
# write it up in docs/debug_log.md before touching the number.

# Golden vectors are included by runtests.jl.

@testset "fft" begin

    @testset "root exponents and roots" begin
        for n in sort(collect(keys(ROOTS_C)))
            e = Falcon.fft_root_exponents(n)
            w = Falcon.fft_roots(n)
            @test length(e) == n
            @test all(isodd, e)                       # roots of x^n+1 are odd powers
            @test all(k -> abs(k) < n, e)             # the invariant the recursion keeps
            @test length(unique(e)) == n
            # w^n = -1, to within rounding
            @test all(z -> abs(z^n + 1) < 1e-10, w)
            # pairs are (w, -w), first child is the principal square root
            for i in 1:2:n
                @test abs(w[i] + w[i + 1]) < 1e-14
            end
            if n > 2
                parent = Falcon.fft_roots(n ÷ 2)
                for i in 1:(n ÷ 2)
                    @test abs(w[2i - 1]^2 - parent[i]) < 1e-13
                    @test abs(w[2i - 1] - sqrt(parent[i])) < 1e-13   # principal
                end
            end
            # conjugate symmetry: the second half conjugates the first
            m = n ÷ 2
            for j in 1:m
                @test abs(w[j + m] - conj(w[j])) < 1e-14
            end
        end
    end

    @testset "roots match the reference table (within its own precision)" begin
        for (n, want) in ROOTS_C
            got = Falcon.fft_roots(n)
            @test _maxabsdiff(got, want) < 1e-13
        end
        # ... and ours is the *more* accurate of the two: every root must sit
        # on the unit circle to a couple of ulp.
        for n in (8, 512, 1024)
            @test all(z -> abs(abs(z) - 1) < 4e-16, Falcon.fft_roots(n))
        end
    end

    @testset "the FFT/NTT orderings are NOT the same" begin
        # Recorded as a test because it is the natural wrong guess: the two
        # tables have the same tree *shape* but order each +- pair by different
        # criteria (smaller residue vs principal square root), so at some nodes
        # index j names opposite roots.  If a future change made them agree,
        # this test failing is the signal to re-read docs/debug_log.md #009
        # rather than to delete the test.
        # Put both exponent lists in the same units.  The NTT's zeta has order
        # 2048, i.e. angle pi/1024, so an NTT exponent e means angle e*pi/1024.
        # An FFT exponent k at degree n means angle k*pi/n, i.e. e = k*1024/n.
        for n in (4, 8, 16, 512)
            fft_units = [mod(k * (1024 ÷ n), 2048) for k in Falcon.fft_root_exponents(n)]
            ntt_units = Falcon.ntt_root_exponents(n)
            # Same set of roots ...
            @test sort(fft_units) == sort(ntt_units)
            # ... but not in the same order.
            @test fft_units != ntt_units
        end
    end

    @testset "forward FFT against the Python reference" begin
        for (f, want) in FFT_KAT
            got = fft(f)
            @test length(got) == length(want)
            @test _maxabsdiff(got, want) < 1e-9 * max(1.0, maximum(abs, want))
        end
    end

    @testset "FFT is evaluation at the roots" begin
        # The definition, by Horner -- the same check as in test_ntt.jl, and
        # for the same reason: a transform evaluating at the wrong points
        # round-trips perfectly.
        for (f, _) in FFT_KAT
            n = length(f)
            w = Falcon.fft_roots(n)
            got = fft(f)
            want = [_evalpoly_c(f, w[j]) for j in 1:n]
            @test _maxabsdiff(got, want) < 1e-9 * max(1.0, maximum(abs, want))
        end
    end

    @testset "conjugate symmetry of a real polynomial's FFT" begin
        for (f, _) in FFT_KAT
            n = length(f)
            m = n ÷ 2
            F = fft(f)
            for j in 1:m
                @test abs(F[j + m] - conj(F[j])) < 1e-9 * max(1.0, abs(F[j]))
            end
        end
    end

    @testset "round trip" begin
        for (f, F) in FFT_KAT
            @test _maxabsdiff(ifft(F), f) < 1e-9 * max(1.0, maximum(abs, f))
            @test _maxabsdiff(ifft(fft(f)), f) < 1e-9 * max(1.0, maximum(abs, f))
            @test ifft(fft(f)) isa Vector{Float64}
        end
    end

    @testset "against the C reference implementation" begin
        # These vectors come from the C reference built with its own config.h,
        # i.e. with FALCON_FPEMU = 1 -- the emulated floating point that the
        # reference forces on for bit-reproducible signing.  The mapping
        # between its representation and ours (half storage, split real/imag,
        # Gray-code permutation) is the interesting part; see from_c_fft.
        for (logn, input, cout) in FFT_C_KAT
            n = 1 << logn
            @test length(input) == n
            @test length(cout) == n

            want = Falcon.from_c_fft(cout)
            got = fft(Float64.(input))
            scale = max(1.0, maximum(abs, want))
            @test _maxabsdiff(got, want) < 1e-10 * scale

            # the packing is invertible
            @test _maxabsdiff(Falcon.to_c_fft(want), cout) < 1e-12 * scale
            @test _maxabsdiff(Falcon.from_c_fft(Falcon.to_c_fft(got)), got) < 1e-12 * scale

            # and the C form really is conjugate-symmetric when unpacked
            m = n ÷ 2
            for j in 1:m
                @test want[j + m] == conj(want[j])
            end
        end
    end

    @testset "split_fft / merge_fft" begin
        for (f, F) in FFT_KAT
            length(F) < 4 && continue
            f0, f1 = Falcon.split_fft(F)
            @test _maxabsdiff(Falcon.merge_fft(f0, f1), F) < 1e-9 * max(1.0, maximum(abs, F))

            # they correspond to polysplit/polymerge in the coefficient domain:
            # this is the identity the whole ffLDL tree rests on.
            g0, g1 = polysplit(f)
            @test _maxabsdiff(f0, fft(g0)) < 1e-9 * max(1.0, maximum(abs, f0))
            @test _maxabsdiff(f1, fft(g1)) < 1e-9 * max(1.0, maximum(abs, f1))
        end
    end

    @testset "the conjugate in split_fft is not a typo" begin
        # split_fft divides by w by multiplying by conj(w). Guard the two
        # plausible corruptions: using w itself, or dropping the factor.
        # Both must break the round trip, or the test above proves nothing.
        for (_, F) in FFT_KAT
            length(F) < 8 && continue
            n = length(F)
            w = Falcon.fft_roots(n)
            m = n ÷ 2
            f0 = ComplexF64[0.5 * (F[2i - 1] + F[2i]) for i in 1:m]
            bad = ComplexF64[0.5 * (F[2i - 1] - F[2i]) * w[2i - 1] for i in 1:m]  # w, not conj(w)
            good = ComplexF64[0.5 * (F[2i - 1] - F[2i]) * conj(w[2i - 1]) for i in 1:m]
            @test _maxabsdiff(Falcon.merge_fft(f0, good), F) < 1e-9 * max(1.0, maximum(abs, F))
            @test _maxabsdiff(Falcon.merge_fft(f0, bad), F) > 1e-6 * max(1.0, maximum(abs, F))
        end
    end

    @testset "FFT-domain arithmetic" begin
        for (f, g, wmul, wdiv, wadj) in FFT_OPS
            @test _maxabsdiff(polymul_fft(f, g), wmul) < 1e-8 * max(1.0, maximum(abs, wmul))
            @test _maxabsdiff(polydiv_fft(f, g), wdiv) < 1e-8 * max(1.0, maximum(abs, wdiv))
            @test _maxabsdiff(polyadj_fft(f), wadj) < 1e-8 * max(1.0, maximum(abs, wadj))

            F, G = fft(f), fft(g)
            @test _maxabsdiff(add_fft(F, G), fft(polyadd(f, g))) < 1e-8 * max(1.0, maximum(abs, F))
            @test _maxabsdiff(sub_fft(F, G), fft(polysub(f, g))) < 1e-8 * max(1.0, maximum(abs, F))
            @test _maxabsdiff(neg_fft(F), fft(polyneg(f))) < 1e-8 * max(1.0, maximum(abs, F))
            # (f/g)*g == f
            @test _maxabsdiff(mul_fft(div_fft(F, G), G), F) < 1e-8 * max(1.0, maximum(abs, F))
        end
    end

    @testset "the FFT product agrees with the exact one" begin
        # The cross-check across modules: poly.jl computes the negacyclic
        # product exactly, this module computes it in floating point.  At
        # these coefficient sizes the FFT must land within 1e-6 of the exact
        # integer -- vastly better than the 0.5 that would be needed to round
        # correctly.  The margin is the point: see the accuracy table in
        # src/fft.jl for where it runs out.
        for (f, g, _) in NEGACYCLIC_Z
            n = length(f)
            n < 4 && continue
            exact = polymul(f, g)
            maximum(abs, exact) > 1e12 && continue     # the big-coefficient case
            approx = polymul_fft(Float64.(f), Float64.(g))
            @test _maxabsdiff(approx, Float64.(exact)) <
                  1e-9 * max(1.0, maximum(abs, Float64.(exact)))
            # and it rounds back to the exact integers
            @test round.(Int, approx) == Int.(exact)
        end
    end

    @testset "adjoint agrees with poly.jl" begin
        for (f, want) in ADJ_Z
            got = polyadj_fft(Float64.(f))
            @test _maxabsdiff(got, Float64.(want)) < 1e-9 * max(1.0, maximum(abs, Float64.(want)))
            @test round.(Int, got) == Int.(want)
        end
        # adj_fft is conjugation, and f*adj(f) is real and non-negative
        for (f, _) in FFT_KAT
            F = fft(f)
            @test adj_fft(F) == conj.(F)
            ffa = mul_fft(F, adj_fft(F))
            @test all(z -> abs(imag(z)) < 1e-9 * max(1.0, abs(z)), ffa)
            @test all(z -> real(z) >= -1e-12, ffa)
        end
    end

    @testset "argument checking" begin
        @test_throws ArgumentError fft(Float64[1.0])
        @test_throws ArgumentError fft(Float64[1.0, 2.0, 3.0])
        @test_throws ArgumentError ifft(ComplexF64[1.0 + 0im])
        @test_throws ArgumentError fft(ComplexF64[1.0 + 0im, 2.0 + 0im])
        @test_throws DimensionMismatch add_fft(ComplexF64[1, 2], ComplexF64[1, 2, 3])
        @test_throws DimensionMismatch mul_fft(ComplexF64[1, 2], ComplexF64[1, 2, 3])
    end
end
