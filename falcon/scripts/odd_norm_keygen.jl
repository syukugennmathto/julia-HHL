#!/usr/bin/env julia
#
# odd_norm_keygen.jl -- is part 2 of ePrint 2024/1709's countermeasure actually
# undeployable, or is it "easily fixable" as that paper says?
#
#     julia --project=falcon falcon/scripts/odd_norm_keygen.jl [tries]
#
# Their section 7.1 needs ||(g,-f)||^2 ODD.  The reference key generator never
# produces one, and the reason is exact rather than accidental:
# `gen_poly_cdt` (src/ntrugen.jl:1179, following keygen.c's poly_small_mkgauss)
# forces the coefficient sum of BOTH f and g to be odd, so
#
#     ||(g,-f)||^2 = sum g_i^2 + sum f_i^2 = g(1) + f(1) = 1 + 1 = 0  (mod 2)
#
# -- always even, by construction, in one line.
#
# The stated reason for the constraint is that it makes Res(f, x^n+1) odd, so
# the binary GCD at the bottom of the Pornin-Prest field-norm descent does not
# fail on a factor of 2.  2024/1709 argues only ONE of the two needs to be odd
# ("there is no reason to require both"), which would give an odd norm.
#
# This script tests that claim directly in a working implementation: sample f
# with an odd coefficient sum and g with an EVEN one, run the real NTRU solver,
# and report how often it succeeds and what the resulting norm parity is.  The
# answer decides how strong the countermeasure-deployability claim may be.

using Falcon
using Printf

const F = Falcon

# A parity-controlled version of the reference's small-polynomial sampler.
Falcon.eval(quote
    const WANT_ODD = Ref(true)      # parity demanded of the NEXT polynomial
    const ALTERNATE = Ref(false)    # if set, flip WANT_ODD after each call
    function gen_poly_cdt(n::Integer, randombytes)
        ni = Int(n); logn = trailing_zeros(ni)
        fi = Vector{Int}(undef, ni); mod2 = 0
        want = WANT_ODD[] ? 1 : 0
        next64 = BufferedU64(randombytes)
        for u in 1:ni
            while true
                s = mkgauss_u64(next64, logn)
                (-127 <= s <= 127) || continue
                if u == ni
                    (mod2 ⊻ (s & 1)) != want && continue
                else
                    mod2 ⊻= (s & 1)
                end
                fi[u] = s; break
            end
        end
        ALTERNATE[] && (WANT_ODD[] = !WANT_ODD[])
        return BigInt[BigInt(c) for c in fi]
    end
end)

function main()
    tries = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 60
    n = 512
    println("# odd_norm_keygen.jl -- can the solver produce an odd ||(g,-f)||^2?")
    println()

    # ---- control: the reference behaviour, both sums odd ------------------
    F.ALTERNATE[] = false; F.WANT_ODD[] = true
    r = chacha20(shake256(codeunits("onk/control"), 56))
    okc = 0; oddc = 0
    for _ in 1:tries
        sk = try falcon_keygen(n, k -> randombytes!(r, k))[1] catch; nothing end
        sk === nothing && continue
        okc += 1
        t = sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f)
        isodd(t) && (oddc += 1)
        @assert isodd(sum(Int.(sk.f))) && isodd(sum(Int.(sk.g)))
    end
    @printf("control  (f odd, g odd -- the reference): %d/%d keys solved, %d with ODD norm\n",
            okc, tries, oddc)

    # ---- treatment: f odd, g even ----------------------------------------
    F.ALTERNATE[] = true; F.WANT_ODD[] = true     # f gets odd, g gets even, then repeats
    r2 = chacha20(shake256(codeunits("onk/treat"), 56))
    okt = 0; oddt = 0; fails = 0
    for _ in 1:tries
        F.WANT_ODD[] = true                        # each attempt starts with f
        sk = try falcon_keygen(n, k -> randombytes!(r2, k))[1]
             catch e; fails += 1; nothing end
        sk === nothing && continue
        okt += 1
        t = sum(x -> Int128(x)^2, sk.g) + sum(x -> Int128(x)^2, sk.f)
        isodd(t) && (oddt += 1)
    end
    @printf("treatment(f odd, g even)                : %d/%d keys solved, %d with ODD norm",
            okt, tries, oddt)
    println(fails > 0 ? "   ($fails attempts threw)" : "")
    println()
    if oddt > 0
        println("=> The solver DOES produce odd-norm keys when only one parity is forced.")
        println("   ePrint 2024/1709's \"easily fixable\" is correct: part 2 of their")
        println("   countermeasure is deployable with a one-line change to key generation.")
        println("   Our claim must therefore be that part 1 ALONE -- the change an")
        println("   implementer makes if they touch only the sampler -- gives no")
        println("   protection at the key-recovery positions, NOT that the countermeasure")
        println("   is undeployable.")
    else
        println("=> No odd-norm key was produced even with one parity relaxed.")
        println("   Report the failure mode before drawing any conclusion.")
    end
    F.ALTERNATE[] = false; F.WANT_ODD[] = true
end

main()
