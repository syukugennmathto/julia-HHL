# ntrugen.jl -- solving the NTRU equation  f*G - g*F = q  in Z[x]/(x^n + 1).
#
# ---------------------------------------------------------------------------
# What is being computed, and why it is the hard part
# ---------------------------------------------------------------------------
#
# A FALCON private key is a short basis of the NTRU lattice:
#
#     B = [[g, -f], [G, -F]]
#
# and for B to be a basis of the *right* lattice -- the one whose public
# description is h = g/f mod q -- the four polynomials must satisfy
#
#     f*G - g*F = q        (in Z[x]/(x^n+1), exactly, over the integers)
#
# Key generation samples the short pair (f, g) from a Gaussian and then has to
# *complete* it to a basis by finding (F, G).  That completion is this module.
#
# Two things make it hard:
#
#   1. It is an equation over Z[x]/(x^n+1), not over a field.  There is no
#      "just invert the matrix".
#   2. The naive solutions are enormous.  (F, G) start out with coefficients
#      thousands of bits long, and must be reduced back down to the size of
#      (f, g) or the key is useless -- a long basis samples wide signatures.
#
# Compare with Dilithium, where key generation is: sample s1, s2, compute
# t = A*s1 + s2, done.  There is no equation to solve, because Module-LWE
# does not ask for a *basis*, only for a hard instance.  FALCON's hash-and-sign
# structure needs the trapdoor to be an actual short basis, and paying for that
# is what this file is.
#
# ---------------------------------------------------------------------------
# The descent: the tower again
# ---------------------------------------------------------------------------
#
# The idea (Pornin-Prest) is to use the tower of subfields one last time.  The
# *field norm*
#
#     N : Q(zeta_2n) -> Q(zeta_n),    N(a) = a(x) * a(-x)
#
# is multiplicative, and a(x)*a(-x) is even, so it can be read as an element of
# the half-degree ring.  Applying N to the equation turns a problem in
# dimension n into one in dimension n/2:
#
#     solve  N(f)*G' - N(g)*F' = q   in dimension n/2,
#     then lift:  F = F'(x^2) * g(-x),   G = G'(x^2) * f(-x)
#
# The lift works because N(f) = f(x) f(-x), so
#
#     f*G - g*F = f * G'(x^2) f(-x) - g * F'(x^2) g(-x)
#               = N(f) G'(x^2) - N(g) F'(x^2)  ... as elements of the big ring
#
# and the bracket is exactly the half-degree equation, lifted.
#
# The recursion bottoms out at n = 1, where Z[x]/(x+1) = Z and the equation is
# just  f0*G0 - g0*F0 = q  -- extended Euclid.  This is the only place the
# solution is *constructed*; everything above is lifting and reducing.
#
# ---------------------------------------------------------------------------
# Where the multiprecision is, and why
# ---------------------------------------------------------------------------
#
# This is the question module 5 was set up to answer.  Sizes, at n = 512:
#
#   * f, g have coefficients of a few bits (sigma_fg ~ 4, so |c| < 30).
#   * Each descent step squares: N(f) has coefficients about twice the bit
#     length of f's.  After 9 levels the bottom values are ~2^9 times longer,
#     i.e. thousands of bits.
#   * At the bottom, extended Euclid multiplies by q and produces a solution
#     with no smallness guarantee at all.
#   * Lifting multiplies those by g(-x) / f(-x), growing them again.
#
# So the intermediate (F, G) genuinely need `BigInt`; there is no clever
# reordering that keeps them in 64 bits.  This is not "the numbers look big",
# it is a doubling per level over nine levels.
#
# The *reduction* is where floating point comes back, and the arrangement is
# worth studying, because it is not the obvious one.

# ---------------------------------------------------------------------------
# Karatsuba multiplication over Z
# ---------------------------------------------------------------------------
#
# The reference uses Karatsuba here rather than the schoolbook product of
# poly.jl, and says so explicitly: these multiplications are "more than 75% of
# the total cost in dimension n = 1024".  We follow it, and test that it agrees
# with `polymul` exactly -- which is the only thing that matters, since the two
# compute the same ring element by different routes.
#
# Note that the FFT is *not* usable here: the coefficients are thousands of
# bits long, far past the point where a Float64 product could be rounded back
# to the right integer (see the accuracy table in fft.jl).

"""
    karatsuba(a, b, n) -> Vector{BigInt}

Karatsuba product of two length-`n` polynomials, returning the full `2n`
coefficients *without* reduction modulo `x^n + 1`.

[Py-ref] scripts/pyref/ntrugen.py:14-38
"""
function karatsuba(a::AbstractVector{BigInt}, b::AbstractVector{BigInt}, n::Int)
    if n == 1
        return BigInt[a[1] * b[1], BigInt(0)]
    end
    m = n ÷ 2
    a0 = @view a[1:m];      a1 = @view a[(m + 1):n]
    b0 = @view b[1:m];      b1 = @view b[(m + 1):n]
    ax = BigInt[a0[i] + a1[i] for i in 1:m]
    bx = BigInt[b0[i] + b1[i] for i in 1:m]
    a0b0 = karatsuba(a0, b0, m)
    a1b1 = karatsuba(a1, b1, m)
    axbx = karatsuba(ax, bx, m)
    @inbounds for i in 1:n
        axbx[i] -= (a0b0[i] + a1b1[i])
    end
    ab = zeros(BigInt, 2n)
    @inbounds for i in 1:n
        ab[i] += a0b0[i]
        ab[i + n] += a1b1[i]
        ab[i + m] += axbx[i]
    end
    return ab
end

"""
    karamul(a, b) -> Vector{BigInt}

Product in `Z[x]/(x^n+1)`: the same ring element `polymul` computes.

[Py-ref] scripts/pyref/ntrugen.py:41-48

## Two routes, chosen by how big the coefficients actually are

The descent's coefficients span an enormous range -- 8 bits at the top,
3152 bits at the bottom (docs/debug_log.md #033) -- so a single representation
is wrong at one end or the other.  `BigInt` everywhere is correct and, at the
top, absurdly expensive: on values that fit an `Int64`, `BigInt` arithmetic
measures ~315x slower and allocates where machine integers allocate nothing.

So the width is decided per call, from the operands:

  * if `bits(a) + bits(b) + log2(n)` fits an `Int64`, multiply in `Int64`;
  * else if it fits an `Int128`, multiply in `Int128`;
  * else fall back to Karatsuba over `BigInt`, unchanged.

The bound is the exact one for a negacyclic convolution: each product needs
`bits(a) + bits(b)`, and at most `n` of them are summed into one coefficient.
`bitsize` rounds up to a whole byte, so the test is conservative in the safe
direction.  Nothing is assumed about the caller.

## Why this is worth a branch

`babai_reduce` is 96% of key generation, essentially all of it in one call at
n = 128, and that call runs 261 iterations of exactly two `karamul`s.  Measured
there, the operands are **24 bits** (`f`, `g`) and **24 bits** (`ki`, which is
`round`ed from a `Float64` and so cannot be wide), for a 56-bit product.  The
loop that dominates key generation was multiplying machine-word values through
GMP.

The fallback is not decoration.  At n = 16, 8 and 4 the coefficients measure
208, 408 and 808 bits and genuinely need `BigInt` -- but those levels run one
or two iterations, not 261.  "Machine words where they fit, multiprecision
where they do not" is not a compromise here; it matches the actual shape of the
problem.

## Why schoolbook and not Karatsuba on the fast path

Karatsuba wins on multiplication count and loses on allocation, and in this
implementation allocation is what costs (`polymulq_ntt` is 7.4x *slower* than
schoolbook at n = 512 for the same reason -- docs/debug_log.md #032).  The
machine-word path accumulates into one preallocated array and allocates nothing
else; Karatsuba would allocate O(n^1.58) temporaries to save multiplications
that are single instructions.
"""
function karamul(a::AbstractVector{BigInt}, b::AbstractVector{BigInt})
    _checklen(a, b)
    n = length(a)
    n == 0 && return BigInt[]

    # bitsize() rounds up to a byte, and log2(n) is rounded up too, so `need`
    # over-estimates: a call that passes the test is safe, and one that fails
    # it merely takes the slow route.
    need = maximum(bitsize, a) + maximum(bitsize, b) + (8 * sizeof(n) - leading_zeros(n))
    if need <= 62
        return BigInt.(_negacyclic_machine(Int64, a, b, n))
    elseif need <= 126
        return BigInt.(_negacyclic_machine(Int128, a, b, n))
    end

    ab = karatsuba(a, b, n)
    return BigInt[ab[i] - ab[i + n] for i in 1:n]     # x^n = -1
end

"""
Negacyclic convolution in a fixed-width integer type, with no allocation beyond
the two converted inputs and the result.  The caller has already established
that no coefficient can overflow `T`; there is deliberately no check here,
because a check per operation is what made the first attempt at this kind of
fast path worthless (docs/debug_log.md #034).
"""
function _negacyclic_machine(::Type{T}, a::AbstractVector{BigInt},
                             b::AbstractVector{BigInt}, n::Int) where {T<:Signed}
    av = T[T(c) for c in a]
    bv = T[T(c) for c in b]
    acc = zeros(T, n)
    @inbounds for i in 1:n
        ai = av[i]
        iszero(ai) && continue          # `ki` in babai_reduce is often sparse
        for j in 1:n
            k = i + j - 2               # degree of the term, zero-indexed
            p = ai * bv[j]
            if k < n
                acc[k + 1] += p
            else
                acc[k - n + 1] -= p     # x^n = -1
            end
        end
    end
    return acc
end

karamul(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer}) =
    karamul(BigInt.(a), BigInt.(b))

# ---------------------------------------------------------------------------
# The tower operations
# ---------------------------------------------------------------------------

"""
    galois_conjugate(a) -> Vector

The Galois conjugate `a(-x)`, i.e. flip the sign of every odd coefficient.

This is the non-trivial automorphism of `Q(zeta_2n)` over `Q(zeta_n)`: it fixes
the even part and negates the odd part, so its fixed field is exactly the
half-degree subfield.  Careful: this is **not** `polyadj` (which is
`a(1/x)`, the complex-conjugation adjoint).  Both are called "conjugate" in
the literature and they are different maps.

[Py-ref] scripts/pyref/ntrugen.py:51-57
"""
galois_conjugate(a::AbstractVector{T}) where {T} =
    T[isodd(i) ? a[i] : -a[i] for i in eachindex(a)]   # i is 1-based: odd i = even degree

"""
    field_norm(a) -> Vector{BigInt}

The relative field norm `N(a) = a(x) * a(-x)`, read as an element of the
half-degree ring `Z[x]/(x^{n/2} + 1)`.

Writing `a = a_e(x^2) + x a_o(x^2)`,

    a(x) a(-x) = a_e(x^2)^2 - x^2 a_o(x^2)^2

which involves only even powers, so substituting `y = x^2` gives

    N(a)(y) = a_e(y)^2 - y * a_o(y)^2      in Z[y]/(y^{n/2} + 1)

The `- y * a_o^2` is a shift by one degree, and the wrap-around picks up the
negacyclic sign -- which is why the last coefficient comes back **added** to
the constant term rather than subtracted.  That single `+` amid a run of `-`
is easy to mistype and produces a norm that is wrong only in one coefficient.

[Py-ref] scripts/pyref/ntrugen.py:60-74
"""
function field_norm(a::AbstractVector{BigInt})
    n2 = length(a) ÷ 2
    ae = BigInt[a[2i - 1] for i in 1:n2]
    ao = BigInt[a[2i] for i in 1:n2]
    res = karamul(ae, ae)
    aos = karamul(ao, ao)
    @inbounds for i in 1:(n2 - 1)
        res[i + 1] -= aos[i]
    end
    res[1] += aos[n2]                # x^{n/2} = -1, hence the sign flip
    return res
end

field_norm(a::AbstractVector{<:Integer}) = field_norm(BigInt.(a))

"""
    lift(a) -> Vector{BigInt}

Lift `a` from the half-degree ring to the full one: `a(x) |-> a(x^2)`, i.e.
spread the coefficients onto the even positions.

The right inverse of `polysplit`'s even part, and the operation that carries a
solution found one level down back up.

[Py-ref] scripts/pyref/ntrugen.py:77-86
"""
function lift(a::AbstractVector{T}) where {T}
    n = length(a)
    res = zeros(T, 2n)
    @inbounds for i in 1:n
        res[2i - 1] = a[i]
    end
    return res
end

# ---------------------------------------------------------------------------
# Extended Euclid at the bottom of the recursion
# ---------------------------------------------------------------------------

"""
    xgcd_floor(b, n) -> (d, u, v)

Extended GCD with `d = u*b + v*n`, using **floor** division and Python-style
remainder.

Julia's `gcdx` would do, except that we want the reference's *exact* Bezout
coefficients so that our `(F, G)` match its vectors byte for byte, and the
reference is written in Python: `b // n` floors and `b % n` takes the sign of
`n`.  Julia's `div`/`rem` truncate toward zero instead, which differs for
negative operands -- and `f0`, `g0` here are very much allowed to be negative.

So this uses `fld`/`mod`, which are Python's `//` and `%`.  Transliterating an
integer algorithm from Python without checking this is a reliable way to get a
correct-looking gcd with different cofactors.

[Py-ref] scripts/pyref/ntrugen.py:153-163
"""
function xgcd_floor(b::BigInt, n::BigInt)
    x0, x1, y0, y1 = BigInt(1), BigInt(0), BigInt(0), BigInt(1)
    while n != 0
        qq = fld(b, n)
        b, n = n, mod(b, n)
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    end
    return (b, x0, y0)
end

# ---------------------------------------------------------------------------
# Babai reduction
# ---------------------------------------------------------------------------

"""
    bitsize(a) -> Int

Bit length of `|a|`, **rounded up to a multiple of 8**; zero for `a = 0`.

The rounding is deliberate imprecision on the reference's part -- it makes the
function a byte count rather than a bit count, which is cheaper and is all the
scaling below needs.

[Py-ref] scripts/pyref/ntrugen.py:89-99

## This was the single most expensive line in key generation

The reference computes it by shifting a byte off at a time:

    val = abs(a);  res = 0
    while val != 0:  res += 8;  val >>= 8

which is fine in Python and a disaster here, for a reason that has nothing to
do with the loop being O(bits) -- it is that **every `>>=` on a `BigInt`
allocates a new one**.  `babai_reduce` calls `bitsize` on all 2n coefficients
of `(F, G)` on every iteration of its loop, and at n = 128 those coefficients
are 6232 bits, i.e. 779 shifts each.  That is ~200,000 `BigInt` allocations per
iteration, 261 iterations, and it accounted for essentially all of the 8 GB
that key generation was allocating (docs/debug_log.md #036).

GMP already knows the answer: `sizeinbase(x, 2)` reads it off the limb count in
constant time.  Same value, no allocation, no loop.

The `Base.GMP.MPZ` path is used only for `BigInt`; for machine integers
`ndigits` is already O(1) and needs no help.
"""
function bitsize(a::BigInt)
    iszero(a) && return 0
    # sizeinbase is exact for a power of two's *upper* bound and can read one
    # too high for other values, so it is not usable directly -- but the answer
    # is rounded up to a byte anyway, and `cld` absorbs the discrepancy except
    # exactly at a byte boundary.  Rather than reason about that, take the
    # exact bit length from Julia's own `ndigits`, which is also O(1) on BigInt.
    return 8 * cld(ndigits(a, base = 2), 8)
end

function bitsize(a::Integer)
    iszero(a) && return 0
    return 8 * cld(ndigits(abs(BigInt(a)), base = 2), 8)
end

"""
    babai_reduce(f, g, F, G) -> (F, G)

Reduce `(F, G)` modulo `(f, g)` by Babai's nearest-plane rounding:

    (F, G) <- (F, G) - k * (f, g),    k = round( (F f* + G g*) / (f f* + g g*) )

where `*` is the adjoint of poly.jl.  Geometrically, `k` is the coordinate of
`(F, G)` along `(f, g)` in the real embedding, so subtracting `k*(f, g)` makes
`(F, G)` as short as one step of nearest-plane can.  This does not change the
value of `f*G - g*F`, because the correction is a multiple of `(f, g)` itself
-- so the NTRU equation is preserved exactly while the coefficients shrink.

Corresponds to algorithm `Reduce` of the specification.
[Py-ref] scripts/pyref/ntrugen.py:102-149

## The floating-point arrangement, and why it is not fragile

`F` and `G` here have coefficients of thousands of bits, and `k` is computed in
Float64.  That looks alarming after fft.jl's accuracy table, which says the FFT
product goes wrong once coefficients pass ~1e7.  Three things make it fine, and
they are worth separating:

  1. **Only the ratio is computed.**  `k` is O(1)-sized even when `F` and `G`
     are astronomically large, because the numerator and denominator grow
     together.  What must be accurate is a quotient, not a product.

  2. **Both sides are pre-scaled to 53 bits.**  `f_adjust = f >> (size - 53)`
     keeps the top 53 bits of each coefficient, which is exactly the Float64
     mantissa, so the conversion to Float64 is *exact*.  The reference does the
     same to `F`, `G` with their own scale, and puts the difference of scales
     back with a shift when subtracting.  Nothing is ever rounded on the way in
     -- only the final `k` is rounded, to an integer, which is what we wanted
     anyway.

  3. **It iterates.**  The `while` loop recomputes `k` until it comes out zero
     or the size stops shrinking.  So an inaccurate `k` is not an error, it is
     a slower reduction: the next pass corrects it.  This is iterative
     refinement, and it is why FALCON can use double precision here without the
     spec having to pin down rounding behaviour at this particular step.

Contrast this with `ffSampling` (module 8), where the floating-point result is
*consumed directly* by a sampler and there is no correcting loop.  That is the
place where precision becomes a specification problem, not here.

CONSTANT TIME: irrelevant here -- key generation is not required to be constant
time in the same way signing is, since it runs once and its timing does not
correlate with a per-message secret.  The reference is nonetheless careful, and
the loop trip count does depend on the key.
"""
const BABAI_STEP = 25

# GMP's fused multiply-accumulate.  `Base.GMP.MPZ` wraps a lot of libgmp but
# not `addmul`/`submul`, and those two are exactly what this file wants: they
# do `r += a*b` and `r -= a*b` *into* `r`, with no temporary, where the Julia
# operators would allocate one GMP object for the product and another for the
# sum.  Measured, allocation of `BigInt`s was 55% of the whole descent
# (docs/debug_log.md #044), so removing the temporaries is the optimisation.
#
# `Ref{BigInt}` is the calling convention Base itself uses for `mpz_t`.
@inline _addmul_ui!(r::BigInt, a::BigInt, b::Culong) =
    ccall((:__gmpz_addmul_ui, :libgmp), Cvoid,
          (Ref{BigInt}, Ref{BigInt}, Culong), r, a, b)

@inline _submul_ui!(r::BigInt, a::BigInt, b::Culong) =
    ccall((:__gmpz_submul_ui, :libgmp), Cvoid,
          (Ref{BigInt}, Ref{BigInt}, Culong), r, a, b)

"""
    _negacyclic_addmul!(acc, f, ki) -> acc

`acc = f * ki` in `Z[x]/(x^n + 1)`, with `ki` a vector of machine integers,
computed in place with no allocation at all.

This is schoolbook, not Karatsuba, and that is deliberate.  `karamul` is
asymptotically better and is the right choice for two vectors of comparable,
large coefficients -- which is what the descent's `lift * conjugate` products
are.  Here one side is a 64-bit correction, so every partial product is a
`bignum * word`, which is exactly `mpz_addmul_ui`: one pass over the limbs,
accumulating in place.  Karatsuba's recursion would buy `n^1.58` instead of
`n^2` word-multiplications and pay for it with a tree of temporary `BigInt`s
and `SubArray`s, and at these sizes the allocation is the whole cost.

`ki` is usually sparse -- most corrections at a given scale are zero -- so the
outer loop skips zeros and the true cost is `nnz(ki) * n`.

## The machine-word tier is not optional

The first version of this was the GMP path alone, and it made the *deep* levels
(n <= 32, coefficients of hundreds of bits) ten times faster and the *top*
level (n = 512, coefficients of 8 bits) seven times slower -- 1.16 ms to 7.79.
The reason is that `karamul` has a machine-integer fast path, and at the top of
the descent it applies: 8-bit `f` times a 25-bit correction over 512 terms is
42 bits, so the whole convolution runs in `Int64` with no GMP at all.  Handing
that case to `mpz_addmul_ui` replaces an add-and-multiply with a library call.
So the same test `karamul` makes is made here, and the answer decides the tier.
Measured both ways at every level: docs/debug_log.md #044.

CONSTANT TIME: the zero skip is a data-dependent branch on a value derived from
the secret basis, and the tier choice is a branch on its magnitude.  Key
generation is not required to be constant time (see `babai_reduce`), and the
reference takes the same liberty, but it is worth naming: this loop's trip
count leaks the sparsity pattern of the corrections.
"""
function _negacyclic_addmul!(acc::Vector{BigInt}, f::Vector{BigInt}, ki::Vector{Int64},
                             fbits::Int = maximum(bitsize, f))
    n = length(f)
    length(acc) == n || throw(DimensionMismatch("acc and f must agree"))
    length(ki) == n || throw(DimensionMismatch("ki and f must agree"))

    kmax = zero(Int64)
    @inbounds for k in ki
        a = abs(k)
        kmax = ifelse(a > kmax, a, kmax)
    end
    # Same over-estimate as `karamul`: `bitsize` rounds up to a byte and the
    # length term is rounded up, so passing this test is safe.
    need = fbits + (64 - leading_zeros(kmax)) + (8 * sizeof(n) - leading_zeros(n))
    if need <= 62
        return _negacyclic_addmul_i64!(acc, f, ki, n)
    end

    @inbounds for i in 1:n
        Base.GMP.MPZ.set_si!(acc[i], 0)
    end
    @inbounds for l in 1:n
        k = ki[l]
        k == 0 && continue
        mag = Culong(abs(k))
        neg = k < 0
        # x^(l-1) * x^(j-1) = x^(l+j-2); the wrap at n costs a sign, since
        # x^n = -1 in this ring.
        base = l - 2
        for j in 1:n
            t = base + j
            if t < n
                if neg
                    _submul_ui!(acc[t + 1], f[j], mag)
                else
                    _addmul_ui!(acc[t + 1], f[j], mag)
                end
            else
                if neg
                    _addmul_ui!(acc[t - n + 1], f[j], mag)
                else
                    _submul_ui!(acc[t - n + 1], f[j], mag)
                end
            end
        end
    end
    return acc
end

"The `Int64` tier of [`_negacyclic_addmul!`](@ref).  No GMP, no allocation
beyond the two working arrays."
function _negacyclic_addmul_i64!(acc::Vector{BigInt}, f::Vector{BigInt},
                                 ki::Vector{Int64}, n::Int)
    fv = Vector{Int64}(undef, n)
    @inbounds for j in 1:n
        fv[j] = Int64(f[j])
    end
    av = zeros(Int64, n)
    @inbounds for l in 1:n
        k = ki[l]
        k == 0 && continue
        base = l - 2
        for j in 1:n
            t = base + j
            p = k * fv[j]
            if t < n
                av[t + 1] += p
            else
                av[t - n + 1] -= p
            end
        end
    end
    @inbounds for i in 1:n
        Base.GMP.MPZ.set_si!(acc[i], av[i])
    end
    return acc
end

function babai_reduce(f::Vector{BigInt}, g::Vector{BigInt},
                      F::Vector{BigInt}, G::Vector{BigInt};
                      step::Integer = BABAI_STEP)
    n = length(f)

    # A *deep* copy: the loop below mutates these in place, and `copy(F)`
    # duplicates the array but not the `BigInt`s in it, so a shallow copy would
    # corrupt the caller's (F, G) (docs/debug_log.md #040).
    F = BigInt[Base.GMP.MPZ.set!(BigInt(), c) for c in F]
    G = BigInt[Base.GMP.MPZ.set!(BigInt(), c) for c in G]

    # ---- the scale budget --------------------------------------------------
    #
    # This is the branch's whole point.  The specification's Reduce, which the
    # Python reference implements and `main` follows, derives the scaling of the
    # correction from the *current* sizes:
    #
    #     size  = max(53, bits(f), bits(g))
    #     shift = Size - size
    #     k     = round( (F f* + G g*) / (f f* + g g*) )   from 53-bit approximations
    #
    # and then stops when `k` rounds to zero.  Working through the magnitudes,
    # that `k` is about `2^(53 - bits(f))`.  Where f and g are *narrower* than
    # 53 bits the correction is large and the reduction grinds along; where they
    # are wider, `k` is about 2^0 and rounding it is a coin flip, so the
    # reduction gives up with (F, G) still enormous.  Measured, that happens at
    # every level below n = 128, and nine levels of reduction collapse into one
    # run on 6240-bit coefficients (docs/debug_log.md #041).
    #
    # The C reference instead carries an explicit *bit budget* and marches it
    # down: `k` is defined as the quotient divided by `2^scale_k`, with
    # `scale_k` starting at `bits(F,G) - bits(f,g)` and dropping by 25 each
    # pass.  Because the scaling comes from the schedule rather than from the
    # current sizes, `k` is a meaningful integer of about 25 bits every time and
    # never rounds away.  [C-ref] keygen.c:3165-3300 (`solve_NTRU_intermediate`)
    #
    # The subtraction is the same one either way -- `(F, G) -= k*(f, g)` scaled
    # -- so `f*G - g*F` is preserved exactly, and the result is still a valid
    # solution of the NTRU equation.  It is a *different* valid solution:
    # Babai reduction is not canonical, and this branch therefore does not
    # reproduce the Python reference's (F, G).  See README.md.
    bits_f = maximum(bitsize, f)
    bits_g = maximum(bitsize, g)
    size_fg = max(bits_f, bits_g)
    scale_fg = max(0, size_fg - 53)
    scale_k = max(0, max(maximum(bitsize, F), maximum(bitsize, G)) - size_fg)

    fa = Float64[Float64(c >> scale_fg) for c in f]
    ga = Float64[Float64(c >> scale_fg) for c in g]
    fa_fft = fft(fa)
    ga_fft = fft(ga)
    den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)),
                      mul_fft(ga_fft, adj_fft(ga_fft)))

    Fa = Vector{Float64}(undef, n)
    Ga = Vector{Float64}(undef, n)
    tmp = BigInt()
    scratch = BigInt()
    # `ki` in machine integers: the correction is bounded by the schedule, and
    # a `Vector{BigInt}` here would allocate n GMP objects per pass and then
    # hand them to `karamul`, which allocates a tree of temporaries of its own.
    ki = Vector{Int64}(undef, n)
    accF = BigInt[BigInt() for _ in 1:n]
    accG = BigInt[BigInt() for _ in 1:n]

    while true
        Size = max(maximum(bitsize, F), maximum(bitsize, G))
        scale_FG = max(0, Size - 53)
        @inbounds for i in 1:n
            Base.GMP.MPZ.fdiv_q_2exp!(scratch, F[i], scale_FG); Fa[i] = Float64(scratch)
            Base.GMP.MPZ.fdiv_q_2exp!(scratch, G[i], scale_FG); Ga[i] = Float64(scratch)
        end
        Fa_fft = fft(Fa)
        Ga_fft = fft(Ga)

        num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)),
                          mul_fft(Ga_fft, adj_fft(ga_fft)))
        ratio = ifft(div_fft(num_fft, den_fft))

        # `ratio` is the true quotient F/f scaled by 2^(scale_fg - scale_FG);
        # we want it scaled by 2^(-scale_k), so correct by 2^(-dc).
        dc = scale_k - scale_FG + scale_fg
        nonzero = false
        @inbounds for i in 1:n
            x = ldexp(real(ratio[i]), -dc)
            # A correction that does not fit an Int64 means the descent has
            # gone wrong for this (f, g); the reference bails out and lets key
            # generation resample, and so do we.
            isfinite(x) && abs(x) < 9.0e18 ||
                throw(NTRUSolveFailure("Babai correction out of range; resample f, g"))
            k = round(Int64, x)
            ki[i] = k
            nonzero |= (k != 0)
        end

        if nonzero
            _negacyclic_addmul!(accF, f, ki, bits_f)
            _negacyclic_addmul!(accG, g, ki, bits_g)
            @inbounds for i in 1:n
                if scale_k == 0
                    Base.GMP.MPZ.sub!(F[i], accF[i])
                    Base.GMP.MPZ.sub!(G[i], accG[i])
                else
                    Base.GMP.MPZ.mul_2exp!(tmp, accF[i], scale_k)
                    Base.GMP.MPZ.sub!(F[i], tmp)
                    Base.GMP.MPZ.mul_2exp!(tmp, accG[i], scale_k)
                    Base.GMP.MPZ.sub!(G[i], tmp)
                end
            end
        end

        # Unlike the specification's loop, an all-zero `k` is *not* a stopping
        # condition -- it just means this pass had nothing to remove at this
        # scale.  The loop is bounded by the schedule instead.
        scale_k <= 0 && break
        scale_k = max(0, scale_k - Int(step))
    end
    return (F, G)
end

# ---------------------------------------------------------------------------
# The recursion
# ---------------------------------------------------------------------------

"""
    NTRUSolveFailure

Raised when the NTRU equation has no solution for the given `(f, g)`.  Key
generation catches this and resamples.
"""
struct NTRUSolveFailure <: Exception
    msg::String
end
Base.showerror(io::IO, e::NTRUSolveFailure) = print(io, "NTRUSolveFailure: ", e.msg)

"""
    ntru_solve(f, g; q = Q) -> (F, G)

Solve `f*G - g*F = q` in `Z[x]/(x^n + 1)`.

Throws `NTRUSolveFailure` when the bottom-level `gcd(f0, g0)` is not 1, in
which case key generation resamples `(f, g)`.  That is a genuine possibility,
not a defensive check: `N^{log n}(f)` and `N^{log n}(g)` are two integers with
no reason to be coprime.

Corresponds to `NTRUSolve` of the specification.
[Py-ref] scripts/pyref/ntrugen.py:166-186
"""
function ntru_solve(f::Vector{BigInt}, g::Vector{BigInt}; q::Integer = Q,
                    step::Integer = BABAI_STEP)
    _checklen(f, g)
    n = length(f)
    if n == 1
        d, u, v = xgcd_floor(f[1], g[1])
        d == 1 || throw(NTRUSolveFailure(
            "gcd(N(f), N(g)) = $d != 1 at the bottom of the descent; resample f, g"))
        # f0*u + g0*v = 1, so f0*(q*u) - g0*(-q*v) = q.
        return (BigInt[-q * v], BigInt[q * u])
    end
    fp = field_norm(f)
    gp = field_norm(g)
    Fp, Gp = ntru_solve(fp, gp; q = q, step = step)
    F = karamul(lift(Fp), galois_conjugate(g))
    G = karamul(lift(Gp), galois_conjugate(f))
    return babai_reduce(f, g, F, G; step = step)
end

ntru_solve(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}; q::Integer = Q,
           step::Integer = BABAI_STEP) =
    ntru_solve(BigInt.(f), BigInt.(g); q = q, step = step)

"""
    ntru_equation_residual(f, g, F, G) -> Vector{BigInt}

`f*G - g*F` computed exactly in `BigInt`.  A correct solution makes this the
constant polynomial `q`.

This exists as a function rather than only as a test because it is the one
check that is worth running in anger: it is exact, it is cheap next to the
solve itself, and a key that fails it is silently broken in a way that only
shows up as "verification always fails" much later.
"""
ntru_equation_residual(f, g, F, G) =
    polysub(karamul(BigInt.(f), BigInt.(G)), karamul(BigInt.(g), BigInt.(F)))

"""
    ntru_equation_holds(f, g, F, G; q = Q) -> Bool

Whether `f*G - g*F = q` exactly.
"""
function ntru_equation_holds(f, g, F, G; q::Integer = Q)
    r = ntru_equation_residual(f, g, F, G)
    return r[1] == q && all(iszero, @view r[2:end])
end

# ---------------------------------------------------------------------------
# The key-generation quality test
# ---------------------------------------------------------------------------

"""
    gs_norm(f, g; q = Q) -> Float64

The squared Gram-Schmidt norm of the NTRU basis `[[g, -f], [G, -F]]`, computed
without needing `F` and `G`.

The Gram-Schmidt vectors of that basis are `(g, -f)` itself and, for the second
row, its component orthogonal to the first.  A short computation (this is the
content of the specification's line 9 of `NTRUGen`) gives the second norm as
`q^2 * || (adj(g), adj(f)) / (f adj(f) + g adj(g)) ||^2`, so the answer is

    max( ||(f, g)||^2,  q^2 * ||(Ft, Gt)||^2 )

with `Ft = adj(g)/(f adj f + g adj g)` and `Gt = adj(f)/(...)`.

Key generation rejects `(f, g)` unless this is at most `1.17^2 * q`
(`gram_schmidt_quality()`), which is the condition that makes the signature
width `sigma` achievable.

[Py-ref] scripts/pyref/ntrugen.py:189-201

Computed in the FFT domain, in Float64.  That is safe here for the same reason
as in `babai_reduce`: the quantity is a ratio of like-sized things, and the
comparison it feeds has a comfortable margin -- a borderline `(f, g)` that
falls on the wrong side of the threshold because of rounding is merely
resampled, not accepted wrongly.  A key that squeaks through is still checked
by the exact `ntru_equation_holds` and by the signature bound at verification.
"""
function gs_norm(f::AbstractVector{<:Real}, g::AbstractVector{<:Real}; q::Integer = Q)
    ff = Float64.(f)
    gg = Float64.(g)
    sqnorm_fg = sum(abs2, ff) + sum(abs2, gg)

    F_fft = fft(ff)
    G_fft = fft(gg)
    ffgg = add_fft(mul_fft(F_fft, adj_fft(F_fft)), mul_fft(G_fft, adj_fft(G_fft)))
    Ft = ifft(div_fft(adj_fft(G_fft), ffgg))
    Gt = ifft(div_fft(adj_fft(F_fft), ffgg))
    sqnorm_FG = float(q)^2 * (sum(abs2, Ft) + sum(abs2, Gt))

    return max(sqnorm_fg, sqnorm_FG)
end

"""
    gs_norm_ok(f, g; q = Q) -> Bool

The key-generation acceptance test: `||B~||^2 <= 1.17^2 * q`.

[Py-ref] scripts/pyref/ntrugen.py:232
"""
gs_norm_ok(f, g; q::Integer = Q) =
    gs_norm(f, g; q = q) <= gram_schmidt_quality()^2 * q

# ---------------------------------------------------------------------------
# Key generation proper
# ---------------------------------------------------------------------------
#
# This is the entry point that ties module 6 to module 7: sample (f, g) from a
# Gaussian, test them, and complete them to a basis.  It lives here rather than
# in samplerz.jl because the interesting part is the *rejection*, not the
# sampling.

"""
    SIGMA_FG_MIN

The `sigmin` argument key generation passes to `samplerz`: `SIGMA_FG_BASE -
0.001`.

[Py-ref] scripts/pyref/ntrugen.py:211 (`samplerz(0, sigma, sigma - 0.001)`)

Not a derived quantity -- it is simply "a hair below sigma", chosen so that the
`sigmin < sigma` precondition of `samplerz` holds with room to spare.  Written
out because a reimplementation that passes `sigma` itself, or `sigma_min` from
the parameter set, would still produce a plausible Gaussian and a different
byte stream.
"""
const SIGMA_FG_MIN = SIGMA_FG_BASE - 0.001

"""
    gen_poly(n, randombytes) -> Vector{BigInt}

Sample a degree-`n` polynomial whose coefficients follow `D_{Z, 0, sigma_fg}`
with `sigma_fg = 1.17 * sqrt(q / (2n))`.

[Py-ref] scripts/pyref/ntrugen.py:204-217

## Why 4096 samples regardless of n

The reference always draws **4096** samples at the fixed width
`SIGMA_FG_BASE = 1.43300980528773`, then folds them in consecutive blocks of
`k = 4096/n`.  Summing `k` independent Gaussians of width `s` gives width
`sqrt(k)*s`, which works out to exactly `sigma_fg`.

Drawing `n` coefficients directly at `sigma_fg` would give the *same
distribution* and a *different byte stream*, so it would pass every statistical
test and fail every KAT.  This is the same class of trap as the reversed KAT
chunks of module 7: a change that is invisible to the mathematics and fatal to
reproducibility.

It also means key generation consumes randomness for 8192 `samplerz` calls
before it has even looked at the result -- at n = 512 the sampling dominates
everything except `ntru_solve`.
"""
function gen_poly(n::Integer, randombytes)
    n < 4096 || throw(ArgumentError("gen_poly requires n < 4096, got $n"))
    4096 % n == 0 || throw(ArgumentError("gen_poly requires n | 4096, got $n"))
    f0 = Vector{Int}(undef, 4096)
    @inbounds for i in 1:4096
        f0[i] = samplerz(0.0, SIGMA_FG_BASE, SIGMA_FG_MIN, randombytes)
    end
    k = 4096 ÷ n
    f = Vector{BigInt}(undef, n)
    @inbounds for i in 1:n
        acc = 0
        for j in 1:k
            acc += f0[(i - 1) * k + j]
        end
        f[i] = acc
    end
    return f
end

# ---------------------------------------------------------------------------
# The C reference's Gaussian sampler for f and g  (THIS BRANCH ONLY)
# ---------------------------------------------------------------------------
#
# `gen_poly` above follows the specification and the Python reference: draw 4096
# values from `samplerz` at `SIGMA_FG_BASE` and fold them 4096/n at a time.  It
# is correct and it is expensive -- 4096 runs of a rejection sampler with an
# exponential in it, per polynomial, and roughly nine candidate (f, g) pairs are
# thrown away per accepted key.
#
# The C reference gets the *same distribution* a different way: a cumulative
# distribution table for the n = 1024 Gaussian, summed `2^(10-logn)` times.  At
# n = 512 that is 1024 table draws instead of 4096 rejection-sampler runs.
#
# The distributions really are the same, which is worth checking rather than
# assuming.  C's table is for `sigma = 1.17*sqrt(q/(2N))` with N = 1024, i.e.
# 2.866; summing `2^(10-logn)` of them gives `2.866 * sqrt(2^(10-logn))`.  Ours
# is `SIGMA_FG_BASE * sqrt(4096/n)` = `1.433 * sqrt(4096/n)`.  At n = 1024 both
# are 2.866; at n = 512 both are 4.0538.  The test suite asserts this.
#
# It also fixes something else.  C forces the sum of the coefficients to be
# **odd**, so that `Res(f, x^n+1)` is odd and the binary GCD at the bottom of
# the descent cannot fail on a common factor of 2.  Measured here, a third of
# all descents were being thrown away on `gcd != 1` without it.
#
# [C-ref] keygen.c:2258-2264 (the table's definition and its sigma)
# [C-ref] keygen.c:2266-2292 (`gauss_1024_12289`, transcribed below)
# [C-ref] keygen.c:4095-4131 (`mkgauss`, `poly_small_mkgauss`)

"""
    GAUSS_1024_12289

Cumulative distribution table for the discrete Gaussian with
`sigma = 1.17*sqrt(q/(2N))`, `q = 12289`, `N = 1024`, scaled by `2^63`.

Entry 0 is `P(x = 0)`.  For `k > 0`, entry `k` is `P(x >= k+1 | x > 0)`.

Transcribed from the C reference, not derived: these are 27 exact 63-bit
integers and re-deriving them would introduce rounding differences.
[C-ref] keygen.c:2266-2292
"""
const GAUSS_1024_12289 = UInt64[
    0x11d137d82df2ab58, 0x590c40f63ff5f974, 0x3898e41d85b975b7,
    0x20a964ef50858ff9, 0x1107d1ae973857eb, 0x07fe1ec29220ea37,
    0x035dafcacd37a439, 0x0144d98306216d42, 0x006d6beeeaf81655,
    0x0020e1a00d6fa84c, 0x0008cdddcd9dda9c, 0x0002192fc3dcdcb4,
    0x000071dfcd3c57e9, 0x00001574938d76eb, 0x000003974b0c33e5,
    0x000000889d3da6fe, 0x0000001204ddc6cb, 0x000000021bd3b27a,
    0x0000000038091f5e, 0x0000000005287db0, 0x00000000006bc528,
    0x000000000007cbfb, 0x0000000000007ffc, 0x0000000000000746,
    0x000000000000005e, 0x0000000000000004, 0x0000000000000000
]

"""
    mkgauss(randombytes, logn) -> Int

One coefficient of `f` or `g`, distributed as the sum of `2^(10-logn)` draws
from [`GAUSS_1024_12289`](@ref).

Each draw consumes two 64-bit words, little-endian, exactly as the C reference
does: the first supplies the sign and decides whether the value is zero, the
second indexes into the table.

CONSTANT TIME: the C reference scans the whole table every draw and combines
with masks precisely so the running time does not depend on the value sampled.
That is reproduced here -- the loop has no early exit -- but Julia gives no
guarantee that `ifelse` compiles to a branch-free select, so this is *shaped*
like constant-time code without being it.  The distinction matters: `f` and `g`
are the secret key.

[C-ref] keygen.c:4095-4131 (`mkgauss`)
"""
function mkgauss(randombytes, logn::Integer)
    return mkgauss_u64(_u64_reader(randombytes), logn)
end

"Adapt a `randombytes`-style callable into a zero-argument 64-bit reader."
_u64_reader(randombytes) = function ()
    b = randombytes(8)
    r = UInt64(0)
    @inbounds for i in 1:8
        r |= UInt64(b[i]) << (8 * (i - 1))
    end
    return r
end

"""
    mkgauss_u64(next64, logn) -> Int

[`mkgauss`](@ref) over a source that hands out 64-bit words directly.

This is the form the sampler actually wants.  Going through
`randombytes(8) -> Vector{UInt8}` allocates twice per draw, and at 1024 draws
per polynomial and roughly fourteen candidate polynomials per accepted key that
came to 58284 allocations per key (docs/debug_log.md #043).
"""
function mkgauss_u64(next64, logn::Integer)
    g = 1 << (10 - Int(logn))
    val = 0
    for _ in 1:g
        r = next64()
        neg = Int(r >> 63)
        r &= ~(UInt64(1) << 63)
        # `fl` becomes true when the value is zero, i.e. when r < table[0].
        fl = r < @inbounds(GAUSS_1024_12289[1])

        r = next64()
        r &= ~(UInt64(1) << 63)
        v = 0
        @inbounds for k in 2:length(GAUSS_1024_12289)
            t = r >= GAUSS_1024_12289[k]
            v = ifelse(t & !fl, k - 1, v)      # first k with r >= table[k]
            fl |= t
        end
        val += neg == 1 ? -v : v
    end
    return val
end

"""
Wrap a `randombytes` callable in a buffer and hand out 64-bit words.

The buffer is refilled `_U64_BLOCK` bytes at a time.  Note that this changes
*when* the underlying source is asked for bytes, though not the values it
produces in sequence -- a replay source with an exact budget will be asked for
more than it holds, which is why this is used only on the CDT path, where no
recorded byte count exists to match.

512 is not a tunable: it is the largest single request the reference ChaCha20
PRNG will serve, because that is its internal buffer (`randombytes!` in
shake.jl throws above it).  Asking for 4096 -- which is what the first version
of this did -- fails outright.  See docs/debug_log.md #043.
"""
const _U64_BLOCK = 512

"""
    BufferedU64(randombytes)

A callable that hands out 64-bit words, refilling from `randombytes` in
[`_U64_BLOCK`](@ref)-byte blocks.

This is a `mutable struct` with a type parameter on the source, and not the
closure it started life as, for a reason worth writing down.  The closure
version reassigned its captured buffer (`buf = randombytes(...)`) from inside
the returned function; Julia answers that by boxing the capture in a
`Core.Box`, whose element type is `Any`.  Every subsequent `buf[i]` is then a
dynamic dispatch returning a boxed value.  Measured, that cost a factor of
about seven on the whole sampler (docs/debug_log.md #043).

A callable struct captures the same state with the same lifetime and stays
concretely typed, so `next()` inlines to a bounds check and a load.
"""
mutable struct BufferedU64{F}
    src::F
    buf::Vector{UInt8}
    pos::Int      # bytes of `buf` already consumed
end

BufferedU64(src) = BufferedU64(src, UInt8[], 0)

"""
Load a little-endian `UInt64` from `buf[i:i+7]`.

On a little-endian machine that is exactly the memory image, so it is one
unaligned 8-byte load; elsewhere it has to be assembled a byte at a time.  The
byte order is the C reference's, not the machine's, so the branch is on
`ENDIAN_BOM` and not left to `reinterpret`.
"""
@inline function _load_u64_le(buf::Vector{UInt8}, i::Int)
    @boundscheck checkbounds(buf, i:(i + 7))
    if ENDIAN_BOM == 0x04030201            # little-endian host
        return GC.@preserve buf unsafe_load(Ptr{UInt64}(pointer(buf, i)))
    else
        r = UInt64(0)
        @inbounds for k in 0:7
            r |= UInt64(buf[i + k]) << (8 * k)
        end
        return r
    end
end

@inline function (b::BufferedU64)()
    p = b.pos
    if p + 8 > length(b.buf)
        b.buf = b.src(_U64_BLOCK)
        p = 0
    end
    b.pos = p + 8
    return @inbounds _load_u64_le(b.buf, p + 1)
end

"Deprecated spelling kept so existing call sites read the same."
_buffered_u64_reader(randombytes) = BufferedU64(randombytes)

"""
    gen_poly_cdt(n, randombytes) -> Vector{BigInt}

`f` or `g`, sampled with [`mkgauss`](@ref) rather than by folding `samplerz`.

Two constraints from the C reference are enforced by resampling the offending
coefficient:

  * `|c| <= 127`, so the coefficient fits the byte the key format stores it in;
  * the **sum of all coefficients is odd**, which makes `Res(f, x^n+1)` odd and
    stops the binary GCD at the bottom of the descent failing on a factor of 2.

The second is the interesting one: without it, a third of all descents here
were being discarded on `gcd != 1`.

[C-ref] keygen.c:4095-4131 (`poly_small_mkgauss`)
"""
function gen_poly_cdt(n::Integer, randombytes)
    ni = Int(n)
    logn = trailing_zeros(ni)
    (1 << logn) == ni || throw(ArgumentError("n must be a power of two, got $ni"))
    logn <= 10 || throw(ArgumentError("gen_poly_cdt is defined for n <= 1024"))
    # Sampled into `Int` and converted at the end.  Writing straight into a
    # `Vector{BigInt}` would allocate a GMP object per coefficient inside the
    # rejection loop, which is the hot loop.
    fi = Vector{Int}(undef, ni)
    mod2 = 0
    # One buffered reader for the whole polynomial: the source is asked for
    # half a kilobyte at a time instead of eight bytes two thousand times over.
    next64 = BufferedU64(randombytes)
    for u in 1:ni
        while true
            s = mkgauss_u64(next64, logn)
            (-127 <= s <= 127) || continue
            if u == ni
                # the last coefficient must make the total odd
                (mod2 ⊻ (s & 1)) == 0 && continue
            else
                mod2 ⊻= (s & 1)
            end
            fi[u] = s
            break
        end
    end
    return BigInt[BigInt(c) for c in fi]
end

"""
    ntru_gen(n, randombytes; q = Q, max_attempts = 1000) -> (f, g, F, G)

Generate a complete FALCON private basis: sample `(f, g)`, reject until they
are good, and solve for `(F, G)`.

Corresponds to `NTRUGen` of the specification.
[Py-ref] scripts/pyref/ntrugen.py:220-245

## The three rejection conditions, in the reference's order

1. **`gs_norm(f, g) > 1.17^2 * q`** -- the basis would be too long, so the
   signature width `sigma` would not be achievable.  This is the condition that
   makes `1.17` appear in three places at once (params.jl).
2. **`f` not invertible mod q** -- then `h = g/f` does not exist and there is
   no public key.  Cheap to test: no NTT coefficient may vanish (ntt.jl).
3. **`ntru_solve` fails** -- the bottom-level gcd is not 1, so the equation has
   no solution at all.

None of the three consumes randomness, so their order does not affect the byte
stream; it is kept as the reference has it for legibility.

`max_attempts` is a guard against an infinite loop on a broken sampler; the
reference has no such bound. Hitting it is a bug, not bad luck: measured
acceptance is roughly one draw in three.

CONSTANT TIME: key generation is the one part of FALCON where variable time is
broadly accepted -- it runs once, and its timing does not correlate with any
per-message secret. The rejection loop above is proudly data-dependent.
"""
function ntru_gen(n::Integer, randombytes; q::Integer = Q, max_attempts::Integer = 1000,
                 sampler::Symbol = :cdt)
    # THIS BRANCH defaults to the C reference's CDT sampler rather than the
    # specification's fold of `samplerz`.  Same distribution, different
    # realisation; see `gen_poly_cdt` and README.md.
    #
    # `sampler = :spec` selects the specification's `gen_poly`, and with it this
    # function still reproduces the Python reference's (f, g, F, G) byte for
    # byte -- which is what the reference-comparison tests use.  Keeping both
    # reachable is the point: the fast path is checked by properties, and the
    # slow path keeps the recorded vectors meaningful.
    gen = if sampler === :cdt
        gen_poly_cdt
    elseif sampler === :spec
        gen_poly
    else
        throw(ArgumentError("sampler must be :cdt or :spec, got :$sampler"))
    end
    for _ in 1:max_attempts
        f = gen(n, randombytes)
        g = gen(n, randombytes)

        # The cheap test first, as the C reference does: if the plain squared
        # norm of (f, g) already exceeds the bound then the Gram-Schmidt norm
        # certainly does, and the FFT that `gs_norm` would run is wasted.  Same
        # predicate, in two stages.  [C-ref] keygen.c:4238-4243
        fi = Int.(f); gi = Int.(g)
        plain = sum(abs2, fi) + sum(abs2, gi)
        plain >= 16823 && continue
        gs_norm(Float64.(f), Float64.(g); q = q) > gram_schmidt_quality()^2 * q && continue
        is_invertible_zq(fi) || continue

        try
            F, G = ntru_solve(f, g; q = q)
            return (f, g, F, G)
        catch e
            e isa NTRUSolveFailure || rethrow()
            continue
        end
    end
    throw(ErrorException(
        "ntru_gen: no acceptable (f, g) in $max_attempts attempts -- " *
        "this indicates a broken sampler, not bad luck"))
end
