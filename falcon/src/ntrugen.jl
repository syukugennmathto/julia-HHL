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
function babai_reduce(f::Vector{BigInt}, g::Vector{BigInt},
                      F::Vector{BigInt}, G::Vector{BigInt})
    n = length(f)
    F = copy(F)
    G = copy(G)
    size = max(53, maximum(bitsize, f), maximum(bitsize, g))

    # Top 53 bits of each coefficient: exactly representable as Float64.
    # `>>` on a negative BigInt is an arithmetic shift (floor division by a
    # power of two), matching Python's `>>`.  Julia agrees with Python here;
    # it is `div` vs `fld` that disagree (see xgcd_floor).
    fa = Float64[Float64(c >> (size - 53)) for c in f]
    ga = Float64[Float64(c >> (size - 53)) for c in g]
    fa_fft = fft(fa)
    ga_fft = fft(ga)
    den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)),
                      mul_fft(ga_fft, adj_fft(ga_fft)))

    while true
        Size = max(53, maximum(bitsize, F), maximum(bitsize, G))
        Size < size && break

        Fa = Float64[Float64(c >> (Size - 53)) for c in F]
        Ga = Float64[Float64(c >> (Size - 53)) for c in G]
        Fa_fft = fft(Fa)
        Ga_fft = fft(Ga)

        num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)),
                          mul_fft(Ga_fft, adj_fft(ga_fft)))
        k = ifft(div_fft(num_fft, den_fft))
        ki = BigInt[BigInt(round(elt)) for elt in k]
        all(iszero, ki) && break

        fk = karamul(f, ki)
        gk = karamul(g, ki)
        shift = Size - size
        @inbounds for i in 1:n
            F[i] -= fk[i] << shift
            G[i] -= gk[i] << shift
        end
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
function ntru_solve(f::Vector{BigInt}, g::Vector{BigInt}; q::Integer = Q)
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
    Fp, Gp = ntru_solve(fp, gp; q = q)
    F = karamul(lift(Fp), galois_conjugate(g))
    G = karamul(lift(Gp), galois_conjugate(f))
    return babai_reduce(f, g, F, G)
end

ntru_solve(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}; q::Integer = Q) =
    ntru_solve(BigInt.(f), BigInt.(g); q = q)

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
function ntru_gen(n::Integer, randombytes; q::Integer = Q, max_attempts::Integer = 1000)
    for _ in 1:max_attempts
        f = gen_poly(n, randombytes)
        g = gen_poly(n, randombytes)

        gs_norm(Float64.(f), Float64.(g); q = q) > gram_schmidt_quality()^2 * q && continue
        is_invertible_zq(Int.(f)) || continue

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
