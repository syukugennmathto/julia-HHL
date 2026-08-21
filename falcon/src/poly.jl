# poly.jl -- arithmetic in the ring R = Z[x]/(x^n + 1) and in R_q = R/qR.
#
# This module is deliberately the *slow, obviously correct* one.  Everything
# here is schoolbook O(n^2); the fast paths live in ntt.jl (exact, mod q) and
# fft.jl (approximate, over C).  Having a naive implementation to compare
# against is not a courtesy to the reader, it is the debugging strategy: when
# the NTT disagrees with the definition, the definition wins.
#
# ---------------------------------------------------------------------------
# The mathematics
# ---------------------------------------------------------------------------
#
# R = Z[x]/(x^n + 1) with n a power of two.  Because x^n + 1 is the 2n-th
# cyclotomic polynomial for such n, R is the ring of integers of the
# cyclotomic field Q(zeta_{2n}), and every structural fact FALCON uses -- the
# tower of subfields, the field norm, the FFT/NTT factorisations -- is a fact
# about that field.
#
# Concretely, reducing modulo x^n + 1 means x^n = -1, so a product wraps around
# *with a sign flip*:
#
#     (f * g)_k = sum_{i+j=k} f_i g_j  -  sum_{i+j=k+n} f_i g_j
#
# This is the "negacyclic" convolution, as opposed to the cyclic convolution of
# Z[x]/(x^n - 1).  That single minus sign is the source of a whole family of
# bugs, because it reappears independently in three places:
#
#     * here, in the schoolbook product;
#     * in ntt.jl, as the requirement for a primitive *2n*-th root of unity
#       (the roots of x^n + 1 are the odd powers of zeta_{2n}, not the n-th
#       roots of unity);
#     * in fft.jl, as the twiddle factors of `splitfft`/`mergefft`.
#
# Fix it in one place and you have fixed one third of it.  This is worth
# saying out loud now, because it is the kind of thing one rediscovers at 2am.
#
# ---------------------------------------------------------------------------
# Comparison with Dilithium / ML-DSA
# ---------------------------------------------------------------------------
#
# Dilithium uses the same *shape* of ring, Z_q[x]/(x^256 + 1), and the same
# negacyclic convolution.  Two differences matter for the implementation:
#
#   * Dilithium lives *entirely* modulo q.  Every polynomial it ever handles is
#     an element of R_q, and `poly.jl` there can be a thin wrapper over the NTT.
#     FALCON needs polynomials over Z with genuinely large coefficients (in the
#     NTRU solving step), over Q, and over C (in the FFT).  So this module is
#     generic in the coefficient type, and that genericity is load-bearing.
#
#   * Dilithium never needs the adjoint.  FALCON does, because its security
#     argument is about a *Gram matrix* and its sampler is a Cholesky/LDL
#     decomposition of that matrix.  `polyadj` below is the ring-level
#     ingredient of that.

# ---------------------------------------------------------------------------
# Basic operations
# ---------------------------------------------------------------------------

"""
    polyadd(f, g)

Coefficient-wise sum in `R`.  Requires `length(f) == length(g)`.
"""
function polyadd(f::AbstractVector{T}, g::AbstractVector{T}) where {T}
    _checklen(f, g)
    return T[f[i] + g[i] for i in eachindex(f)]
end

"""
    polysub(f, g)

Coefficient-wise difference in `R`.
"""
function polysub(f::AbstractVector{T}, g::AbstractVector{T}) where {T}
    _checklen(f, g)
    return T[f[i] - g[i] for i in eachindex(f)]
end

"""
    polyneg(f)

Additive inverse in `R`.
"""
polyneg(f::AbstractVector{T}) where {T} = T[-c for c in f]

function _checklen(f, g)
    length(f) == length(g) ||
        throw(DimensionMismatch("polynomials have different degrees: " *
                                "$(length(f)) and $(length(g))"))
    return nothing
end

"""
    polymul(f, g)

Negacyclic product in `R = Z[x]/(x^n + 1)`, computed by the definition in
O(n^2) operations.

## Overflow

The accumulator is widened to `widen(T)` (so `Int64` accumulates in `Int128`),
because this is precisely where a careless implementation silently wraps
around.  With `n = 512` and coefficients bounded by `B`, the accumulator can
reach `n * B^2`; for `B ~ 10^9` -- an entirely realistic size for the
intermediate `F`, `G` of the NTRU solving step -- that is about `5 * 10^20`,
comfortably past `typemax(Int64) ~ 9.2 * 10^18`.

The result is then narrowed back to `T`, and the narrowing is *checked*: if a
coefficient does not fit, an `OverflowError` is raised rather than a wrong
answer returned.  Callers that expect large results should pass `BigInt`
coefficients (for which `widen` is the identity and nothing can overflow).

This is the first appearance of the multiprecision question that dominates
ntrugen.jl; see docs/math/03_poly.md for the size analysis.
"""
function polymul(f::AbstractVector{T}, g::AbstractVector{T}) where {T<:Integer}
    _checklen(f, g)
    n = length(f)
    A = widen(T)
    acc = zeros(A, n)
    @inbounds for i in 1:n, j in 1:n
        k = i + j - 2                      # degree of the term, zero-indexed
        p = A(f[i]) * A(g[j])
        if k < n
            acc[k + 1] += p
        else
            acc[k - n + 1] -= p            # x^n = -1
        end
    end
    return T[_checked_narrow(T, c) for c in acc]
end

function _checked_narrow(::Type{T}, c) where {T<:Integer}
    (typemin(T) <= c <= typemax(T)) || throw(OverflowError(
        "coefficient $c does not fit in $T; use BigInt coefficients here"))
    return T(c)
end
_checked_narrow(::Type{BigInt}, c) = BigInt(c)

"""
    polyadj(f)

The adjoint (Galois conjugate) of `f` in `R`:

    adj(f)(x) = f(x^{-1}) = f_0 - f_{n-1} x - f_{n-2} x^2 - ... - f_1 x^{n-1}

i.e. `adj(f)_0 = f_0` and `adj(f)_i = -f_{n-i}` for `i >= 1`.

## Why the ring has an adjoint at all, and why FALCON needs it

Embed `R` into `C^n` by evaluating at the `n` primitive `2n`-th roots of unity
(that is exactly what `fft.jl` will do).  Under that embedding, multiplication
by `f` becomes multiplication by the vector of its evaluations, and `adj(f)`
becomes the *complex conjugate* of that vector -- because the roots come in
conjugate pairs and `conj(zeta) = zeta^{-1}`.

So `adj` is the Hermitian adjoint of "multiply by `f`":

    <f * a, b> = <a, adj(f) * b>

and consequently `f * adj(f)` is a self-adjoint, positive element -- the
polynomial analogue of `|f|^2`.  Every Gram matrix in FALCON is built out of
such products, which is why LDL* over `R` makes sense at all.  Getting the
sign of the `i >= 1` branch wrong gives a Gram matrix that is not Hermitian,
and then the LDL decomposition produces a "sigma" that is negative, and then
`sqrt` produces a `NaN`, and only then do you find out.
"""
function polyadj(f::AbstractVector{T}) where {T}
    n = length(f)
    out = Vector{T}(undef, n)
    out[1] = f[1]
    @inbounds for i in 2:n
        out[i] = -f[n - i + 2]
    end
    return out
end

"""
    sqnorm(vs...) -> BigInt

Squared Euclidean norm of the concatenation of the given coefficient vectors,
computed **exactly** and returned as a `BigInt` whatever the input type.

Exact arithmetic is not paranoia here: the key-generation rejection test
compares a Gram-Schmidt norm against `1.17^2 * q`, and the verification test
compares `||(s1, s2)||^2` against `beta^2`.  Both are *decisions*, so a
rounding error changes the answer rather than perturbing it.

[Py-ref] scripts/pyref/common.py:38-44 (`sqnorm`)

## Why there are two loops

The obvious implementation -- accumulate straight into a `BigInt` -- was
measured at **46% of the entire cost of `falcon_verify`** and 14% of signing
(docs/debug_log.md #033).  That is not because the numbers are big.  In
verification they are tiny: the coefficients are centred mod `q`, so `|c| <=
6144`, and the whole sum is under `2^38`.  It is because `BigInt` in Julia is a
*heap object* wrapping a GMP `mpz_t`.  Every `+` and every `*` allocates, so a
loop over 2n coefficients allocates 4n times and then pays the collector.
Measured on values that fit an `Int64` either way, `BigInt` arithmetic is ~315x
slower than `Int64` and ~40x slower than `Int128`.

So the fast path accumulates in `Int128` with *checked* arithmetic and the slow
path is kept, unchanged, for anything that does not fit.  The fallback is not
decoration: `ntrugen.jl` legitimately handles coefficients thousands of bits
long, and `Int128(c)` on one of those throws `InexactError` before any wrong
answer can be produced.

## Why the guard is a range test and not checked arithmetic

The first version of this used `Base.Checked.checked_mul`/`checked_add` on the
`Int128` accumulator, which is the obvious way to be safe.  It gave back only
3.4x of the 300x, and the reason is worth knowing: **`Int128` has no hardware
overflow flag**, so Julia's checked operations on it are emulated in software.
Measured, checked `Int128` arithmetic is ~83x slower than unchecked.  (A
`try`/`catch` around the loop, the other suspect, costs nothing at all here.)

So the overflow argument is made *once, about the inputs*, instead of on every
operation.  If every coefficient satisfies `|c| <= 2^40` then each square is at
most `2^80`, and the sum of fewer than `2^47` of them stays under `2^127`.  Two
comparisons per coefficient replace two emulated 128-bit checked operations,
and the bound is checked rather than assumed -- anything outside it takes the
`BigInt` path and is still exact.

`2^40` is not a tuned number; it is simply far above anything a signature can
contain (`|c| <= q/2 = 6144`, about `2^12.6`) and far below where the
accumulator could be troubled.

CONSTANT TIME: the fast path is taken or not depending on the magnitude of the
coefficients, so the running time depends on the data.  For verification that
is harmless -- everything here is public -- but the same function is called on
a *secret* basis during key generation, and there the branch is a leak.  A
production implementation would fix the width by parameter set and never
branch.
"""
function sqnorm(vs::AbstractVector...)
    acc = _sqnorm_i128(vs)
    acc === nothing || return BigInt(acc)
    # Fallback: exact, unconditional, and the only path for coefficients that
    # are genuinely large (ntrugen.jl's descent reaches thousands of bits).
    slow = BigInt(0)
    for v in vs, c in v
        slow += BigInt(c) * BigInt(c)
    end
    return slow
end

"""
    sqnorm_machine(v1, v2) -> Int128

`||v1||^2 + ||v2||^2` for two vectors of machine integers, with **no branch on
the data and no allocation**.

Used by signing's acceptance test.  `sqnorm` above dispatches between an
`Int128` accumulator and a `BigInt` one by *scanning the coefficients*, which
is correct in general -- the descent in ntrugen.jl really does reach thousands
of bits -- and is a per-coefficient branch on secret data when the caller is
`falcon_sign`.  It also returns a `BigInt`, so the accept/reject comparison
allocates.

Neither is needed there.  A signature candidate is a vector of `Int`s obtained
by rounding, and its coefficients cannot approach `2^40`: `s2` is bounded by
the encoder's format and `s1` by `c - s2*h` over centred residues.  So the
accumulator is provably safe *from the types and lengths alone*, which are
public, and the loop can run unconditionally.

CONSTANT TIME: this is one of the few places in this implementation where the
constant-time version was also the fast version, so there was no trade to make.
`@simd` is applied because the loop is now free of the early exit that
prevented vectorisation.  What this does NOT fix is that the *number of
rejections* in `falcon_sign` still depends on the key and the message; see
docs/constant_time.md.

Measured effect: docs/debug_log.md #051.
"""
function sqnorm_machine(v1::AbstractVector{<:Integer}, v2::AbstractVector{<:Integer})
    acc = Int128(0)
    @inbounds @simd for i in eachindex(v1)
        x = Int128(v1[i])
        acc += x * x
    end
    @inbounds @simd for i in eachindex(v2)
        x = Int128(v2[i])
        acc += x * x
    end
    return acc
end

"Largest coefficient magnitude the `Int128` accumulator is proved safe for."
const _SQNORM_LIM = Int128(1) << 40

"Total number of coefficients the same proof allows (`2^47`), checked so that
the bound is a fact about this call rather than an assumption about callers."
const _SQNORM_MAXLEN = 1 << 47

"""
Sum of squares in `Int128`, or `nothing` if the inputs are outside the range
that makes the accumulator provably safe.  Split out from [`sqnorm`](@ref) so
the inner loop is a plain typed loop over one vector.
"""
function _sqnorm_i128(vs::Tuple)
    sum(length, vs; init = 0) < _SQNORM_MAXLEN || return nothing
    acc = Int128(0)
    for v in vs
        a = _sqnorm_i128_one(acc, v)
        a === nothing && return nothing
        acc = a
    end
    return acc
end

function _sqnorm_i128_one(acc::Int128, v::AbstractVector)
    @inbounds for c in v
        # Compared in the coefficient's own type: a BigInt too wide to convert
        # fails this test rather than throwing inside the conversion.
        (-_SQNORM_LIM <= c <= _SQNORM_LIM) || return nothing
        x = Int128(c)
        acc += x * x
    end
    return acc
end

# ---------------------------------------------------------------------------
# Arithmetic modulo q
# ---------------------------------------------------------------------------
#
# Representatives are kept in [0, q).  Note that `mod` in Julia already returns
# a non-negative result for a positive modulus (unlike C's `%`), so the usual
# "add q if negative" fixup is not needed -- which is one fewer place to get it
# wrong than in the C reference.

"""
    polyaddq(f, g, q = Q)

Sum in `R_q`, with representatives in `[0, q)`.
"""
polyaddq(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}, q::Integer = Q) =
    (_checklen(f, g); Int[mod(Int(f[i]) + Int(g[i]), q) for i in eachindex(f)])

"""
    polysubq(f, g, q = Q)

Difference in `R_q`, with representatives in `[0, q)`.
"""
polysubq(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}, q::Integer = Q) =
    (_checklen(f, g); Int[mod(Int(f[i]) - Int(g[i]), q) for i in eachindex(f)])

"""
    polymulq(f, g, q = Q)

Negacyclic product in `R_q`, by the O(n^2) definition.

No widening is needed: with `q = 12289` and `n = 1024`, the accumulator is
bounded by `n * (q-1)^2 < 1.6 * 10^11`, which fits an `Int64` with a factor of
`5 * 10^7` to spare.  (Stated explicitly because the analogous bound is what
one must actually check before writing `Int32` in the C version.)
"""
function polymulq(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer},
                  q::Integer = Q)
    _checklen(f, g)
    n = length(f)
    acc = zeros(Int, n)
    @inbounds for i in 1:n, j in 1:n
        k = i + j - 2
        p = Int(f[i]) * Int(g[j])
        if k < n
            acc[k + 1] += p
        else
            acc[k - n + 1] -= p
        end
    end
    return Int[mod(c, q) for c in acc]
end

"""
    centered(f, q = Q)

Lift representatives from `[0, q)` to the *centred* range `(-q/2, q/2]`.

The distinction matters whenever a norm is taken: `||s||` of a signature must
be measured with centred representatives, since `q - 1` is a *small* number in
`R_q` (namely `-1`) but a large one in `[0, q)`.  Forgetting this is the
canonical way to make `verify` reject every honest signature -- and, because
the resulting squared norm is enormous rather than merely slightly too big, it
is at least a loud failure.
"""
centered(f::AbstractVector{<:Integer}, q::Integer = Q) =
    Int[c > q ÷ 2 ? Int(c) - Int(q) : Int(c) for c in f]

# ---------------------------------------------------------------------------
# Split and merge in the coefficient domain
# ---------------------------------------------------------------------------
#
# These are the ring-level shadow of the tower of subfields
#
#     Q(zeta_2) subset Q(zeta_4) subset ... subset Q(zeta_{2n})
#
# Writing f(x) = f0(x^2) + x * f1(x^2) splits an element of Z[x]/(x^n+1) into
# two elements of Z[x]/(x^{n/2}+1).  Every recursive algorithm in FALCON --
# the FFT, the NTRU solving descent, the ffLDL tree, ffSampling -- is a
# recursion on this decomposition.  `splitfft`/`mergefft` in fft.jl are the
# same operation seen through the FFT, and keeping the two in step is the
# whole game.
#
# [Py-ref] scripts/pyref/common.py:8-36 (`split`, `merge`)

"""
    polysplit(f) -> (f0, f1)

Split `f` of length `n` into its even and odd parts, each of length `n/2`,
such that `f(x) = f0(x^2) + x * f1(x^2)`.
"""
function polysplit(f::AbstractVector{T}) where {T}
    n = length(f)
    iseven(n) || throw(ArgumentError("cannot split an odd-length polynomial"))
    m = n ÷ 2
    return (T[f[2i - 1] for i in 1:m], T[f[2i] for i in 1:m])
end

"""
    polymerge(f0, f1) -> f

Inverse of [`polysplit`](@ref): interleave `f0` and `f1` back into `f`.
"""
function polymerge(f0::AbstractVector{T}, f1::AbstractVector{T}) where {T}
    _checklen(f0, f1)
    m = length(f0)
    f = Vector{T}(undef, 2m)
    @inbounds for i in 1:m
        f[2i - 1] = f0[i]
        f[2i] = f1[i]
    end
    return f
end
