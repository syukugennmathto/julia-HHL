# ntt.jl -- the number-theoretic transform over Z_q, q = 12289.
#
# ---------------------------------------------------------------------------
# What the NTT is, stated in the way that makes the code obvious
# ---------------------------------------------------------------------------
#
# The NTT is *multipoint evaluation*.  Over Z_q, the polynomial x^n + 1 splits
# into n distinct linear factors,
#
#     x^n + 1 = prod_{j=1}^{n} (x - w_j)      (mod q)
#
# so the Chinese Remainder Theorem gives a ring isomorphism
#
#     Z_q[x]/(x^n + 1)  ~=  Z_q x Z_q x ... x Z_q     (n copies)
#          f            |->  (f(w_1), ..., f(w_n))
#
# and *that map is the NTT*.  Multiplication becomes coordinatewise
# multiplication because evaluation is a ring homomorphism.  Everything else
# -- butterflies, twiddle factors, bit reversal -- is bookkeeping about how to
# evaluate at all n points in O(n log n) instead of O(n^2).
#
# Our tests check exactly this statement: `ntt(f)[j] == f(roots[j]) mod q`,
# evaluated by Horner.  If that holds, the transform is right, whatever the
# butterflies look like.
#
# ---------------------------------------------------------------------------
# Why a *2n*-th root of unity, not an n-th one
# ---------------------------------------------------------------------------
#
# This is the negacyclic sign of poly.jl showing up for the second time.
#
# For the *cyclic* ring Z_q[x]/(x^n - 1) the evaluation points are the n-th
# roots of unity, and an n-th root suffices.  For our ring the factorisation is
#
#     x^n + 1 = prod_{j=0}^{n-1} (x - zeta^{2j+1}),    zeta of order 2n
#
# i.e. the roots of x^n + 1 are the *odd* powers of a primitive 2n-th root of
# unity.  They are 2n-th roots that are not n-th roots.  Hence the requirement
#
#     2n | q - 1
#
# and hence q = 12289 = 3 * 2^12 + 1, which supports 2n up to 4096.
#
# ---------------------------------------------------------------------------
# Why this NTT does not look like Dilithium's
# ---------------------------------------------------------------------------
#
# A Dilithium implementation writes the NTT as an iterative Cooley-Tukey loop
# nest with a bit-reversed output order and Montgomery-form twiddles.  That is
# the fast way, and Dilithium can afford it because the NTT is the *only*
# transform it has.
#
# FALCON has two transforms -- this one over Z_q and the FFT over C -- and they
# must agree structurally, because the tower of subfields is the skeleton of
# `ffLDL` and `ffSampling`.  So the reference implementation writes both
# recursively, on the same `split`/`merge` decomposition:
#
#     f(x) = f0(x^2) + x * f1(x^2)
#
# We follow that.  It is slower than a tuned iterative NTT and that is a
# deliberate trade: the point of this module is to be *the same shape* as
# fft.jl, so that a discrepancy between them is visible as a discrepancy in
# one shared idea rather than as two unrelated bugs.
#
# NOTE ON THE C REFERENCE: the C implementation's `Zf(to_ntt_monty)` uses a
# different (bit-reversed, Montgomery-domain) coefficient order.  Cross-checking
# this module against C therefore requires a permutation and a Montgomery
# factor; cross-checking it against the Python reference does not.  That is why
# the vectors in test/vectors/ntt_kat.jl come from Python.

# ---------------------------------------------------------------------------
# The root of unity and the root tables
# ---------------------------------------------------------------------------

"""
    NTT_ORDER

`2 * 1024 = 2048`, the order of the root of unity we need to support every
degree up to `n = 1024`.
"""
const NTT_ORDER = 2048

"""
    NTT_ZETA

A primitive `NTT_ORDER`-th root of unity modulo `q`, namely `7`.

[derived] `7` is the *smallest* integer of multiplicative order exactly 2048
modulo 12289; `find_ntt_zeta()` recomputes it, and `test_ntt.jl` checks that
the search still returns this value.  It is written as a literal so that the
root tables are reproducible without running a search, and verified rather
than trusted.

Note that `zeta` itself is not a root of `x^n + 1` for any of our `n`: it has
order 2048, so `zeta^1024 = -1`, and the roots of `x^n + 1` are the odd
multiples of `1024/n` in the exponent.
"""
const NTT_ZETA = 7

"""
    find_ntt_zeta(q = Q, order = NTT_ORDER) -> Int

Search for the smallest primitive `order`-th root of unity mod `q`.

Exists so that `NTT_ZETA` is a *checked* constant rather than folklore.  Not
used at run time.
"""
function find_ntt_zeta(q::Integer = Q, order::Integer = NTT_ORDER)
    (q - 1) % order == 0 ||
        throw(ArgumentError("no order-$order root exists mod $q: $order does not divide $(q-1)"))
    for g in 2:(q - 1)
        powermod(g, order, q) == 1 || continue
        # order exactly `order`: g^(order/p) != 1 for every prime p | order.
        # Here order = 2^11, so the only prime is 2.
        powermod(g, order ÷ 2, q) == 1 && continue
        return g
    end
    throw(ErrorException("no primitive order-$order root of unity found mod $q"))
end

"""
    invmod_q(x, q = Q) -> Int

Inverse of `x` modulo the prime `q`, by Fermat's little theorem
(`x^(q-2) = x^-1`).

The reference implementation instead carries a 12289-entry lookup table
(`inv_mod_q` in `ntt_constants.py`).  A table is faster, but a table indexed
by a *secret* value is a textbook cache-timing leak, so it would have to go
in a real implementation anyway.

CONSTANT TIME: `powermod` is not constant time either (its addition chain
depends on the exponent -- though here the exponent `q-2` is public, so this
particular call is fine).  What is *not* fine, and cannot be fixed here, is
that division is only used where the divisor is secret; see `polydivq`.
"""
invmod_q(x::Integer, q::Integer = Q) =
    (m = mod(Int(x), q);
     m == 0 && throw(DivideError());
     powermod(m, q - 2, q))

# Root tables are built once per degree and cached.  They are pure functions of
# (q, NTT_ZETA, n), so the cache is a memoisation, not hidden state.
const _ROOT_CACHE = Dict{Int,Vector{Int}}()
const _ROOT_EXP_CACHE = Dict{Int,Vector{Int}}()

"""
    ntt_root_exponents(n) -> Vector{Int}

Exponents `e` such that the `j`-th evaluation point is `NTT_ZETA^e`.

## The recursion, and why the ordering is not arbitrary

At the bottom, `x^2 + 1` has the two roots `zeta^{512}` and `zeta^{1536}`
(their squares are `zeta^{1024} = -1`).

Going up one level, each root `w` of `x^{n} + 1` has exactly two square roots
`+-sqrt(w)`, and those are the two roots of `x^{2n} + 1` sitting above it.  In
exponents: `e |-> {e/2, e/2 + 1024}`.  This is precisely the statement that the
tower of subfields

    Q(zeta_4) subset Q(zeta_8) subset ... subset Q(zeta_{2n})

is a tower of quadratic extensions, read modulo q.

Within each pair the reference orders the *smaller residue first*; the pair is
always `(r, q - r)`.  That convention is what makes `split_ntt` and `merge_ntt`
below line up with `polysplit`/`polymerge`, so it is load-bearing, not
cosmetic.  `test_ntt.jl` checks the generated tables against the reference's
hard-coded ones for every `n` from 2 to 1024.
"""
function ntt_root_exponents(n::Integer)
    _check_degree(n)
    return get!(_ROOT_EXP_CACHE, Int(n)) do
        if n == 2
            return _order_pair(NTT_ORDER ÷ 4, 3 * (NTT_ORDER ÷ 4))
        end
        parent = ntt_root_exponents(n ÷ 2)
        out = Int[]
        for e in parent
            iseven(e) || throw(ErrorException(
                "odd exponent $e at degree $n: the tower cannot be descended further"))
            append!(out, _order_pair(e ÷ 2, e ÷ 2 + NTT_ORDER ÷ 2))
        end
        return out
    end
end

# Order a pair of exponents so that the smaller *residue* comes first.
function _order_pair(a::Int, b::Int)
    va = powermod(NTT_ZETA, a, Q)
    vb = powermod(NTT_ZETA, b, Q)
    return va <= vb ? [a, b] : [b, a]
end

"""
    ntt_roots(n) -> Vector{Int}

The `n` roots of `x^n + 1` modulo `q`, in the reference's recursive order.
`ntt(f)[j]` is `f` evaluated at `ntt_roots(n)[j]`.
"""
ntt_roots(n::Integer) =
    get!(_ROOT_CACHE, Int(n)) do
        [powermod(NTT_ZETA, e, Q) for e in ntt_root_exponents(n)]
    end

function _check_degree(n::Integer)
    (n >= 2 && ispow2(n)) ||
        throw(ArgumentError("degree must be a power of two and at least 2, got $n"))
    n <= NTT_ORDER ÷ 2 ||
        throw(ArgumentError("degree $n exceeds what a $(NTT_ORDER)-th root of unity supports"))
    return nothing
end

"""
    INV2_Q

The inverse of 2 modulo q, `6145`.  [Py-ref] scripts/pyref/ntt.py:15 (`i2`).
"""
const INV2_Q = 6145

# ---------------------------------------------------------------------------
# The transform
# ---------------------------------------------------------------------------

"""
    merge_ntt(f0_ntt, f1_ntt) -> f_ntt

Combine the NTTs of the even and odd halves into the NTT of the whole.

If `f(x) = f0(x^2) + x f1(x^2)` and `w` is a root of `x^n + 1`, then
`f(w) = f0(w^2) + w f1(w^2)`.  The two roots sitting above a given root `v` of
`x^{n/2} + 1` are `+-w` with `w^2 = v`, so

    f( w) = f0(v) + w f1(v)
    f(-w) = f0(v) - w f1(v)

which is the butterfly below.  The `+-` pairing is exactly why the root table
must be ordered in pairs.

[Py-ref] scripts/pyref/ntt.py:41-56 (`merge_ntt`)
"""
function merge_ntt(f0_ntt::AbstractVector{<:Integer}, f1_ntt::AbstractVector{<:Integer})
    _checklen(f0_ntt, f1_ntt)
    m = length(f0_ntt)
    n = 2m
    w = ntt_roots(n)
    out = Vector{Int}(undef, n)
    @inbounds for i in 1:m
        t = mod(w[2i - 1] * Int(f1_ntt[i]), Q)
        out[2i - 1] = mod(Int(f0_ntt[i]) + t, Q)
        out[2i]     = mod(Int(f0_ntt[i]) - t, Q)
    end
    return out
end

"""
    split_ntt(f_ntt) -> (f0_ntt, f1_ntt)

Inverse of [`merge_ntt`](@ref): recover the NTTs of the halves.

Solving the butterfly for `f0(v)` and `f1(v)`:

    f0(v) = (f(w) + f(-w)) / 2
    f1(v) = (f(w) - f(-w)) / (2w)

so this needs the inverse of 2 (that is `INV2_Q`) and the inverse of the root.

[Py-ref] scripts/pyref/ntt.py:22-39 (`split_ntt`)
"""
function split_ntt(f_ntt::AbstractVector{<:Integer})
    n = length(f_ntt)
    _check_degree(n)
    m = n ÷ 2
    w = ntt_roots(n)
    f0 = Vector{Int}(undef, m)
    f1 = Vector{Int}(undef, m)
    @inbounds for i in 1:m
        a = Int(f_ntt[2i - 1])
        b = Int(f_ntt[2i])
        f0[i] = mod(INV2_Q * (a + b), Q)
        f1[i] = mod(INV2_Q * (a - b) % Q * invmod_q(w[2i - 1]), Q)
    end
    return (f0, f1)
end

"""
    ntt(f) -> f_ntt

Forward NTT: evaluate `f` at all `n` roots of `x^n + 1`.

`ntt(f)[j] == f(ntt_roots(n)[j]) mod q`, and that identity is what
`test_ntt.jl` checks directly.

[Py-ref] scripts/pyref/ntt.py:58-75 (`ntt`)

CONSTANT TIME: the transform itself is fine -- its control flow depends only
on `n`, which is public.  `mod` on a signed value is a conditional in
disguise on some targets, which is why the C reference carries its own
`mq_add`/`mq_sub` written as branchless arithmetic.
"""
function ntt(f::AbstractVector{<:Integer})
    n = length(f)
    _check_degree(n)
    if n > 2
        f0, f1 = polysplit(f)
        return merge_ntt(ntt(f0), ntt(f1))
    else
        s = ntt_roots(2)[1]              # a square root of -1; 1479
        a, b = Int(f[1]), Int(f[2])
        return Int[mod(a + s * b, Q), mod(a - s * b, Q)]
    end
end

"""
    intt(f_ntt) -> f

Inverse NTT: interpolate from the `n` evaluations back to coefficients.

[Py-ref] scripts/pyref/ntt.py:77-95 (`intt`)
"""
function intt(f_ntt::AbstractVector{<:Integer})
    n = length(f_ntt)
    _check_degree(n)
    if n > 2
        f0n, f1n = split_ntt(f_ntt)
        return polymerge(intt(f0n), intt(f1n))
    else
        s = ntt_roots(2)[1]
        a, b = Int(f_ntt[1]), Int(f_ntt[2])
        return Int[mod(INV2_Q * (a + b), Q),
                   mod(INV2_Q * invmod_q(s) % Q * (a - b), Q)]
    end
end

# ---------------------------------------------------------------------------
# Operations in the NTT domain
# ---------------------------------------------------------------------------
#
# These are all coordinatewise, which is the entire point of the transform.

"Coordinatewise sum in the NTT domain.  [Py-ref] ntt.py:130-132"
ntt_add(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}) =
    (_checklen(f, g); Int[mod(Int(f[i]) + Int(g[i]), Q) for i in eachindex(f)])

"Coordinatewise difference in the NTT domain.  [Py-ref] ntt.py:135-137"
ntt_sub(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}) =
    (_checklen(f, g); Int[mod(Int(f[i]) - Int(g[i]), Q) for i in eachindex(f)])

"Coordinatewise product in the NTT domain.  [Py-ref] ntt.py:140-145"
ntt_mul(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}) =
    (_checklen(f, g); Int[mod(Int(f[i]) * Int(g[i]), Q) for i in eachindex(f)])

"""
    ntt_div(f_ntt, g_ntt) -> Vector{Int}

Coordinatewise quotient in the NTT domain.  Throws `DivideError` if any
coordinate of `g_ntt` is zero, i.e. if `g` is not invertible in `R_q`.

[Py-ref] scripts/pyref/ntt.py:148-155 (`div_ntt`)
"""
function ntt_div(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer})
    _checklen(f, g)
    any(iszero, g) && throw(DivideError())
    return Int[mod(Int(f[i]) * invmod_q(g[i]), Q) for i in eachindex(f)]
end

# ---------------------------------------------------------------------------
# Coefficient-domain wrappers
# ---------------------------------------------------------------------------

"""
    polymulq_ntt(f, g) -> Vector{Int}

Product in `R_q` computed through the NTT: `intt(ntt(f) .* ntt(g))`.

Must agree with the schoolbook `polymulq` from poly.jl for every input; that
agreement, between two genuinely different algorithms, is the main correctness
argument for this module.

[Py-ref] scripts/pyref/ntt.py:118-120 (`mul_zq`)
"""
polymulq_ntt(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}) =
    intt(ntt_mul(ntt(f), ntt(g)))

"""
    polydivq(f, g) -> Vector{Int}

Quotient in `R_q`.  Throws `DivideError` if `g` is not invertible.

This is how the public key is computed: `h = g / f mod q`.

[Py-ref] scripts/pyref/ntt.py:123-127 (`div_zq`)

CONSTANT TIME: this is a genuinely awkward one.  The divisor here is the
*secret* `f`, so both the inversion and the "is it invertible" test act on
secret data.  A real implementation cannot branch on the answer (the C
reference has `Zf(is_invertible)`, written branchlessly, returning a mask
rather than a bool).  We branch, and throw.
"""
polydivq(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer}) =
    intt(ntt_div(ntt(f), ntt(g)))

"""
    is_invertible_zq(f) -> Bool

Whether `f` is invertible in `R_q = Z_q[x]/(x^n+1)`.

Because `x^n + 1` splits into distinct linear factors mod q, `R_q` is a
product of `n` copies of the field `Z_q`, and an element is invertible exactly
when *no* NTT coordinate vanishes.  There is no "partially invertible" case,
which is what makes the test this cheap.

Key generation needs it: `f` must be invertible for `h = g/f` to exist, and
the reference resamples `f, g` when it is not
([Py-ref] scripts/pyref/ntrugen.py:236-238).

CONSTANT TIME: returning a `Bool` derived from secret data is exactly the leak
discussed in `polydivq`.  Noted rather than fixed; see the C reference's
`Zf(is_invertible)` (inner.h:628) for the branchless form.
"""
is_invertible_zq(f::AbstractVector{<:Integer}) = !any(iszero, ntt(f))
