# fft.jl -- the FFT over C for R[x]/(x^n + 1), and arithmetic in that domain.
#
# ---------------------------------------------------------------------------
# What this is, and why it is not FFTW
# ---------------------------------------------------------------------------
#
# Structurally this module is ntt.jl with `Z_q` replaced by `C`.  The same
# statement holds: the transform is *multipoint evaluation* at the n roots of
# x^n + 1, which over C are the 2n-th roots of unity of odd index,
#
#     w_k = exp(i * pi * k / n),   k odd
#
# and the same recursion `f(x) = f0(x^2) + x f1(x^2)` gives the same
# split/merge butterflies.
#
# We do not use FFTW.  Not because FFTW is slow -- because FFTW gives *an*
# FFT, and what FALCON needs is *this* FFT: the one whose recursion tree is
# the tower of subfields
#
#     Q(zeta_4) subset Q(zeta_8) subset ... subset Q(zeta_2n)
#
# because `ffLDL` builds a binary tree over exactly that recursion and
# `ffSampling` walks it.  `splitfft`/`mergefft` are not implementation details
# of the transform; they are operations the sampler calls directly, on
# subtrees, in the middle of sampling.  A cyclic FFT of a different radix or a
# different output order would compute the same polynomial values and still be
# useless here.
#
# ---------------------------------------------------------------------------
# Why FALCON needs floating point at all
# ---------------------------------------------------------------------------
#
# Dilithium never leaves Z_q.  FALCON has to sample from a Gaussian *on a
# lattice*, and the sampler is an LDL* (Cholesky-like) decomposition of the
# Gram matrix followed by a walk down the resulting tree.  LDL* needs square
# roots and division; those do not exist in Z_q in any useful sense, because
# "short" is a statement about the *real* embedding of the ring, not about
# residues.
#
# So FALCON computes in the complex embedding, in Float64, and the standard
# has to pin down the arithmetic well enough that everyone agrees.  That is
# the reason FN-DSA was the last of the four NIST selections to be
# standardised, and the C reference makes the stakes unusually explicit: the
# reference `config.h` *forces* the software floating-point emulation on, with
# the comment that native FPUs "may yield slight discrepancies that could
# affect determinism", and that non-determinism in signing "can lead to a
# CATASTROPHIC SECURITY FAILURE" -- two different signatures on one message
# under one key leak the key.
#
# We are not doing constant time, but this is worth stating plainly: for
# FALCON, *bit-reproducibility of Float64 across platforms is a security
# property*, not a nicety.
#
# ---------------------------------------------------------------------------
# Accuracy, measured
# ---------------------------------------------------------------------------
#
# Empirical, at n = 512, using this recursion in double precision (numbers
# measured with the Python reference, which uses the same algorithm):
#
#     coefficient size    round-trip error     product error
#     |c| <= 10           2.5e-14              5.0e-12
#     |c| <= 6144         1.7e-11              1.5e-06
#     |c| <= 1e6          2.6e-09              4.3e-02
#
# The round-trip error is ~2.5e-15 relative, i.e. a handful of ulp times
# sqrt(log n) -- textbook FFT behaviour.  The interesting column is the
# product: FALCON only ever needs the product to be correct *after rounding to
# the nearest integer*, so the tolerance is 0.5 in absolute terms.  At
# signature-sized coefficients there are five orders of magnitude of headroom;
# at coefficients of size 1e6 -- entirely realistic for the intermediate F, G
# of the NTRU solving step -- the margin is down to a factor of ten, and at
# 1e7 it is gone.
#
# That is the quantitative answer to "where do we need multiprecision": not
# "wherever the numbers look big", but specifically where the FFT product's
# error would exceed 1/2.  ntrugen.jl is where that bites.
#
# (The specification has its own analysis of this; the session that wrote this
# file could not reach the PDF -- see docs/debug_log.md #002 -- so the numbers
# above are measured rather than quoted.  They are not a substitute for the
# spec's argument and should be replaced by it when the PDF is available.)

# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------

const _FFT_ROOT_CACHE = Dict{Int,Vector{ComplexF64}}()
const _FFT_EXP_CACHE = Dict{Int,Vector{Int}}()

"""
    fft_root_exponents(n) -> Vector{Int}

Odd integers `k` with `|k| < n` such that the `j`-th evaluation point is
`exp(i pi k_j / n)`.

## The recursion

At the bottom, `x^2 + 1` has roots `+i` and `-i`, i.e. `k = +1, -1`.

Each root `w` of `x^n + 1` has two square roots `+-sqrt(w)`, and those are the
two roots of `x^{2n} + 1` above it -- the same quadratic-tower statement as in
ntt.jl, read in `C` instead of `Z_q`.  In exponents (note the denominator
doubles as we go up, so `k` itself is unchanged for the principal root):

    k  |->  { k,  k -+ n }

with the sign chosen to keep `|k| < 2n`.  The first child is the *principal*
square root (argument in `(-pi/2, pi/2]`), which is the reference's
convention.

## This is NOT the same ordering as the NTT

`ntt.jl` orders each `+-` pair by *smaller residue in [0, q)*; here we order by
*principal square root*.  Those criteria disagree at some nodes, so
`ntt_roots(n)[j]` and `fft_roots(n)[j]` are in general **not** the same root of
the tower -- one is the negation of the other.

That is harmless, because nothing in FALCON ever requires index `j` of the NTT
to mean the same thing as index `j` of the FFT: they are applied to different
objects (the NTT to the public key mod q, the FFT to the secret basis over R).
What *is* shared is the tree *shape* -- pairs, and parent = square -- and that
is what `split`/`merge` recurse on.

It is still a trap worth naming, because "the two tables are the same table
read in two rings" is the natural guess and it is false.  See
docs/debug_log.md #009.
"""
function fft_root_exponents(n::Integer)
    _check_degree(n)
    return get!(_FFT_EXP_CACHE, Int(n)) do
        n == 2 && return [1, -1]
        parent = fft_root_exponents(n ÷ 2)
        out = Int[]
        for k in parent
            push!(out, k)                        # principal square root
            push!(out, k > 0 ? k - n : k + n)    # its negation, kept in (-n, n)
        end
        return out
    end
end

"""
    fft_roots(n) -> Vector{ComplexF64}

The `n` roots of `x^n + 1` in `C`, in the reference's order.
`fft(f)[j]` is `f` evaluated at `fft_roots(n)[j]`.

Computed as `cispi(k/n)`.  Since `n` is a power of two, `k/n` is exact in
binary, so each root is correctly rounded to within an ulp or so.

Note that this is *more* accurate than the Python reference, whose
hard-coded table carries decimal constants of only about 15 significant
digits -- up to ~370 ulp of error at n = 512.  The C reference's table is
given to 27 digits and is correctly rounded.  Our tests therefore compare
against the Python vectors with a tolerance, and against the C vectors much
more tightly.  (docs/debug_log.md #010.)
"""
fft_roots(n::Integer) =
    get!(_FFT_ROOT_CACHE, Int(n)) do
        [cispi(k / n) for k in fft_root_exponents(n)]
    end

# ---------------------------------------------------------------------------
# The transform
# ---------------------------------------------------------------------------

"""
    merge_fft(f0_fft, f1_fft) -> f_fft

Combine the FFTs of the even and odd halves.  With `v = w^2`,

    f( w) = f0(v) + w f1(v)
    f(-w) = f0(v) - w f1(v)

Corresponds to algorithm `mergefft_2` of the specification.
[Py-ref] scripts/pyref/fft.py:33-51 (`merge_fft`)
"""
function merge_fft(f0_fft::AbstractVector{ComplexF64}, f1_fft::AbstractVector{ComplexF64})
    _checklen(f0_fft, f1_fft)
    m = length(f0_fft)
    n = 2m
    w = fft_roots(n)
    out = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:m
        t = w[2i - 1] * f1_fft[i]
        out[2i - 1] = f0_fft[i] + t
        out[2i]     = f0_fft[i] - t
    end
    return out
end

"""
    split_fft(f_fft) -> (f0_fft, f1_fft)

Inverse of [`merge_fft`](@ref):

    f0(v) = (f(w) + f(-w)) / 2
    f1(v) = (f(w) - f(-w)) / (2w)

Note the `/ w`: the reference multiplies by `conj(w)` instead, which is exact
here because `|w| = 1`, and is both faster and *more accurate* than a complex
division.  If you write `/ w[2i-1]` the tests still pass, at slightly worse
accuracy; if you write `* w[2i-1]` -- forgetting that the inverse rotation is
the conjugate -- everything downstream silently produces the wrong tree.  That
is the classic FALCON sign bug and it is why this line has a comment.

Corresponds to algorithm `splitfft_2` of the specification.
[Py-ref] scripts/pyref/fft.py:14-31 (`split_fft`)
"""
function split_fft(f_fft::AbstractVector{ComplexF64})
    n = length(f_fft)
    _check_degree(n)
    m = n ÷ 2
    w = fft_roots(n)
    f0 = Vector{ComplexF64}(undef, m)
    f1 = Vector{ComplexF64}(undef, m)
    @inbounds for i in 1:m
        a = f_fft[2i - 1]
        b = f_fft[2i]
        f0[i] = 0.5 * (a + b)
        # conj(w), not 1/w and emphatically not w: |w| = 1 so conj(w) = w^-1.
        f1[i] = 0.5 * (a - b) * conj(w[2i - 1])
    end
    return (f0, f1)
end

"""
    fft(f) -> Vector{ComplexF64}

Forward FFT: evaluate the real polynomial `f` at all `n` roots of `x^n + 1`.

`fft(f)[j] == f(fft_roots(n)[j])`, which is what test_fft.jl checks directly.

[Py-ref] scripts/pyref/fft.py:53-72 (`fft`)

Note the redundancy: for real `f` the values come in conjugate pairs, so half
of the output is determined by the other half.  We keep all `n` values, as the
Python reference does; the C reference keeps only `n/2` (see
docs/math/05_fft.md for the mapping).
"""
function fft(f::AbstractVector{<:Real})
    n = length(f)
    _check_degree(n)
    if n > 2
        f0, f1 = polysplit(f)
        return merge_fft(fft(f0), fft(f1))
    else
        a, b = Float64(f[1]), Float64(f[2])
        return ComplexF64[complex(a, b), complex(a, -b)]   # f(+i), f(-i)
    end
end

fft(f::AbstractVector{ComplexF64}) =
    throw(ArgumentError("fft expects coefficients (real); got a complex vector -- " *
                        "did you mean ifft?"))

"""
    ifft(f_fft) -> Vector{Float64}

Inverse FFT: interpolate back to real coefficients.

[Py-ref] scripts/pyref/fft.py:74-92 (`ifft`)

At the base case the reference reads only `f_fft[1]`, taking its real and
imaginary parts as the two coefficients, and ignores `f_fft[2]` entirely --
because for a real polynomial the second value is the conjugate of the first
and carries no new information.  We reproduce that: it is not an optimisation
we are free to change, since it is what makes `ifft` a left inverse of `fft`
*on real inputs* while quietly discarding any non-real component.
"""
function ifft(f_fft::AbstractVector{ComplexF64})
    n = length(f_fft)
    _check_degree(n)
    if n > 2
        f0f, f1f = split_fft(f_fft)
        return polymerge(ifft(f0f), ifft(f1f))
    else
        return Float64[real(f_fft[1]), imag(f_fft[1])]
    end
end

# ---------------------------------------------------------------------------
# Arithmetic in the FFT domain
# ---------------------------------------------------------------------------
#
# Coordinatewise, as always -- that is what the transform buys.  These are the
# operations ffLDL and ffSampling actually call; they spend almost all their
# time in the FFT domain and only transform back at the very end.

"Coordinatewise sum.  [Py-ref] fft.py:120-122"
add_fft(f::AbstractVector{ComplexF64}, g::AbstractVector{ComplexF64}) =
    (_checklen(f, g); ComplexF64[f[i] + g[i] for i in eachindex(f)])

"Coordinatewise difference.  [Py-ref] fft.py:125-127"
sub_fft(f::AbstractVector{ComplexF64}, g::AbstractVector{ComplexF64}) =
    (_checklen(f, g); ComplexF64[f[i] - g[i] for i in eachindex(f)])

"Coordinatewise negation."
neg_fft(f::AbstractVector{ComplexF64}) = ComplexF64[-c for c in f]

"Coordinatewise product.  [Py-ref] fft.py:130-133"
mul_fft(f::AbstractVector{ComplexF64}, g::AbstractVector{ComplexF64}) =
    (_checklen(f, g); ComplexF64[f[i] * g[i] for i in eachindex(f)])

"""
    div_fft(f_fft, g_fft)

Coordinatewise quotient.  [Py-ref] fft.py:136-141

Unlike `ntt_div`, there is no invertibility test: over `C` the only obstruction
is an exact zero, and a coordinate of `g_fft` being *near* zero is the real
danger -- it amplifies the rounding error of everything downstream without
raising anything.  In FALCON this matters in exactly one place, the
`gs_norm` computation of key generation, where dividing by `f*adj(f) + g*adj(g)`
is safe precisely because that quantity is bounded away from zero by the
rejection test.
"""
div_fft(f::AbstractVector{ComplexF64}, g::AbstractVector{ComplexF64}) =
    (_checklen(f, g); ComplexF64[f[i] / g[i] for i in eachindex(f)])

"""
    adj_fft(f_fft)

The adjoint in the FFT domain: coordinatewise complex conjugation.
[Py-ref] scripts/pyref/fft.py:144-146 (`adj_fft`)

This is the payoff of the `polyadj` discussion in poly.jl.  In the coefficient
domain the adjoint is the fiddly `adj(f)_i = -f_{n-i}` with its sign; in the
FFT domain it is just `conj`.  Multiplication by `f` is diagonal in this basis,
so its Hermitian adjoint is the conjugate diagonal -- there is nothing left to
prove.

Hence `f * adj(f)` becomes `|f_j|^2` coordinatewise: real, non-negative, and
manifestly so.  That is what makes the Gram matrix of ffsampling.jl positive
definite by construction rather than by argument.
"""
adj_fft(f::AbstractVector{ComplexF64}) = ComplexF64[conj(c) for c in f]

# ---------------------------------------------------------------------------
# Coefficient-domain wrappers
# ---------------------------------------------------------------------------

"Product in `R` via the FFT (approximate).  [Py-ref] fft.py:108-110"
polymul_fft(f::AbstractVector{<:Real}, g::AbstractVector{<:Real}) =
    ifft(mul_fft(fft(f), fft(g)))

"Quotient in `R` via the FFT (approximate).  [Py-ref] fft.py:113-115"
polydiv_fft(f::AbstractVector{<:Real}, g::AbstractVector{<:Real}) =
    ifft(div_fft(fft(f), fft(g)))

"Adjoint in `R` via the FFT.  Agrees with `polyadj` up to rounding.  [Py-ref] fft.py:118-120"
polyadj_fft(f::AbstractVector{<:Real}) = ifft(adj_fft(fft(f)))

# ---------------------------------------------------------------------------
# Interoperating with the C reference
# ---------------------------------------------------------------------------

"""
    from_c_fft(v) -> Vector{ComplexF64}

Convert the C reference's FFT representation into ours.

The C reference stores the FFT of a real polynomial of degree `n` as `n`
doubles: `n/2` real parts followed by `n/2` imaginary parts, keeping only the
independent half of the conjugate-symmetric spectrum.  Its ordering within
that half differs from ours by a **Gray-code permutation**: C index `k`
(0-based) is our index `k XOR (k >> 1)`.

This mapping was determined experimentally by building the C reference and
comparing outputs at n = 4 .. 1024; see docs/debug_log.md #011.  It is not
documented in `inner.h`, which says only "see falcon-fft.c for details on the
internal representation".

The returned vector has all `n` entries, the second half filled in by
conjugate symmetry.
"""
function from_c_fft(v::AbstractVector{<:Real})
    n = length(v)
    _check_degree(n)
    m = n ÷ 2
    out = Vector{ComplexF64}(undef, n)
    @inbounds for k in 0:(m - 1)
        j = k ⊻ (k >> 1)                     # Gray code, 0-based
        z = complex(Float64(v[k + 1]), Float64(v[m + k + 1]))
        out[j + 1] = z
        out[j + m + 1] = conj(z)             # the redundant half
    end
    return out
end

"""
    to_c_fft(f_fft) -> Vector{Float64}

Inverse of [`from_c_fft`](@ref): pack our `n`-element complex representation
into the C reference's `n`-double layout.  The redundant half of `f_fft` is
discarded (it is not checked for consistency).
"""
function to_c_fft(f_fft::AbstractVector{ComplexF64})
    n = length(f_fft)
    _check_degree(n)
    m = n ÷ 2
    out = Vector{Float64}(undef, n)
    @inbounds for k in 0:(m - 1)
        j = k ⊻ (k >> 1)
        out[k + 1] = real(f_fft[j + 1])
        out[m + k + 1] = imag(f_fft[j + 1])
    end
    return out
end

# ---------------------------------------------------------------------------
# Choosing the root table: accuracy versus reference-compatibility
# ---------------------------------------------------------------------------
#
# By default `fft_roots` computes its roots with `cispi`, which is correctly
# rounded to within an ulp and agrees with the *C* reference, whose table is
# given to 27 decimal digits.  The Python reference's table carries only ~15
# significant digits and is the outlier (docs/debug_log.md #010).
#
# For most of this project that difference is invisible: it is 1e-16 relative,
# and every FFT test passes against either table.
#
# It stops being invisible in `ffsampling` (module 8).  There the FFT result
# becomes the *centre* of a discrete Gaussian, the sampler's rejection loop
# compares against random bytes, and a last-ulp difference in the centre flips
# a comparison and returns a different integer.  Measured on the six recorded
# ffSampling vectors: with our roots, one of six reproduces the Python
# reference; with the Python table installed, six of six.
#
# That is the whole FN-DSA floating-point problem in one measurement, and it is
# why the C reference emulates floating point in integer arithmetic instead of
# trusting the FPU.  See docs/math/08_ffsampling.md.
#
# So the table is made switchable -- not to paper over the difference, but so
# that the difference can be *demonstrated*.  The inaccurate table is not
# shipped in `src/`; it lives in `test/vectors/fft_kat.jl` as what it is,
# recorded reference data.

"""
    set_fft_roots!(table)

Install an explicit root table, replacing the computed one.  `table` maps a
degree `n` to the `n` roots of `x^n + 1` in that degree's expected order.

Intended for reproducing another implementation's floating-point results
exactly; see the discussion above.  Call [`reset_fft_roots!`](@ref) to go back
to the computed roots.
"""
function set_fft_roots!(table::AbstractDict{Int,Vector{ComplexF64}})
    empty!(_FFT_ROOT_CACHE)
    for (n, v) in table
        _check_degree(n)
        length(v) == n || throw(ArgumentError(
            "root table for degree $n has $(length(v)) entries"))
        _FFT_ROOT_CACHE[n] = copy(v)
    end
    return nothing
end

"""
    reset_fft_roots!()

Discard any installed root table and go back to computing roots with `cispi`.
"""
reset_fft_roots!() = (empty!(_FFT_ROOT_CACHE); nothing)

"""
    with_fft_roots(f, table)

Run `f()` with `table` installed as the root table, restoring the previous
state afterwards even if `f` throws.
"""
function with_fft_roots(f, table::AbstractDict{Int,Vector{ComplexF64}})
    saved = copy(_FFT_ROOT_CACHE)
    try
        set_fft_roots!(table)
        return f()
    finally
        empty!(_FFT_ROOT_CACHE)
        merge!(_FFT_ROOT_CACHE, saved)
    end
end
