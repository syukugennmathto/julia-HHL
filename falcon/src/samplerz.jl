# samplerz.jl -- the discrete Gaussian sampler over Z.
#
# ===========================================================================
# This is the module the whole book is about
# ===========================================================================
#
# Everything so far has been machinery.  This is the part that made FN-DSA the
# last of the NIST selections to be standardised, and the reason is visible in
# the code below rather than in any theorem.
#
# The job: given a centre `mu` and a width `sigma`, return an integer `z`
# distributed as D_{Z, mu, sigma}.  Two properties are required *at once*:
#
#   1. the output distribution must be correct to within a statistical distance
#      of about 2^-45 (that is what `sigma_min` was derived from, params.jl);
#   2. the running time, the memory access pattern, and the number of random
#      bytes consumed must not depend on `mu`, on `sigma`, or on `z`.
#
# Requirement 2 is where it hurts.  The natural way to sample a Gaussian is
# rejection sampling, and rejection sampling *by construction* runs for a
# secret-dependent number of iterations.  And `mu` here is not public: it comes
# from the secret basis, one coordinate at a time, inside `ffSampling`.  A
# timing signal that leaks which `z` was returned leaks the lattice point,
# which leaks the basis.
#
# ---------------------------------------------------------------------------
# How the reference squares that circle
# ---------------------------------------------------------------------------
#
# Three tricks, stacked:
#
#   * **A fixed base distribution.**  Instead of sampling at width `sigma`,
#     sample from a *fixed* half-Gaussian of width `MAX_SIGMA = 1.8205` using a
#     table lookup that touches every table entry every time (`basesampler`).
#     `sigma` never enters the table, so the table access is independent of it.
#     This is why params.jl has a `sigma_max` that is a strange decimal: it is
#     the width the table was built for, not a derived quantity.
#
#   * **Rejection against the *ratio* of the two Gaussians.**  The acceptance
#     probability is `ccs * exp(-x)` where `x` is the log-ratio.  Because the
#     base distribution is fixed and slightly *wider* than any admissible
#     `sigma`, the ratio is always a valid probability -- this is what the
#     constraint `sigma_min < sigma < sigma_max` buys.
#
#   * **A polynomial approximation to exp, evaluated in integer arithmetic.**
#     `approxexp` computes `2^63 * ccs * exp(-x)` with a degree-12 polynomial
#     in fixed point, entirely in 64-bit integers.  No `exp()` call, no
#     floating-point branch, no lookup indexed by anything secret.  The
#     coefficients are lifted from FACCT (doi:10.1109/TC.2019.2940949).
#
# What is *not* fixed: the number of rejection iterations still depends on the
# samples drawn.  The reference accepts that leak and mitigates it elsewhere;
# a fully constant-time implementation has to do more work still.  We do not
# attempt it -- see the CONSTANT TIME notes below for what a real one must do.
#
# ---------------------------------------------------------------------------
# Why this is the module without a safety net
# ---------------------------------------------------------------------------
#
# ntrugen.jl also used floating point (Babai reduction), and it was safe
# because the computation *iterates*: an inaccurate `k` costs a slower
# reduction, not a wrong answer (docs/math/06_ntrugen.md 6.5).
#
# Here there is no correcting loop.  The Float64 value `x` is fed straight into
# `berexp`, whose single bit decides whether this `z` is the output.  A
# different rounding of `x` gives a different bit gives a different signature.
# So *this* is the place where "the same message under the same key must give
# the same signature on every platform" turns into a constraint on floating
# point -- and non-determinism here is the catastrophic failure that the C
# reference's config.h warns about (docs/math/05_fft.md 5.7).
#
# CONSTANT TIME (whole module): out of scope here, and comprehensively so.
# What a real implementation must change:
#   * `basesampler` must not early-exit its table scan (ours does not, but it
#     also does not force the comparisons to be branchless);
#   * `berexp`'s byte-comparison loop breaks as soon as the bytes differ,
#     which leaks;
#   * `samplerz`'s rejection loop count is data-dependent;
#   * every `Float64` operation on `x` must be shown to round identically on
#     every target, which is why the C reference emulates floating point in
#     integer arithmetic rather than trusting the FPU.

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

"""
    RCDT_PREC

Precision of the reverse cumulative distribution table, in bits: 72.
[Py-ref] scripts/pyref/samplerz.py:13
"""
const RCDT_PREC = 72

"""
    INV_2SIGMA2

`1 / (2 * MAX_SIGMA^2)`, the exponent scale of the fixed base distribution.
[Py-ref] scripts/pyref/samplerz.py:11
"""
const INV_2SIGMA2 = 1 / (2 * (1.8205^2))

"""
    LN2, ILN2

`ln(2)` and `1/ln(2)`, **as the reference spells them** -- truncated decimals,
not the correctly-rounded Float64 values.

    LN2  = 0.69314718056     vs  log(2)  = 0.6931471805599453
    ILN2 = 1.44269504089     vs  1/log(2) = 1.4426950408889634

The difference is about 1e-12 relative.  It does not matter for the
*distribution*, because `berexp` only uses these to split `x` into an integer
number of halvings plus a remainder, and any consistent split works.  It
matters enormously for **reproducing the reference's byte stream**, because a
different split gives a different `s`, a different shift, and eventually a
different accept/reject decision on the same random bytes.

So these are transcribed exactly and not "corrected".  This is the clearest
small example in the project of an implementation detail that has become
normative.

[Py-ref] scripts/pyref/samplerz.py:16-17
"""
const LN2 = 0.69314718056
const ILN2 = 1.44269504089

"""
    RCDT

Reverse cumulative distribution table of a distribution very close to the
half-Gaussian `D_{Z+, 0, MAX_SIGMA}`, at 72 bits of precision.

`RCDT[i]` is `2^72 * P(z0 >= i)` for `i = 1 .. 18`, so `basesampler` can draw a
72-bit uniform `u` and count how many entries exceed it.

Note which half-Gaussian: `D_{Z+, 0, sigma}` weights **every** non-negative
integer once, `P(z) ∝ exp(-z^2/(2 sigma^2))` for `z = 0, 1, 2, ...`.  It is not
the "folded" distribution that counts `z > 0` twice; that would give
`P(z0 = 0) = 0.219` where the table says `0.359`.  The sign is supplied
separately by `samplerz`, which is why no folding happens here.

Eighteen entries of ~72 bits: these do **not** fit in `UInt64`, and Python's
arbitrary-precision integers hide that fact completely.  `UInt128` here.

[Py-ref] scripts/pyref/samplerz.py:24-43

The values are transcribed, not derived -- but `test_samplerz.jl` checks them
against the half-Gaussian they are supposed to represent, which is what would
catch a transposed digit.
"""
const RCDT = UInt128[
    3024686241123004913666,
    1564742784480091954050,
    636254429462080897535,
    199560484645026482916,
    47667343854657281903,
    8595902006365044063,
    1163297957344668388,
    117656387352093658,
    8867391802663976,
    496969357462633,
    20680885154299,
    638331848991,
    14602316184,
    247426747,
    3104126,
    28824,
    198,
    1,
]

"""
    EXP_COEFFS

Coefficients of the degree-12 polynomial approximating `exp(-x)` in fixed
point, from FACCT (doi:10.1109/TC.2019.2940949):

    (2^-63) * sum(EXP_COEFFS[13 - i] * x^i)  ~=  exp(-x)

[Py-ref] scripts/pyref/samplerz.py:48-62 (the list `C`)

`UInt64`, not `Int64`: the last coefficient is `0x8000000000000000 = 2^63`,
which does not fit in a signed 64-bit integer.  In Python this is invisible.
"""
const EXP_COEFFS = UInt64[
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
]

# ---------------------------------------------------------------------------
# The base sampler
# ---------------------------------------------------------------------------

"""
    basesampler(randombytes) -> Int

Sample `z0` in `{0, ..., 18}` from a distribution very close to the
half-Gaussian `D_{Z+, 0, 1.8205}`.

Consumes exactly `RCDT_PREC/8 = 9` bytes, always, whatever it returns.  It
draws a 72-bit uniform and counts how many table entries it falls below; the
count *is* the sample, because the table is a reverse CDF.

[Py-ref] scripts/pyref/samplerz.py:65-75

CONSTANT TIME: the loop visits all 18 entries unconditionally, which is the
right shape.  It is still not constant time in Julia: `<` on `UInt128` and the
`+` of a `Bool` are not guaranteed branchless. The C reference writes this
comparison as arithmetic on the sign bit.
"""
function basesampler(randombytes)
    bytes = randombytes(RCDT_PREC >> 3)
    length(bytes) == RCDT_PREC >> 3 ||
        throw(ArgumentError("basesampler needs $(RCDT_PREC >> 3) bytes, got $(length(bytes))"))
    u = UInt128(0)
    for i in 1:(RCDT_PREC >> 3)
        u |= UInt128(bytes[i]) << (8 * (i - 1))       # little-endian
    end
    z0 = 0
    for elt in RCDT
        z0 += Int(u < elt)
    end
    return z0
end

# ---------------------------------------------------------------------------
# exp, in integers
# ---------------------------------------------------------------------------

"""
    approxexp(x, ccs) -> UInt64

An integral approximation of **`2^64 * ccs * exp(-x)`**.  Both `x` and `ccs`
must be positive; `x` is expected in `[0, ln 2)` and `ccs` in `(0, 1)`.

[Py-ref] scripts/pyref/samplerz.py:78-98

## The scale is 2^64, not 2^63

The reference's own docstring says "an integral approximation of
`2^63 * ccs * exp(-x)`".  That is off by a factor of two: the line
`z = int(ccs * (1 << 63)) << 1` doubles it, and measurement confirms the ratio
is exactly 2.000000 across the domain.

`2^64` is the *correct* scale for what the value is used for, so this is a slip
in the comment rather than in the code: `berexp` compares the result against a
uniform 64-bit integer, eight bytes at a time, so the probability has to be
expressed as a fraction of `2^64`.  Recorded because copying the reference's
claim into our own docstring is exactly how a wrong statement propagates --
see docs/debug_log.md #023.

## Widths

This is a Horner evaluation in fixed point, and the intermediate products are
what Python hides:

  * `y` and the coefficients live in `UInt64` (the top coefficient is exactly
    `2^63`, so `Int64` is already wrong);
  * `z * y` reaches `2^126`, so the product must be formed in `UInt128` before
    the `>> 63`;
  * the final `z = 2 * floor(ccs * 2^63)` reaches just under `2^64`, and its
    product with `y` reaches just under `2^127` -- again `UInt128`.

Transliterating this from Python without choosing those types deliberately
gives silent wraparound, and the resulting sampler still produces
plausible-looking integers.  Only a KAT catches it.
"""
function approxexp(x::Float64, ccs::Float64)
    y = EXP_COEFFS[1]
    # `int(x * 2^63)` in Python truncates toward zero; x >= 0 here so trunc
    # and floor agree.  x < ln2 keeps this below 2^63.
    z = UInt64(trunc(UInt64, x * 9223372036854775808.0))   # 2^63
    @inbounds for i in 2:length(EXP_COEFFS)
        y = EXP_COEFFS[i] - UInt64((UInt128(z) * UInt128(y)) >> 63)
    end
    z = UInt64(trunc(UInt64, ccs * 9223372036854775808.0)) << 1
    y = UInt64((UInt128(z) * UInt128(y)) >> 63)
    return y
end

"""
    berexp(x, ccs, randombytes) -> Bool

Return `true` with probability approximately `ccs * exp(-x)`.

The trick: `exp(-x)` for large `x` is out of range for the polynomial, so `x`
is split as `x = s*ln2 + r` with `r in [0, ln2)`; then
`exp(-x) = exp(-r) / 2^s`, and dividing by `2^s` is a shift.

The acceptance test then compares the 64-bit value `z` against a uniform
64-bit value, one byte at a time from the top, stopping at the first byte that
differs.  Comparing from the top is what makes the *expected* number of random
bytes small (a little over one).

[Py-ref] scripts/pyref/samplerz.py:101-118

CONSTANT TIME: the early `break` is a direct timing leak of how many leading
bytes matched -- which correlates with the acceptance probability, hence with
`x`, hence with the secret centre.  A real implementation runs all eight
iterations and combines the results branchlessly.  The C reference does
exactly that.
"""
function berexp(x::Float64, ccs::Float64, randombytes)
    s = trunc(Int, x * ILN2)          # Python int(): truncation toward zero
    r = x - s * LN2
    s = min(s, 63)
    z = (approxexp(r, ccs) - UInt64(1)) >> s
    w = 0
    for i in 56:-8:0
        p = randombytes(1)
        length(p) == 1 || throw(ArgumentError("berexp needs 1 byte per iteration"))
        w = Int(p[1]) - Int((z >> i) & 0xFF)
        w != 0 && break
    end
    return w < 0
end

# ---------------------------------------------------------------------------
# The sampler
# ---------------------------------------------------------------------------

"""
    samplerz(mu, sigma, sigmin, randombytes) -> Int

Sample an integer from the discrete Gaussian `D_{Z, mu, sigma}`.

The inputs must satisfy `1 < sigmin < sigma < 1.8205`; `ffSampling` guarantees
this by construction, and `params.jl` explains where the bounds come from.

`randombytes(k)` must return `k` pseudorandom bytes.  The order and *number* of
bytes consumed is part of the specification: 9 for the base sampler, 1 for the
sign, then 1 to 8 inside `berexp`, per rejection round.  Reproducing a KAT
means reproducing that consumption exactly.

[Py-ref] scripts/pyref/samplerz.py:121-155

Corresponds to `SamplerZ` of the specification.
"""
function samplerz(mu::Float64, sigma::Float64, sigmin::Float64, randombytes)
    s = Int(floor(mu))
    r = mu - s
    dss = 1 / (2 * sigma * sigma)
    ccs = sigmin / sigma

    while true
        z0 = basesampler(randombytes)
        b = Int(randombytes(1)[1]) & 1
        z = b + (2b - 1) * z0
        x = ((z - r)^2) * dss
        x -= (z0^2) * INV_2SIGMA2
        if berexp(x, ccs, randombytes)
            return z + s
        end
    end
end

samplerz(mu::Real, sigma::Real, sigmin::Real, randombytes) =
    samplerz(Float64(mu), Float64(sigma), Float64(sigmin), randombytes)

# ---------------------------------------------------------------------------
# Byte sources
# ---------------------------------------------------------------------------

"""
    bytesource(rng::ChaCha20) -> Function

Adapt the reference PRNG of shake.jl into the `randombytes` callable that
`samplerz` expects.
"""
bytesource(rng::ChaCha20) = k -> randombytes!(rng, k)

"""
    ReplayBytes(bytes; reversed_chunks = false)

A `randombytes` source that hands out a fixed byte string, for KAT replay.
Throws if more bytes are requested than were supplied -- which is itself
informative: it means the sampler consumed more randomness than the reference
did on the same input, i.e. it took a different path.

## `reversed_chunks`

The official `samplerz` KAT vectors must be replayed with `reversed_chunks =
true`.  Their harness (`test.py`, `KAT_randbytes`) takes the next `k` bytes
from the recorded string and returns them **reversed**:

    oc = octets[:2k];  octets = octets[2k:];  return bytes.fromhex(oc)[::-1]

That reversal is real, unlike the pair of reversals in the reference PRNG which
cancel exactly (docs/debug_log.md #005).  Feeding the recorded bytes in
their natural order gives a different sample -- not an error, just a different
draw, so nothing complains.  See docs/debug_log.md #022.
"""
mutable struct ReplayBytes
    bytes::Vector{UInt8}
    pos::Int
    reversed_chunks::Bool
end

ReplayBytes(bytes::AbstractVector{UInt8}; reversed_chunks::Bool = false) =
    ReplayBytes(collect(bytes), 0, reversed_chunks)

function (rb::ReplayBytes)(k::Integer)
    rb.pos + k <= length(rb.bytes) || throw(ArgumentError(
        "replay exhausted: asked for $k more bytes at offset $(rb.pos), " *
        "only $(length(rb.bytes)) available -- the sampler took a different " *
        "path from the reference"))
    out = rb.bytes[(rb.pos + 1):(rb.pos + k)]
    rb.pos += k
    return rb.reversed_chunks ? reverse(out) : out
end

"Number of bytes consumed so far."
consumed(rb::ReplayBytes) = rb.pos
