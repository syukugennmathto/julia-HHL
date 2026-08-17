# params.jl -- FALCON / FN-DSA parameter sets.
#
# PROVENANCE POLICY
# -----------------
# No constant in this file was written from memory.  Every value carries a tag
# saying where it was read:
#
#   [C-ref]  the FALCON reference implementation in C, file:line
#            (https://github.com/algorand/falcon, a mirror of the round-3
#            reference code; the macros there are the normative sizes)
#   [Py-ref] the Python reference implementation, file:line
#            (https://github.com/tprest/falcon.py, MIT, Thomas Prest;
#            vendored under scripts/pyref/)
#   [derived] reproduced from a closed-form formula, and checked numerically
#            against the two references above by test/test_params.jl
#
#   [SPEC?]  <-- READ THIS.  The session that wrote this file had no network
#            route to falcon-sign.info, nvlpubs.nist.gov or eprint.iacr.org,
#            so the specification PDF itself could not be opened and the
#            *table numbers* could not be cited.  Every value below is
#            nevertheless pinned by two mutually independent implementations
#            (C and Python) that agree exactly.  The `spec_ref` fields are
#            deliberately left as "TODO" rather than guessed: filling them in
#            requires the PDF, and an invented table number is worse than an
#            absent one.

"""
    FalconParams

One FALCON parameter set.  Field-by-field provenance is given in the
constructor calls below (`FALCON_512`, `FALCON_1024`).

Distinguish carefully between the three standard deviations, because they
belong to three different Gaussians and mixing them up is a classic bug:

  * `sigma_fg`  -- used at *key generation* time, to sample the short
                   polynomials `f` and `g` coefficient-wise over Z.
  * `sigma`     -- the standard deviation of the *signature* Gaussian over the
                   lattice; it is the target width of `ffSampling`.
  * `sigma_min`
    /`sigma_max`-- the admissible range for the per-node standard deviation
                   that `ffSampling` hands to the integer sampler `samplerz`.
                   `sigma_min` is a smoothing parameter of Z; `sigma_max` is
                   the width of the fixed half-Gaussian that the CDT base
                   sampler implements.
"""
struct FalconParams
    "ring degree n; the ring is R = Z[x]/(x^n + 1)"
    n::Int
    "log2(n); the reference code indexes all its tables by this"
    logn::Int
    "modulus q of the NTRU lattice"
    q::Int

    "std. dev. of the signature Gaussian over the lattice"
    sigma::Float64
    "lower bound on the per-node std. dev. passed to samplerz"
    sigma_min::Float64
    "upper bound on the per-node std. dev.; width of the CDT base sampler"
    sigma_max::Float64
    "effective std. dev. of the key-generation polynomials f, g over Z"
    sigma_fg::Float64

    "acceptance bound beta^2 for ||(s1,s2)||^2; *inclusive*, equals floor(beta^2)"
    sig_bound::Int

    "public key length in bytes (exact)"
    pubkey_bytes::Int
    "private key length in bytes (exact)"
    privkey_bytes::Int
    "signature length in bytes in the PADDED format (exact)"
    sig_bytes::Int

    "placeholder for the specification table reference; see PROVENANCE POLICY"
    spec_ref::String
end

# ---------------------------------------------------------------------------
# The modulus
# ---------------------------------------------------------------------------

"""
    Q

The FALCON modulus, `q = 12289 = 12*1024 + 1 = 3*2^12 + 1`.

[Py-ref] scripts/pyref/common.py:5  (`q = 12 * 1024 + 1`)

The shape `3*2^12 + 1` is what makes the NTT possible: `q - 1 = 3*2^12` is
divisible by `2^12`, so Z_q contains a primitive `2^12`-th root of unity, hence
a primitive `2n`-th root for every `n <= 2^11`.  A primitive *2n*-th root
(rather than an n-th root) is exactly what the negacyclic ring `x^n + 1`
requires -- see the discussion in ntt.jl.

Contrast with Dilithium/ML-DSA, where `q = 8380417 = 2^23 - 2^13 + 1` is far
larger.  FALCON can afford a 14-bit modulus because its security comes from
the NTRU lattice being *short*, not from a large modulus-to-noise ratio.
"""
const Q = 12289

# Sanity: the two shapes must agree.  (Cheap, and catches a fat-fingered digit.)
@assert Q == 12 * 1024 + 1
@assert Q == 3 * 2^12 + 1

# ---------------------------------------------------------------------------
# Lengths that do not depend on the degree
# ---------------------------------------------------------------------------

"Length in bytes of the header byte prepended to a signature.  [Py-ref] falcon.py:41"
const HEAD_LEN = 1

"Length in bytes of the signing salt r.  [Py-ref] falcon.py:42"
const SALT_LEN = 40

"Length in bytes of the seed of the ChaCha20 PRNG.  [Py-ref] falcon.py:43"
const SEED_LEN = 56

# ---------------------------------------------------------------------------
# Closed forms (checked against the references in test/test_params.jl)
# ---------------------------------------------------------------------------

"""
    smoothing_eta(eps) -> Float64

The smoothing parameter of the integer lattice Z at parameter `eps`:

    eta_eps(Z) = (1/pi) * sqrt( (1/2) * ln(2 * (1 + 1/eps)) )

[derived]  Reproduces the reference `sigmin` values to a relative error of
about 3e-13 when evaluated at `eps = falcon_eps(n)` (see below), which is the
agreement one expects between a Float64 evaluation here and however the
reference authors evaluated it.  The *authoritative* values remain the
literals in `FALCON_512` / `FALCON_1024`; this function exists to document
where they come from and to catch a mistyped digit.
"""
smoothing_eta(eps::Float64) = (1 / pi) * sqrt(0.5 * log(2 * (1 + 1 / eps)))

"""
    falcon_eps(n) -> Float64

The statistical-distance budget `eps = 1 / sqrt(2^64 * n^3)`.

[derived]  Recovered by inverting the reference `sigmin` values: solving
`smoothing_eta(eps) == sigmin` gives `log2(1/eps) = 45.5` for `n = 512` and
`47.0` for `n = 1024`, i.e. `1/eps = 2^32 * n^(3/2)`.  The `2^64` is the usual
bound on the attacker's work and `n^3` the usual union bound over queries; the
exact wording belongs in [SPEC?] and is left for the PDF.
"""
falcon_eps(n::Integer) = 1 / sqrt(2.0^64 * float(n)^3)

"""
    SIGMA_FG_BASE

The *fixed* standard deviation `1.43300980528773` at which key generation
samples its raw Gaussian coefficients, independently of `n`.

[Py-ref] scripts/pyref/ntrugen.py:210-211
         (comment `# 1.17 * sqrt(12289 / 8192)`, then `sigma = 1.43300980528773`)

This constant is easy to misread, so it is worth stating precisely what the
reference does (`gen_poly`, ntrugen.py:204-217):

  * it draws **4096** samples from `D_{Z, 0, SIGMA_FG_BASE}`, whatever `n` is;
  * it then sums them in consecutive blocks of `k = 4096 / n`.

Summing `k` independent Gaussians of width `s` gives a Gaussian of width
`sqrt(k) * s`, so the effective width is

    sqrt(4096/n) * 1.17 * sqrt(q/8192) = 1.17 * sqrt(q / (2n))

which is the `sigma_fg` field of `FalconParams`.  Reproducing this *exact*
procedure -- 4096 draws, then folding -- matters for KAT agreement: sampling
`n` coefficients directly at `sigma_fg` gives the same distribution but a
different byte-for-byte output from the same PRNG stream.
"""
const SIGMA_FG_BASE = 1.43300980528773

"""
    gram_schmidt_quality() -> Float64

The quality factor `1.17`: key generation is rejected unless the Gram-Schmidt
norm of the NTRU basis satisfies `||B~||^2 <= 1.17^2 * q`.

[Py-ref] scripts/pyref/ntrugen.py:232  (`if gs_norm(f, g, q) > (1.17 ** 2) * q`)

This single constant is what ties the three sigmas together: the signature
width is `sigma = 1.17 * sqrt(q) * sigma_min`, i.e. (Gram-Schmidt norm bound)
times (smoothing parameter of Z).  That is the standard "Klein/GPV sampler is
safe at width ||B~|| * eta_eps(Z)" statement, specialised to FALCON's basis.
"""
gram_schmidt_quality() = 1.17

# ---------------------------------------------------------------------------
# The parameter sets
# ---------------------------------------------------------------------------

"""
    FALCON_512

FALCON-512 / FN-DSA-512.

Byte lengths cross-checked against the C reference macros, evaluated at
`logn = 9`:

    FALCON_PUBKEY_SIZE(9)  = (7 << 7) + 1                = 897
    FALCON_PRIVKEY_SIZE(9) = ((10 - 4) << 7) + 512 + 1   = 1281
    FALCON_SIG_PADDED_SIZE(9)                            = 666

[C-ref] falcon.h:308 (PRIVKEY), falcon.h:317 (PUBKEY), falcon.h:335 (SIG_PADDED)
[C-ref] common.c:241-253 (`l2bound[]`, indexed by logn; `l2bound[9] = 34034726`)
[Py-ref] falcon.py:129-135 (`params[512]`: sigma, sigmin, sig_bound, sig_bytelen)

The two references agree exactly on `sig_bound = 34034726` and
`sig_bytes = 666`, which is the cross-check that matters most: those are the
two constants a wrong value of which produces a "verification always fails"
or "signature never fits" bug that is miserable to bisect.
"""
const FALCON_512 = FalconParams(
    512,                    # n
    9,                      # logn
    Q,                      # q
    165.7366171829776,      # sigma       [Py-ref] falcon.py:131
    1.2778336969128337,     # sigma_min   [Py-ref] falcon.py:132
    1.8205,                 # sigma_max   [Py-ref] samplerz.py:10 (MAX_SIGMA)
    1.17 * sqrt(Q / (2 * 512)),  # sigma_fg  [Py-ref] ntrugen.py:208 (effective width)
    34034726,               # sig_bound   [C-ref] common.c:250 == [Py-ref] falcon.py:133
    897,                    # pubkey_bytes   [C-ref] falcon.h:317 @ logn=9
    1281,                   # privkey_bytes  [C-ref] falcon.h:308 @ logn=9
    666,                    # sig_bytes   [C-ref] falcon.h:335 == [Py-ref] falcon.py:134
    "TODO: spec table number (PDF unreachable from the authoring session)",
)

"""
    FALCON_1024

FALCON-1024 / FN-DSA-1024.  Same provenance as `FALCON_512`, at `logn = 10`:

    FALCON_PUBKEY_SIZE(10)  = (7 << 8) + 1                = 1793
    FALCON_PRIVKEY_SIZE(10) = ((10 - 5) << 8) + 1024 + 1  = 2305
    FALCON_SIG_PADDED_SIZE(10)                            = 1280

[C-ref] common.c:252 (`l2bound[10] = 70265242`)
[Py-ref] falcon.py:137-143 (`params[1024]`)
"""
const FALCON_1024 = FalconParams(
    1024,                   # n
    10,                     # logn
    Q,                      # q
    168.38857144654395,     # sigma       [Py-ref] falcon.py:139
    1.298280334344292,      # sigma_min   [Py-ref] falcon.py:140
    1.8205,                 # sigma_max   [Py-ref] samplerz.py:10 (MAX_SIGMA)
    1.17 * sqrt(Q / (2 * 1024)),  # sigma_fg [Py-ref] ntrugen.py:208 (effective width)
    70265242,               # sig_bound   [C-ref] common.c:252 == [Py-ref] falcon.py:141
    1793,                   # pubkey_bytes   [C-ref] falcon.h:317 @ logn=10
    2305,                   # privkey_bytes  [C-ref] falcon.h:308 @ logn=10
    1280,                   # sig_bytes   [C-ref] falcon.h:335 == [Py-ref] falcon.py:141
    "TODO: spec table number (PDF unreachable from the authoring session)",
)

"""
    params(n) -> FalconParams

Look up the parameter set for degree `n`.  Only the two standardised degrees
are provided; the reference implementations also define toy degrees 2..256,
which we will add if and when the test suite needs them (they are useful for
debugging `ffsampling`, where n = 512 is far too big to read by eye).
"""
function params(n::Integer)
    n == 512 && return FALCON_512
    n == 1024 && return FALCON_1024
    throw(ArgumentError("unsupported FALCON degree n = $n (expected 512 or 1024)"))
end

# ---------------------------------------------------------------------------
# The C reference size macros, transcribed
# ---------------------------------------------------------------------------
# These are kept as executable Julia so the test suite can check the literals
# above against the macro at *every* logn, not just at 9 and 10.  A transcription
# error in a shift then shows up immediately instead of at KAT time.

"[C-ref] falcon.h:317, FALCON_PUBKEY_SIZE"
c_pubkey_size(logn::Integer) = (logn <= 1 ? 4 : (7 << (logn - 2))) + 1

"[C-ref] falcon.h:308, FALCON_PRIVKEY_SIZE"
c_privkey_size(logn::Integer) =
    (logn <= 3 ? (3 << logn) : ((10 - (logn >> 1)) << (logn - 2)) + (1 << logn)) + 1

"[C-ref] falcon.h:335, FALCON_SIG_PADDED_SIZE"
c_sig_padded_size(logn::Integer) =
    44 + 3 * (256 >> (10 - logn)) + 2 * (128 >> (10 - logn)) +
    3 * (64 >> (10 - logn)) + 2 * (16 >> (10 - logn)) -
    2 * (2 >> (10 - logn)) - 8 * (1 >> (10 - logn))

"[C-ref] falcon.h:327, FALCON_SIG_COMPRESSED_MAXSIZE"
c_sig_compressed_maxsize(logn::Integer) =
    ((((11 << logn) + (101 >> (10 - logn))) + 7) >> 3) + 41
