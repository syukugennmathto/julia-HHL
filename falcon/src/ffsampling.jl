# ffsampling.jl -- LDL* over the ring, the Falcon tree, and fast Fourier sampling.
#
# ===========================================================================
# What this module is for
# ===========================================================================
#
# Signing has to solve this problem: given a target point `c` (the hash of the
# message) and the secret short basis `B`, find a lattice point close to `c` --
# close enough that `||c - v||` lands under `beta^2`, and *randomised* so that
# the signature leaks nothing about `B`.
#
# The classical answer is Klein/GPV sampling: walk the Gram-Schmidt
# orthogonalisation of `B` from the last coordinate to the first, sampling one
# integer per coordinate from a Gaussian centred on the current residual.  That
# is O(n^2) per signature and needs the full 2n x 2n Gram-Schmidt basis, which
# for n = 512 is a megabyte of doubles.
#
# FALCON's contribution -- the "F" in the name -- is doing the same thing in
# O(n log n) time and O(n) space, by exploiting the tower of subfields one last
# time.  Instead of orthogonalising a 2n x 2n real matrix, orthogonalise a
# 2 x 2 matrix *over the ring*, then recurse into the two diagonal entries by
# splitting them down the tower.  The result is a binary tree with n leaves,
# and sampling is a walk down it.
#
# So: `ffLDL` is Gram-Schmidt over the ring, and `ffSampling` is Klein's
# algorithm over the ring.  Everything that made the previous seven modules
# necessary -- `polysplit`, `splitfft`, the adjoint, the FFT that had to have
# *this* recursion tree -- exists so that these two functions can be written.
#
# ---------------------------------------------------------------------------
# Why LDL* and not Cholesky
# ---------------------------------------------------------------------------
#
# The Gram matrix `G = B * adj(B)^T` is Hermitian and positive definite, so it
# has a Cholesky factorisation `G = C * adj(C)^T`.  But Cholesky needs square
# roots, and `sqrt` of a ring element is not a ring element.
#
# LDL* factors instead as `G = L * D * adj(L)^T` with `L` unit lower triangular
# and `D` diagonal.  No square roots: `L10 = G10 / G00` and
# `D11 = G11 - L10 adj(L10) G00` are ring operations.  The square roots appear
# only at the very end, on the *leaves*, where the entries have become plain
# positive reals -- that is what `normalize_tree!` does, turning `||b_i||^2`
# into `sigma / ||b_i||`.
#
# This is the whole reason the tree bottoms out at degree 1: it recurses until
# the diagonal entries are scalars, and only then takes a square root.
#
# ---------------------------------------------------------------------------
# This is the consumer of everything module 7 warned about
# ---------------------------------------------------------------------------
#
# `ffSampling` walks the tree, and at each of the `n` leaves it calls `samplerz`
# with a centre `mu` that is a *Float64 computed from the secret basis*.  There
# is no correcting loop (docs/math/06_ntrugen.md 6.5, 07_samplerz.md 7.4): the
# sampled integer goes straight into the signature.
#
# So this is where the two halves of the FN-DSA story meet:
#
#   * the centre `mu` must be computed to bit-identical Float64 on every
#     platform, or signing is non-deterministic, or the key leaks (module 5);
#   * and the sampling at that centre must not leak `mu` through timing,
#     because `mu` *is* the secret basis, one coordinate at a time (module 7).
#
# CONSTANT TIME: out of scope, and here the difficulty is structural rather
# than local.  The recursion's shape is public (it depends only on `n`), which
# helps; but every arithmetic operation on the way down handles secret data in
# floating point, and the reference implementation's answer is to emulate
# floating point in integers -- for the whole of this module, not just for
# `samplerz`.

# ---------------------------------------------------------------------------
# The tree
# ---------------------------------------------------------------------------

"""
    FalconTree

A node of the ffLDL decomposition tree.  Either an [`FFLDLNode`](@ref) (an
internal node, carrying the `L10` entry of a 2x2 LDL* factorisation and two
subtrees) or an [`FFLDLLeaf`](@ref) (a degree-1 diagonal entry).

The reference represents this as nested Python lists, distinguishing internal
nodes from leaves by `len(tree) == 3`.  We use types instead: the shape is
identical, and a leaf can no longer be indexed as if it were a node.
"""
abstract type FalconTree end

"""
    FFLDLNode(l10, left, right)

Internal node: `l10` is the `L[1][0]` entry of the LDL* factorisation at this
level, `left` and `right` are the subtrees for `D00` and `D11` after splitting
them down the tower.

[Py-ref] scripts/pyref/ffsampling.py:113-137 (`ffldl_fft`), the 3-element list
"""
struct FFLDLNode <: FalconTree
    l10::Vector{ComplexF64}
    left::FalconTree
    right::FalconTree
end

"""
    FFLDLLeaf

A leaf of the tree: a degree-1 diagonal entry of the decomposition.

Before normalisation `value` holds the (length-2, effectively real) diagonal
entry `D`; after [`normalize_tree!`](@ref), `sigma` holds
`sigma_signature / sqrt(real(D[1]))`, which is the width `ffSampling` hands to
`samplerz` at this leaf.

[Py-ref] scripts/pyref/falcon.py:175-190 (`normalize_tree`)
"""
mutable struct FFLDLLeaf <: FalconTree
    value::Vector{ComplexF64}
    sigma::Float64
    "`sqrt(D) * (1/sigma)`, the reciprocal of `sigma`, stored as the C reference stores it"
    isigma::Float64
end

FFLDLLeaf(value::Vector{ComplexF64}) = FFLDLLeaf(value, NaN, NaN)

"Number of leaves, i.e. the degree the tree was built for."
nleaves(t::FFLDLLeaf) = 1
nleaves(t::FFLDLNode) = nleaves(t.left) + nleaves(t.right)

"Depth of the tree."
treedepth(t::FFLDLLeaf) = 0
treedepth(t::FFLDLNode) = 1 + max(treedepth(t.left), treedepth(t.right))

"""
    leaf_sigmas(t) -> Vector{Float64}

The leaves' `sigma` values in left-to-right order.  Used by the tests: it is
the compact fingerprint of a normalised tree.
"""
leaf_sigmas(t::FFLDLLeaf) = [t.sigma]
leaf_sigmas(t::FFLDLNode) = vcat(leaf_sigmas(t.left), leaf_sigmas(t.right))

"""
    node_l10s(t) -> Vector{Vector{ComplexF64}}

The internal nodes' `l10` vectors in pre-order.  The other half of the
fingerprint.
"""
node_l10s(t::FFLDLLeaf) = Vector{Vector{ComplexF64}}()
node_l10s(t::FFLDLNode) = vcat([t.l10], node_l10s(t.left), node_l10s(t.right))

# ---------------------------------------------------------------------------
# Gram matrix and LDL*
# ---------------------------------------------------------------------------

"""
    gram_fft(B) -> Matrix{Vector{ComplexF64}}

The Gram matrix `G[i,j] = sum_k B[i,k] * adj(B[j,k])`, in the FFT domain.

`B` is a 2x2 matrix whose entries are FFT-domain polynomials.  For FALCON,
`B = [g -f; G -F]`.

[Py-ref] scripts/pyref/ffsampling.py:15-32 (`gram`, coefficient domain)

`G` is Hermitian and positive definite by construction: in the FFT domain
`adj` is coordinatewise conjugation (fft.jl), so `G[i,i] = sum_k |B[i,k]|^2` is
real and positive, and `G[j,i] = conj(G[i,j])`.  That is what makes the LDL*
below well defined, and it is the payoff of getting `polyadj`'s sign right back
in module 3.
"""
function gram_fft(B::Matrix{Vector{ComplexF64}})
    size(B) == (2, 2) || throw(ArgumentError("gram_fft expects a 2x2 matrix"))
    n = length(B[1, 1])
    G = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    for i in 1:2, j in 1:2
        acc = zeros(ComplexF64, n)
        for k in 1:2
            acc = add_fft(acc, mul_fft(B[i, k], adj_fft(B[j, k])))
        end
        G[i, j] = acc
    end
    return G
end

"""
    ldl_fft(G) -> (L10, D00, D11)

LDL* factorisation of a 2x2 Hermitian matrix over the ring, in the FFT domain:

    G = L * D * adj(L)^T,   L = [1 0; L10 1],   D = diag(D00, D11)

with

    L10 = G10 / G00
    D00 = G00
    D11 = G11 - L10 * adj(L10) * G00

Corresponds to algorithm `LDL*` of the specification.
[Py-ref] scripts/pyref/ffsampling.py:59-84 (`ldl_fft`)

No square roots appear -- that is the point of LDL* over Cholesky, since the
ring has no square roots.  They are deferred to the leaves, where the entries
are scalars (`normalize_tree!`).
"""
function ldl_fft(G::Matrix{Vector{ComplexF64}})
    size(G) == (2, 2) || throw(ArgumentError("ldl_fft expects a 2x2 matrix"))
    D00 = copy(G[1, 1])

    # THIS BRANCH computes D11 the way the C reference does, which is not the
    # way the specification writes it.  Algorithm 8 (and the Python reference,
    # ffsampling.py:80-82) says
    #
    #     L10 = G10 / G00
    #     D11 = G11 - L10 * adj(L10) * G00           -- three multiplications
    #
    # while C keeps the quotient the other way round and folds:
    #
    #     mu  = G01 / G00
    #     D11 = G11 - mu * adj(G01)                  -- one multiplication
    #     L10 = adj(mu)
    #
    # [C-ref] scripts/cref/fft.c, `Zf(poly_LDL_fft)`
    #
    # They are equal in exact arithmetic -- `mu*adj(G01) = |G01|^2/G00` and
    # `L10*adj(L10)*G00 = |G01|^2/|G00|^2 * G00`, and G00 is real because it is
    # Hermitian -- and they differ in the last bit in floating point.  Since the
    # sampler *rounds* what comes out of here, a last-bit difference becomes a
    # different signature (docs/debug_log.md #048).
    #
    # `G01` is recovered as `adj(G[2,1])` rather than read from `G[1,2]`: the
    # conjugate is a sign flip and therefore exact, and this way the routine
    # does not depend on the caller having filled the redundant upper corner.
    if !LDL_CREF[]
        # The specification's own spelling, for the ablation.  Algorithm 8 and
        # ffsampling.py:80-82.  Three multiplications where C does one.
        L10 = div_fft(G[2, 1], G[1, 1])
        D11 = sub_fft(G[2, 2], mul_fft(mul_fft(L10, adj_fft(L10)), G[1, 1]))
        return (L10, D00, D11)
    end
    G01 = adj_fft(G[2, 1])
    mu = div_fft(G01, G[1, 1])
    L10 = adj_fft(mu)
    D11 = sub_fft(G[2, 2], mul_fft(mu, G[2, 1]))
    return (L10, D00, D11)
end

"""
    LDL_CREF

Whether [`ldl_fft`](@ref) computes `D11` the way the C reference does (`true`,
the default) or the way Algorithm 8 of the specification writes it (`false`).

As with [`CDIV_CREF`](@ref), this exists for the ablation in
scripts/divergence_rate.jl, not for production use.
"""
const LDL_CREF = Ref(true)

"""
    ffldl_fft(G) -> FalconTree

Build the ffLDL decomposition tree of the Gram matrix `G`.

Corresponds to algorithm `ffLDL` of the specification.
[Py-ref] scripts/pyref/ffsampling.py:113-137 (`ffldl_fft`)

## The recursion

Factor `G` as `L*D*adj(L)`, keep `L10`, then descend into the two diagonal
entries.  Descending means: split `D00` down the tower with `split_fft`, giving
`(d00, d01)`, and rebuild a 2x2 Gram matrix one level down:

    G0 = [d00  d01;  adj(d01)  d00]

That matrix is the Gram matrix of the sublattice generated by `D00` viewed over
the smaller field -- the multiplication-by-`D00` operator written in the basis
`{1, x}` of the quadratic extension.  Its Hermitian shape (equal diagonal,
conjugate off-diagonal) is exactly the shape of multiplication by a ring
element, which is why the recursion closes.

The recursion stops at degree 2, where the two diagonal entries are (real)
scalars and become leaves.
"""
function ffldl_fft(G::Matrix{Vector{ComplexF64}})
    n = length(G[1, 1])
    L10, D00, D11 = ldl_fft(G)
    if n > 2
        d00, d01 = split_fft(D00)
        d10, d11 = split_fft(D11)
        G0 = Matrix{Vector{ComplexF64}}(undef, 2, 2)
        G0[1, 1] = d00; G0[1, 2] = d01
        G0[2, 1] = adj_fft(d01); G0[2, 2] = d00
        G1 = Matrix{Vector{ComplexF64}}(undef, 2, 2)
        G1[1, 1] = d10; G1[1, 2] = d11
        G1[2, 1] = adj_fft(d11); G1[2, 2] = d10
        return FFLDLNode(L10, ffldl_fft(G0), ffldl_fft(G1))
    elseif n == 2
        return FFLDLNode(L10, FFLDLLeaf(D00), FFLDLLeaf(D11))
    else
        throw(ArgumentError("ffldl_fft needs degree at least 2, got $n"))
    end
end

"""
    normalize_tree!(tree, sigma) -> tree

Turn the leaves from `||b_i||^2` into the sampling widths `sigma / ||b_i||`.

[Py-ref] scripts/pyref/falcon.py:175-190 (`normalize_tree`)

This is where the square roots that LDL* avoided finally happen, and they are
harmless because by now the entries are positive reals rather than ring
elements.  The value stored is what `ffSampling` passes to `samplerz` as its
`sigma`; the constraint `sigma_min < that < sigma_max` from params.jl is what
the key-generation Gram-Schmidt bound `1.17^2 q` was chosen to guarantee.
"""
function normalize_tree!(tree::FalconTree, sigma::Real)
    # `logn` here is the degree of the *whole* tree, not of the current node --
    # the C reference calls it `orig_logn` and threads it down for exactly the
    # same reason: `fpr_inv_sigma` is indexed by it.
    logn = round(Int, log2(nleaves(tree)))
    return _normalize_tree!(tree, sigma, logn)
end

function _normalize_tree!(tree::FFLDLNode, sigma::Real, logn::Int)
    _normalize_tree!(tree.left, sigma, logn)
    _normalize_tree!(tree.right, sigma, logn)
    return tree
end

function _normalize_tree!(leaf::FFLDLLeaf, sigma::Real, logn::Int)
    x = real(leaf.value[1])
    # The specification's value: the width handed to `SamplerZ`.
    leaf.sigma = sigma / sqrt(x)
    # THIS BRANCH also stores the C reference's value, which is its reciprocal
    # and is *not* computed as one:
    #
    #     [C-ref] scripts/cref/sign.c, `ffLDL_binary_normalize`:
    #         tree[0] = fpr_mul(fpr_sqrt(tree[0]), fpr_inv_sigma[orig_logn]);
    #
    # "We actually store in the tree leaf the inverse of the value mandated by
    # the specification: this saves a division both here and in the sampler."
    # `fpr_inv_sigma[logn]` is a *table constant*, and -- this is the part that
    # is easy to get wrong -- it is **not** the correctly rounded reciprocal of
    # the sigma in Table 3.3.  At logn = 9 the table entry and `1/sigma`
    # differ by one ulp; at logn = 10, by two.  Presumably the table was
    # computed from a sigma with more digits than Table 3.3 prints.  Whatever
    # the reason, a reciprocal computed here does not reproduce it, so the
    # literals are transcribed (`INV_SIGMA_CREF`).
    #
    # And C multiplies by it where the specification divides, so
    # `sqrt(x) * inv_sigma` is not bit-equal to `1 / (sigma / sqrt(x))`
    # either (docs/debug_log.md #050).
    leaf.isigma = sqrt(x) * _inv_sigma_cref(sigma, logn)
    leaf.value[2] = 0
    return leaf
end

"""
    falcon_tree(B_fft, sigma) -> FalconTree

Build and normalise the Falcon tree from a secret basis in the FFT domain.
Convenience wrapper: `normalize_tree!(ffldl_fft(gram_fft(B)), sigma)`.
"""
falcon_tree(B::Matrix{Vector{ComplexF64}}, sigma::Real) =
    normalize_tree!(ffldl_fft(gram_fft(B)), sigma)

# ---------------------------------------------------------------------------
# Nearest plane (deterministic) and sampling (randomised)
# ---------------------------------------------------------------------------

"""
    ffnp_fft(t, T) -> (z0, z1)

Fast Fourier nearest plane: the *deterministic* walk down the tree, rounding to
the nearest integer at each leaf instead of sampling.

[Py-ref] scripts/pyref/ffsampling.py:163-186 (`ffnp_fft`)

Not used in signing -- rounding deterministically is exactly what leaks the
basis, which is why GPV requires the randomised version.  It is kept because it
is the same recursion with the randomness removed, so it can be tested without
a byte stream and its output can be checked against the plain (non-fast)
nearest plane on small degrees.  When `ffSampling` misbehaves, this is the
first thing to compare against.
"""
function ffnp_fft(t::NTuple{2,Vector{ComplexF64}}, T::FalconTree)
    n = length(t[1])
    if n > 1
        node = T::FFLDLNode
        t1a, t1b = split_fft(t[2])
        z1 = merge_fft(ffnp_fft((t1a, t1b), node.right)...)
        t0b = add_fft(t[1], mul_fft(sub_fft(t[2], z1), node.l10))
        t0a, t0bb = split_fft(t0b)
        z0 = merge_fft(ffnp_fft((t0a, t0bb), node.left)...)
        return (z0, z1)
    else
        return (ComplexF64[round(real(t[1][1]))], ComplexF64[round(real(t[2][1]))])
    end
end


# ---------------------------------------------------------------------------
# The bottom of the recursion, as the C reference writes it
# ---------------------------------------------------------------------------

"""
    FFSAMPLING_CREF

Whether `ffsampling_fft` uses the C reference's spelling of the bottom two
levels (`_ffsampling_c4`, with the leaves' reciprocal widths) or the
specification's generic recursion.

**Defaults to `true` on this branch**, which is what makes signing reproduce
the C reference's bytes.  Set it to `false` -- see [`with_spec_ffsampling`](@ref)
-- to get the specification's route back; the recorded Python-reference vectors
are replayed that way, since they were produced by it.
"""
const FFSAMPLING_CREF = Ref(true)

"""
    with_spec_ffsampling(f)

Run `f()` with `ffsampling_fft` on the specification's generic recursion rather
than the C reference's spelling, restoring the previous setting afterwards.
"""
function with_spec_ffsampling(f)
    old = FFSAMPLING_CREF[]
    FFSAMPLING_CREF[] = false
    try
        return f()
    finally
        FFSAMPLING_CREF[] = old
    end
end

"""
    with_spec_spelling(f; cdiv = false, ldl = false, ffsampling = false)

Run `f()` with any subset of the three algebraically-equal-but-differently-
rounded formulations switched from the C reference's spelling to the
specification's, restoring all three afterwards.

The three are independent, and the point of separating them is that they are
not equally novel:

  * `ffsampling` -- the hand-unrolled bottom levels.  This is the SAME
    difference ePrint 2024/1709 §6.1 identifies between the reference's
    `sign_dyn` and `sign_tree` modes, down to the folded `1/(2*sqrt(2))`
    constant, and §7.2 of that paper gives a countermeasure for it.  Measuring
    it reproduces a published result; it does not establish a new one.
  * `cdiv`, `ldl` -- complex division and the `D11` entry of LDL*.  These are
    differences between the specification (with the Python reference) and the C
    reference, present in BOTH of C's signing modes, and they do not appear in
    2024/1709.  They perturb the tree's `l10` and its leaf widths, so they
    reach the sampler's centres by a different route than `ffsampling` does.

So `with_spec_spelling(f; cdiv = true, ldl = true)` -- holding `ffsampling` at
the C spelling -- is the ablation that decides whether anything here is not
already 2024/1709 §6.1.  See scripts/divergence_rate.jl.
"""
function with_spec_spelling(f; cdiv::Bool = false, ldl::Bool = false,
                            ffsampling::Bool = false)
    old = (CDIV_CREF[], LDL_CREF[], FFSAMPLING_CREF[])
    cdiv && (CDIV_CREF[] = false)
    ldl && (LDL_CREF[] = false)
    ffsampling && (FFSAMPLING_CREF[] = false)
    try
        return f()
    finally
        CDIV_CREF[], LDL_CREF[], FFSAMPLING_CREF[] = old
    end
end

"""
    INV_SIGMA_CREF

`1/sigma` per degree, indexed by `logn + 1`, transcribed from the C reference's
`fpr_inv_sigma[]` table.

[C-ref] scripts/cref/fpr.h, `fpr_inv_sigma[]` (the FALCON_FPNATIVE branch)

**These are not the reciprocals of Table 3.3's sigmas.**  Computing
`1/165.7366171829776` gives a Float64 one ulp below the entry below, and at
`logn = 10` two ulp below.  Table 3.3 prints sigma to twelve significant
figures; the table here was evidently derived from something more precise.
Since the reference multiplies the tree leaves by *this* number, reproducing
its signatures means carrying *this* number (docs/debug_log.md #050).

Index 1 (`logn = 0`) is the reference's unused zero slot, kept so the indexing
matches the C source line for line.
"""
const INV_SIGMA_CREF = (
    0.0,                                     # logn = 0, unused in the reference
    0.0069054793295940891952143765991630516, # logn = 1
    0.0068102267767177975961393730687908629,
    0.0067188101910722710707826117910434131,
    0.0065883354370073665545865037227681924,
    0.0064651781207602900738053897763485516,
    0.0063486788828078995327741182928037856,
    0.0062382586529084374473367528433697537,
    0.0061334065020930261548984001431770281,
    0.0060336696681577241031668062510953022, # logn = 9  (FALCON-512)
    0.0059386453095331159950250124336477482, # logn = 10 (FALCON-1024)
)

"""
The `fpr_inv_sigma` entry the C reference would use for this width.

Indexing by the *tree's* degree is not enough: the toy trees the recorded
Python vectors use are built at degree 8 with FALCON-512's `sigma`, and the C
table is indexed by the degree whose `sigma` it is.  So the lookup is by the
value of `sigma`, and anything that is not a standardised width falls back to
the reciprocal -- which is what the specification's route would do anyway.
"""
function _inv_sigma_cref(sigma::Real, logn::Int)
    sig = Float64(sigma)
    sig == FALCON_512.sigma  && return INV_SIGMA_CREF[10]   # logn = 9
    sig == FALCON_1024.sigma && return INV_SIGMA_CREF[11]   # logn = 10
    1 <= logn <= 10 && return INV_SIGMA_CREF[logn + 1]
    return 1 / sig
end

"1/sqrt(2), as the C reference spells it.  [C-ref] scripts/cref/fpr.h, `fpr_invsqrt2`"
const INVSQRT2 = 0.707106781186547524400844362105

"1/sqrt(8).  [C-ref] scripts/cref/fpr.h, `fpr_invsqrt8`"
const INVSQRT8 = 0.353553390593273762200422181052

"""
    _ffsampling_c4(t0, t1, node, sigmin, randombytes) -> (z0, z1)

The last two levels of `ffsampling_fft`, at degree 4, transcribed operation by
operation from the C reference's hand-unrolled block.

[C-ref] scripts/cref/sign.c, `ffSampling_fft`, the `logn == 2` case.

## NOT WIRED IN.  Read this before changing that.

`ffsampling_fft` does **not** call this.  It is reachable, tested against the C
reference, and deliberately left disconnected, because what has been verified
is not what would be needed to connect it.

What is verified: with a *deterministic* stand-in for the sampler (plain
`floor`, which consumes no randomness), substituting this block makes the whole
recursion agree with C bit for bit at every degree from 4 to 512, where the
generic recursion disagrees from degree 4 upward.  So the **arithmetic** is
right.

What is not verified: the **order and count of the `samplerz` calls**.  Each
call consumes bytes, so a transcription that computes the same numbers while
asking for randomness in a different order still produces a different
signature -- and the `floor` measurement is blind to exactly that.  Wiring
this in made `FFSAMPLING_KAT` (recorded from the Python reference) fail, which
is expected on this branch and therefore says nothing either way about whether
the four calls below are in C's order.

Settling it needs C driven with the real sampler off a controlled PRNG, which
is the next step and has not been done.  Until then this stays disconnected:
shipping it would be trading a divergence that is understood for one that is
not (docs/debug_log.md #048).

## Why it exists at all

THIS BRANCH is trying to reproduce the C reference's signature *bytes*.  Every
floating-point primitive underneath was checked against C and found bit-exact
-- `fft`, `ifft`, `split_fft`, `merge_fft`, `add_fft`, `sub_fft`, `mul_fft`,
`adj_fft`, `div_fft` (once spelled C's way, see fft.jl) and `ldl_fft` (ditto).
The recursion built from them still disagreed, and bisecting by degree put the
first disagreement at degree 4 -- exactly where C stops calling its own
`poly_split_fft` / `poly_merge_fft` and runs a flat block of scalar
operations instead (docs/debug_log.md #048).

The block is not doing different mathematics.  It is doing the same split and
merge with the constants pre-folded:

  * `split`: the generic routine computes `0.5*(a - b)` and then multiplies by
    `conj(w)` where `w = exp(i*pi/4)`.  C folds the halving and the twiddle
    into one constant, `1/sqrt(8)`, and does one multiplication where the
    generic path does two.
  * `merge`: the generic routine multiplies by `w`; C computes
    `(b_re - b_im)/sqrt(2)` and `(b_re + b_im)/sqrt(2)`, two real
    multiplications instead of a complex one.

Fewer roundings, different roundings.  Substituting this block made the whole
recursion agree with C at every degree from 4 to 512.

CONSTANT TIME: unchanged from the generic path -- the data-dependent work is
all inside `samplerz`.

NOTE ON THE ORDER OF SAMPLER CALLS: the four `samplerz` calls below are in C's
order, and that order matters because each call consumes bytes from
`randombytes`.  The bit-exactness measurement that motivated this block used a
deterministic stand-in for the sampler, so it verified the *arithmetic* and
not the *consumption order*; the order here is transcribed from the C source,
not measured.
"""
function _ffsampling_c4(t0::Vector{ComplexF64}, t1::Vector{ComplexF64},
                        node::FFLDLNode, sigmin::Real, randombytes)
    L = node.l10
    tree0 = node.left::FFLDLNode
    tree1 = node.right::FFLDLNode
    l0 = tree0.l10[1]
    l1 = tree1.l10[1]
    # the C reference's leaves: 1/sigma, not sigma
    s0a = (tree0.left::FFLDLLeaf).isigma;  s0b = (tree0.right::FFLDLLeaf).isigma
    s1a = (tree1.left::FFLDLLeaf).isigma;  s1b = (tree1.right::FFLDLLeaf).isigma
    for s in (s0a, s0b, s1a, s1b)
        isnan(s) && throw(ArgumentError(
            "ffsampling_fft on an un-normalised tree: call normalize_tree! first"))
    end

    # --- split t1, sample, merge into z1
    a_re = real(t1[1]); a_im = imag(t1[1]); b_re = real(t1[2]); b_im = imag(t1[2])
    c_re = a_re + b_re; c_im = a_im + b_im
    w0 = 0.5 * c_re;    w1 = 0.5 * c_im
    c_re = a_re - b_re; c_im = a_im - b_im
    w2 = (c_re + c_im) * INVSQRT8
    w3 = (c_im - c_re) * INVSQRT8
    x0 = w2; x1 = w3
    w2 = Float64(samplerz_isigma(x0, s1b, sigmin, randombytes))
    w3 = Float64(samplerz_isigma(x1, s1b, sigmin, randombytes))
    a_re = x0 - w2; a_im = x1 - w3
    b_re = real(l1); b_im = imag(l1)
    c_re = a_re * b_re - a_im * b_im
    c_im = a_re * b_im + a_im * b_re
    x0 = c_re + w0; x1 = c_im + w1
    w0 = Float64(samplerz_isigma(x0, s1a, sigmin, randombytes))
    w1 = Float64(samplerz_isigma(x1, s1a, sigmin, randombytes))
    a_re = w0; a_im = w1; b_re = w2; b_im = w3
    c_re = (b_re - b_im) * INVSQRT2
    c_im = (b_re + b_im) * INVSQRT2
    z1u = complex(a_re + c_re, a_im + c_im)
    z1v = complex(a_re - c_re, a_im - c_im)
    z1 = ComplexF64[z1u, z1v, conj(z1u), conj(z1v)]

    # --- tb0 = t0 + (t1 - z1) * L
    w0 = real(t1[1]) - real(z1[1]); w1 = real(t1[2]) - real(z1[2])
    w2 = imag(t1[1]) - imag(z1[1]); w3 = imag(t1[2]) - imag(z1[2])
    a_re = w0; a_im = w2; b_re = real(L[1]); b_im = imag(L[1])
    w0 = a_re * b_re - a_im * b_im
    w2 = a_re * b_im + a_im * b_re
    a_re = w1; a_im = w3; b_re = real(L[2]); b_im = imag(L[2])
    w1 = a_re * b_re - a_im * b_im
    w3 = a_re * b_im + a_im * b_re
    w0 += real(t0[1]); w1 += real(t0[2]); w2 += imag(t0[1]); w3 += imag(t0[2])

    # --- split tb0, sample, merge into z0
    a_re = w0; a_im = w2; b_re = w1; b_im = w3
    c_re = a_re + b_re; c_im = a_im + b_im
    w0 = 0.5 * c_re;    w1 = 0.5 * c_im
    c_re = a_re - b_re; c_im = a_im - b_im
    w2 = (c_re + c_im) * INVSQRT8
    w3 = (c_im - c_re) * INVSQRT8
    x0 = w2; x1 = w3
    w2 = Float64(samplerz_isigma(x0, s0b, sigmin, randombytes))
    w3 = Float64(samplerz_isigma(x1, s0b, sigmin, randombytes))
    a_re = x0 - w2; a_im = x1 - w3
    b_re = real(l0); b_im = imag(l0)
    c_re = a_re * b_re - a_im * b_im
    c_im = a_re * b_im + a_im * b_re
    x0 = c_re + w0; x1 = c_im + w1
    w0 = Float64(samplerz_isigma(x0, s0a, sigmin, randombytes))
    w1 = Float64(samplerz_isigma(x1, s0a, sigmin, randombytes))
    a_re = w0; a_im = w1; b_re = w2; b_im = w3
    c_re = (b_re - b_im) * INVSQRT2
    c_im = (b_re + b_im) * INVSQRT2
    z0u = complex(a_re + c_re, a_im + c_im)
    z0v = complex(a_re - c_re, a_im - c_im)
    z0 = ComplexF64[z0u, z0v, conj(z0u), conj(z0v)]

    return (z0, z1)
end

"""
    ffsampling_fft(t, T, sigmin, randombytes) -> (z0, z1)

Fast Fourier sampling: walk the tree, sampling one integer per leaf from
`D_{Z, mu, sigma_leaf}` where `mu` is the current residual and `sigma_leaf` is
the normalised leaf value.

Corresponds to algorithm `ffSampling` of the specification.
[Py-ref] scripts/pyref/ffsampling.py:189-214 (`ffsampling_fft`)

## The order matters

The recursion does `z1` **first**, then computes the corrected centre

    t0b = t0 + (t1 - z1) * L10

and only then does `z0`.  That is Klein's algorithm: the second coordinate is
sampled against its own Gram-Schmidt width, and the first is then re-centred by
the part of `z1`'s error that lies along it.  Doing them in the other order, or
forgetting the correction, gives a sampler that still returns short-ish vectors
and still verifies -- and leaks the basis.  There is no test that catches that
except agreement with the reference.

CONSTANT TIME: see the module header.  `mu` here is secret, and it is a
Float64.
"""
function ffsampling_fft(t::NTuple{2,Vector{ComplexF64}}, T::FalconTree,
                        sigmin::Real, randombytes)
    n = length(t[1])
    if n == 4 && FFSAMPLING_CREF[]
        return _ffsampling_c4(t[1], t[2], T::FFLDLNode, sigmin, randombytes)
    elseif n > 1
        node = T::FFLDLNode
        t1a, t1b = split_fft(t[2])
        z1 = merge_fft(ffsampling_fft((t1a, t1b), node.right, sigmin, randombytes)...)
        t0b = add_fft(t[1], mul_fft(sub_fft(t[2], z1), node.l10))
        t0a, t0bb = split_fft(t0b)
        z0 = merge_fft(ffsampling_fft((t0a, t0bb), node.left, sigmin, randombytes)...)
        return (z0, z1)
    else
        leaf = T::FFLDLLeaf
        isnan(leaf.sigma) && throw(ArgumentError(
            "ffsampling_fft on an un-normalised tree: call normalize_tree! first"))
        z0 = samplerz(real(t[1][1]), leaf.sigma, sigmin, randombytes)
        z1 = samplerz(real(t[2][1]), leaf.sigma, sigmin, randombytes)
        return (ComplexF64[z0], ComplexF64[z1])
    end
end
