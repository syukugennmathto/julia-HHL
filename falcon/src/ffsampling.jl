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
end

FFLDLLeaf(value::Vector{ComplexF64}) = FFLDLLeaf(value, NaN)

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
    L10 = div_fft(G[2, 1], G[1, 1])
    D11 = sub_fft(G[2, 2], mul_fft(mul_fft(L10, adj_fft(L10)), G[1, 1]))
    return (L10, D00, D11)
end

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
function normalize_tree!(tree::FFLDLNode, sigma::Real)
    normalize_tree!(tree.left, sigma)
    normalize_tree!(tree.right, sigma)
    return tree
end

function normalize_tree!(leaf::FFLDLLeaf, sigma::Real)
    leaf.sigma = sigma / sqrt(real(leaf.value[1]))
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
    if n > 1
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
