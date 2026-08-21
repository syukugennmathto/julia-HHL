# falcon.jl -- keygen, sign, verify.
#
# ---------------------------------------------------------------------------
# The scheme, in one page
# ---------------------------------------------------------------------------
#
# FALCON is hash-and-sign over an NTRU lattice.
#
#   * The public key is a single polynomial `h = g/f mod q`.  The lattice it
#     describes is `{ (s1, s2) : s1 + s2*h = 0 mod q }`.
#   * A signature of a message is a *short* pair `(s1, s2)` with
#
#         s1 + s2*h = c  (mod q),      ||(s1, s2)||^2 <= beta^2
#
#     where `c = HashToPoint(salt || message)`.  Only `s2` is transmitted:
#     `s1` is recomputed by the verifier as `c - s2*h`.
#   * Anyone can find *a* solution to the congruence (set `s2 = 0`, `s1 = c`);
#     the hard part is finding a *short* one, and that needs the secret short
#     basis.
#
# Verification is therefore trivial arithmetic -- one NTT product, one
# subtraction, one norm -- and signing is everything in modules 6 to 8.
#
# ---------------------------------------------------------------------------
# What this file can and cannot reproduce
# ---------------------------------------------------------------------------
#
# Verification is integer-only, so it must agree with the C reference exactly,
# and `test/vectors/cref_kat.jl` holds a genuine C-generated FALCON-512
# signature for us to accept.  That is real interoperability.
#
# Signing is *not* reproducible.  It runs through the FFT, and module 8
# measured that a 1-ulp difference in the FFT root table changes the output
# (docs/debug_log.md #025).  We therefore do not attempt to match anyone's
# signature bytes; we check that our signatures verify, that they are short,
# and that the reference's verify would accept them (by construction, since
# verification is shared).
#
# The C reference goes further than we do here: `falcon.c` wraps its signing
# call in `set_fpu_cw(2)`, forcing the x87 control word to double precision, on
# top of emulating floating point in integer arithmetic.  Two independent
# belt-and-braces measures against exactly the divergence #025 exhibits.

# ---------------------------------------------------------------------------
# Hashing a message to a lattice point
# ---------------------------------------------------------------------------

"""
    hash_to_point(message, salt, n; q = Q) -> Vector{Int}

Map `salt || message` to a point of `Z_q^n`, by reading SHAKE256 two bytes at a
time and rejecting values that would bias the result.

[Py-ref] scripts/pyref/falcon.py:252-278 (`__hash_to_point__`)
[C-ref] inner.h:536 (`Zf(hash_to_point_vartime)`)

## The rejection bound

With `k = floor(2^16 / q) = 5`, a 16-bit value is accepted only if it is below
`k*q = 61445`, and then reduced mod q.  Accepting all 65536 values and reducing
would make the residues `0 .. 3651` very slightly more likely than the rest --
a bias of about `2^-13` per coefficient, which is enormous next to the `2^-45`
statistical distance the parameters were chosen for.

CONSTANT TIME: this is the *vartime* variant, and that is fine here: the salt
and message are public, so the number of rejections leaks nothing.  The C
reference also carries `hash_to_point_ct` for the case where the message is
secret.
"""
function hash_to_point(message::AbstractVector{UInt8}, salt::AbstractVector{UInt8},
                       n::Integer; q::Integer = Q)
    q < (1 << 16) || throw(ArgumentError("hash_to_point needs q < 2^16"))
    k = (1 << 16) ÷ q
    bound = k * q
    xof = shake256_xof(salt, message)
    out = Vector{Int}(undef, n)

    # The XOF is squeezed in blocks rather than two bytes at a time -- the same
    # byte stream read differently, since a sponge does not care where the
    # caller puts the boundaries.  Two bytes at a time meant one allocation and
    # one call per coefficient (docs/debug_log.md #038).
    #
    # The block sizes are chosen in whole Keccak permutations, because that is
    # the only granularity the sponge actually has.  A draw is accepted with
    # probability 61445/65536 ~ 0.9375, so the expected need is 2n/0.9375 =
    # 2.133n bytes -- 1092 at n = 512, or 8.03 rate blocks.  The first fill is
    # therefore 8 blocks and the top-ups are one block each, which costs on
    # average about 8.5 permutations per call.  Filling a single generous
    # buffer instead cost 10, and sizing it at exactly the expected need cost
    # 16 half the time (#039).
    rate = SHAKE256_RATE
    first = max(rate, (2 * Int(n) * 17) >> 4)
    first = ((first + rate - 1) ÷ rate) * rate            # whole permutations
    buf = Vector{UInt8}(undef, first)
    squeeze!(xof, buf, 0, first)
    have = first
    i = 1
    pos = 1
    @inbounds while i <= n
        if pos + 1 > have
            # Everything drawn so far is consumed, so the buffer is reused from
            # the start.  `first` is even and pairs are consumed two at a time,
            # so this can never split a pair across the refill.
            squeeze!(xof, buf, 0, rate)
            have = rate
            pos = 1
        end
        elt = (Int(buf[pos]) << 8) + Int(buf[pos + 1])   # big-endian, per the reference
        pos += 2
        if elt < bound
            # `elt % q` would be an integer division -- `q` arrives as a keyword
            # argument, so it is not a compile-time constant and LLVM cannot
            # turn it into a multiply-shift.  But `elt < bound = k*q` with
            # k = 5, so the remainder is at most four conditional subtractions
            # away, and this says the same thing the rejection bound says.
            r = elt
            while r >= q
                r -= q
            end
            out[i] = r
            i += 1
        end
    end
    return out
end

hash_to_point(message::AbstractString, salt, n; q::Integer = Q) =
    hash_to_point(codeunits(message), salt, n; q = q)

# ---------------------------------------------------------------------------
# Keys
# ---------------------------------------------------------------------------

"""
    FalconPublicKey

Just `h`, plus the parameter set it belongs to.
"""
struct FalconPublicKey
    params::FalconParams
    h::Vector{Int}
end

"""
    FalconPrivateKey

The short basis `(f, g, F, G)`, together with the two derived objects signing
needs: the basis in the FFT domain and the normalised ffLDL tree.

Both derivations are deterministic functions of `(f, g, F, G)`, so they are
*caches*, not key material -- but they are expensive (the tree is the
`O(n log n)` part of module 8), which is why a real implementation expands the
key once and signs many times.
"""
struct FalconPrivateKey
    params::FalconParams
    f::Vector{Int}
    g::Vector{Int}
    F::Vector{Int}
    G::Vector{Int}
    B0_fft::Matrix{Vector{ComplexF64}}
    tree::FalconTree
end

"""
    expand_privkey(f, g, F, G, params) -> FalconPrivateKey

Build the FFT basis and the normalised Falcon tree from a short basis.

    B0 = [ g   -f ]
         [ G   -F ]

[Py-ref] scripts/pyref/falcon.py:374-388 (inside `keygen`)
"""
function expand_privkey(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer},
                        F::AbstractVector{<:Integer}, G::AbstractVector{<:Integer},
                        params::FalconParams)
    n = params.n
    (length(f) == length(g) == length(F) == length(G) == n) ||
        throw(ArgumentError("basis polynomials must all have length $n"))
    B = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    B[1, 1] = fft(Float64.(g)); B[1, 2] = fft(Float64.(-f))
    B[2, 1] = fft(Float64.(G)); B[2, 2] = fft(Float64.(-F))
    tree = falcon_tree(B, params.sigma)
    return FalconPrivateKey(params, Int.(f), Int.(g), Int.(F), Int.(G), B, tree)
end

"""
    falcon_keygen(n, randombytes) -> (sk, pk)

Generate a key pair: solve the NTRU equation for a short basis, expand it, and
publish `h = g/f mod q`.

[Py-ref] scripts/pyref/falcon.py:357-390 (`keygen`)
"""
function falcon_keygen(n::Integer, randombytes; sampler::Symbol = :cdt)
    p = params(n)
    f, g, F, G = ntru_gen(n, randombytes; q = p.q, sampler = sampler)
    sk = expand_privkey(Int.(f), Int.(g), Int.(F), Int.(G), p)
    fq = Int[mod(c, p.q) for c in f]
    gq = Int[mod(c, p.q) for c in g]
    h = polydivq(gq, fq, )
    return (sk, FalconPublicKey(p, h))
end

"""
    public_key(sk) -> FalconPublicKey

Recompute the public key from a private key: `h = g/f mod q`.
"""
function public_key(sk::FalconPrivateKey)
    q = sk.params.q
    h = polydivq(Int[mod(c, q) for c in sk.g], Int[mod(c, q) for c in sk.f])
    return FalconPublicKey(sk.params, h)
end

# ---------------------------------------------------------------------------
# Signing
# ---------------------------------------------------------------------------

"""
    sample_preimage(sk, point, randombytes) -> (s1, s2)

Find a short `(s1, s2)` with `s1 + s2*h = point (mod q)`.

[Py-ref] scripts/pyref/falcon.py:316-355 (`__sample_preimage__`)

## What the algebra is doing

We want a lattice point near `(point, 0)`.  The target expressed in the basis'
coordinates is `(point, 0) * B0^{-1}`, and because `det(B0) = q` and `B0` has
its particular NTRU shape, that inverse is available in closed form -- which is
why the code is two multiplications and a division by `q` rather than a matrix
solve:

    t0 = point * d / q,    t1 = -point * b / q      with B0 = [a b; c d]

`ffSampling` then rounds `(t0, t1)` to a nearby *integral* vector `z`, and
`v = z * B0` is the lattice point.  The signature is the residual
`s = (point, 0) - v`, which is short exactly because `z` was sampled close to
`t`.

The congruence `s1 + s2*h = point` holds because `v` is in the lattice, where
that expression vanishes -- so it is inherited, not enforced.  That is worth
noticing: nothing in signing ever computes `h`.
"""
function sample_preimage(sk::FalconPrivateKey, point::AbstractVector{<:Integer},
                         randombytes)
    p = sk.params
    n = p.n
    a, b = sk.B0_fft[1, 1], sk.B0_fft[1, 2]
    c, d = sk.B0_fft[2, 1], sk.B0_fft[2, 2]

    point_fft = fft(Float64.(point))
    qf = Float64(p.q)
    t0 = ComplexF64[(point_fft[i] * d[i]) / qf for i in 1:n]
    t1 = ComplexF64[(-point_fft[i] * b[i]) / qf for i in 1:n]

    z0, z1 = ffsampling_fft((t0, t1), sk.tree, p.sigma_min, randombytes)

    v0 = ifft(add_fft(mul_fft(z0, a), mul_fft(z1, c)))
    v1 = ifft(add_fft(mul_fft(z0, b), mul_fft(z1, d)))
    s1 = Int[Int(point[i]) - round(Int, v0[i]) for i in 1:n]
    s2 = Int[-round(Int, v1[i]) for i in 1:n]
    return (s1, s2)
end

"""
    falcon_sign(sk, message, randombytes; max_attempts = 1000) -> Vector{UInt8}

Sign `message`, returning a PADDED-format signature of exactly
`sk.params.sig_bytes` bytes.

[Py-ref] scripts/pyref/falcon.py:392-420 (`sign`)
[C-ref] falcon.c:433-494 (`falcon_sign_dyn_finish`, the `for (;;)` loop)

## The retry loop -- both of its conditions

Signing restarts when either check fails:

1. **`||(s1, s2)||^2 > beta^2`.**  The sampler is a *Gaussian*, so an
   occasional long vector is expected; it must be discarded, not shipped.
2. **The compressed `s2` does not fit** in `sig_bytes - 41` bytes.  The
   Golomb-Rice encoding is variable-length (module 9) and the format is fixed,
   so this genuinely happens.

Note what is *not* redrawn: **the salt stays the same across retries.**  Only
the sampling randomness is fresh.  Both references agree on this -- the C code
draws the nonce once in `falcon_sign_start` and copies the same 40 bytes into
the buffer inside the loop.  It matters because a fresh salt would mean a fresh
hash point, i.e. a different problem instance; keeping it means we are
re-rolling the *sampler* on one instance.

`max_attempts` guards against an infinite loop on a broken sampler; the
references have no such bound. Measured acceptance is the overwhelming
majority of first attempts.

CONSTANT TIME: the number of attempts depends on the key and the message.  The
reference accepts that leak.
"""
function falcon_sign(sk::FalconPrivateKey, message, randombytes;
                     max_attempts::Integer = 1000)
    p = sk.params
    msg = message isa AbstractString ? collect(codeunits(message)) : collect(message)
    salt = randombytes(SALT_LEN)
    length(salt) == SALT_LEN || throw(ArgumentError("randombytes returned a short salt"))
    point = hash_to_point(msg, salt, p.n; q = p.q)

    for _ in 1:max_attempts
        s1, s2 = sample_preimage(sk, point, randombytes)
        nrm = sqnorm(s1, s2)
        nrm <= p.sig_bound || continue                 # too long: resample
        sig = encode_signature(salt, s2, p.logn, p.sig_bytes)
        sig === nothing && continue                    # does not fit: resample
        return sig
    end
    throw(ErrorException(
        "falcon_sign: no acceptable signature in $max_attempts attempts -- " *
        "this indicates a broken sampler or key, not bad luck"))
end

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

"""
    falcon_verify(pk, message, signature) -> Bool

Verify a signature.  Returns `false` for any malformed input rather than
throwing, because a verifier is fed attacker-controlled bytes.

[Py-ref] scripts/pyref/falcon.py:422-455 (`verify`)
[C-ref] inner.h:592 (`Zf(verify_raw)`)

## The three checks

1. the signature parses, and canonically (module 9);
2. `s1 = c - s2*h mod q`, taken with **centred** representatives -- this is
   where `q-1` has to mean `-1` and not a huge number (poly.jl's `centered`);
3. `||(s1, s2)||^2 <= beta^2`, the bound being *inclusive*
   (params.jl, `[C-ref] common.c:239`).

All integer arithmetic.  No floating point appears anywhere in verification,
which is why this function -- unlike signing -- agrees with the C reference
exactly, and why `test/vectors/cref_kat.jl` can be used as a real
interoperability test.
"""
function falcon_verify(pk::FalconPublicKey, message, signature::AbstractVector{UInt8})
    p = pk.params
    # `collect` rather than passing `message` through: it costs one small copy
    # and buys a concretely typed local.  Returning `codeunits(...)` for a
    # string and the argument itself otherwise leaves `msg` a union, and the
    # dynamic dispatch that follows costs more than the copy
    # (docs/debug_log.md #039).
    msg = message isa AbstractString ? collect(codeunits(message)) : collect(message)

    logn, salt, s2 = try
        decode_signature(signature)
    catch e
        e isa ArgumentError && return false
        rethrow()
    end
    logn == p.logn || return false
    length(signature) == p.sig_bytes || return false

    point = hash_to_point(msg, salt, p.n; q = p.q)

    # `polymulq_fast!` is the iterative in-place NTT of ntt.jl rather than the
    # schoolbook convolution.  It computes the same ring element -- the test
    # suite asserts that at every degree -- and it is the difference between
    # 0.157 ms and 0.021 ms at n = 512, which is most of verification
    # (docs/debug_log.md #038).
    #
    # The subtraction, the centring and the norm are folded into one pass here
    # instead of three array-returning calls.  That is not micro-optimisation
    # for its own sake: at this point the remaining cost of verification is
    # almost entirely the per-array traffic, so each intermediate that is not
    # materialised is a measurable fraction of the whole.
    n = p.n
    prod = Vector{UInt32}(undef, n)
    scratch = Vector{UInt32}(undef, n)
    polymulq_fast!(prod, s2, pk.h, scratch)

    q = p.q
    half = q >> 1
    acc = 0
    @inbounds for i in 1:n
        d = point[i] - Int(prod[i])                 # in (-q, q)
        d < 0 && (d += q)                           # now in [0, q)
        d > half && (d -= q)                        # centred: (-q/2, q/2]
        acc += d * d + Int(s2[i]) * Int(s2[i])
    end
    # `acc` cannot overflow: |d|, |s2[i]| <= q/2 = 6144, so each term is under
    # 2^26 and 2n of them stay under 2^38.  `sig_bound` is a BigInt, and the
    # comparison promotes.
    return acc <= p.sig_bound
end

"""
    signature_norm(pk, message, signature) -> Union{BigInt,Nothing}

The squared norm `||(s1, s2)||^2` that [`falcon_verify`](@ref) compares against
`beta^2`, or `nothing` if the signature does not parse.

Exposed because "how much headroom did that signature have" is the useful
diagnostic when the bound is the thing under suspicion, and because the tests
want to assert the margin rather than just the verdict.
"""
function signature_norm(pk::FalconPublicKey, message, signature::AbstractVector{UInt8})
    p = pk.params
    msg = message isa AbstractString ? collect(codeunits(message)) : collect(message)
    logn, salt, s2 = try
        decode_signature(signature)
    catch e
        e isa ArgumentError && return nothing
        rethrow()
    end
    logn == p.logn || return nothing
    point = hash_to_point(msg, salt, p.n; q = p.q)
    s2q = Int[mod(c, p.q) for c in s2]
    s1 = centered(polysubq(point, polymulq(s2q, pk.h, p.q), p.q), p.q)
    return sqnorm(s1, s2)
end

# ---------------------------------------------------------------------------
# Serialised keys
# ---------------------------------------------------------------------------

"""
    privkey_from_bytes(bytes) -> FalconPrivateKey

Decode a serialised private key and expand it for signing.  `G` is recovered
from the NTRU equation, since the format omits it (module 9).
"""
function privkey_from_bytes(bytes::AbstractVector{UInt8})
    logn, f, g, F = decode_privkey(bytes)
    p = params(1 << logn)
    G = Int.(recover_G(f, g, F; q = p.q))
    return expand_privkey(f, g, F, G, p)
end

"""
    pubkey_from_bytes(bytes) -> FalconPublicKey

Decode a serialised public key.
"""
function pubkey_from_bytes(bytes::AbstractVector{UInt8})
    logn, h = decode_pubkey(bytes)
    return FalconPublicKey(params(1 << logn), h)
end

"Serialise a public key."
pubkey_bytes(pk::FalconPublicKey) = encode_pubkey(pk.h, pk.params.logn)

"Serialise a private key (omitting `G`, as the format does)."
privkey_bytes(sk::FalconPrivateKey) = encode_privkey(sk.f, sk.g, sk.F, sk.params.logn)
