# shake.jl -- the two pseudorandom primitives FALCON needs.
#
#   1. SHAKE256 as an extendable-output function (XOF), used to derive the
#      "hashed message" point and to expand seeds.
#   2. The ChaCha20-based PRNG of the reference implementation, used for the
#      randomness consumed by key generation and by signing.
#
# Why both?  SHAKE256 is the *specified* hash; the ChaCha20 PRNG is the
# *specified* way of turning a 56-byte seed into the long stream that
# `samplerz` consumes.  KAT agreement requires reproducing the second one
# byte-for-byte, including its slightly surprising buffering behaviour.

# ---------------------------------------------------------------------------
# SHAKE256, implemented here rather than taken from a library
# ---------------------------------------------------------------------------
#
# The first draft of this module used `SHA.shake256`.  That was a mistake made
# by reading SHA.jl's *master branch* on GitHub instead of the version actually
# shipped as a stdlib: SHA v0.7.0, which is what Julia 1.11 bundles, exports no
# SHAKE at all (`isdefined(SHA, :shake256)` is `false`).  See
# docs/debug_log.md #014.
#
# Writing Keccak ourselves costs about ninety lines and buys three things:
#
#   1. it works on whatever Julia the reader has, with no dependency whose API
#      might move underneath us;
#   2. it gives a *real* incremental squeeze, so the awkward
#      "regenerate a longer digest and rely on the prefix property" workaround
#      of the first draft (debug_log #004) disappears entirely;
#   3. the round constants and rotation offsets can be **generated from their
#      definitions** rather than transcribed from a table.  That is the same
#      discipline params.jl applies to the FALCON constants, and it matters
#      more here: a mistyped rotation offset produces a hash that is wrong but
#      looks perfectly random, so only a KAT would ever catch it.
#
# CONSTANT TIME: Keccak is branch-free and table-free by construction, so this
# is constant time without effort.  As in the ChaCha20 note below, the
# side-channel difficulty in FALCON is downstream, in samplerz and ffsampling.

"The rate of SHAKE256 in bytes: 200 - 2*32 = 136."
const SHAKE256_RATE = 136

"""
    _keccak_round_constants() -> Vector{UInt64}

The 24 round constants of Keccak-f[1600], generated from the LFSR definition
in FIPS 202 rather than copied from a table.

`rc(t)` is bit `t mod 255` of the sequence produced by the LFSR with feedback
polynomial `x^8 + x^6 + x^5 + x^4 + 1`; round constant `i` has bit
`2^j - 1` set to `rc(j + 7i)` for `j = 0..6`.
"""
function _keccak_round_constants()
    rc = Vector{UInt64}(undef, 24)
    # NOTE the type: this must be wider than 8 bits.  Written as `0x01` it is a
    # UInt8, `lfsr <<= 1` silently drops the bit we are about to test for, the
    # feedback never fires, and every round constant comes out wrong -- while
    # still looking like plausible random-ish values.  debug_log #015.
    lfsr = 1
    for i in 1:24
        c = UInt64(0)
        for j in 0:6
            bit = lfsr & 0x01
            # advance: shift left, and reduce by the feedback polynomial
            lfsr <<= 1
            if (lfsr & 0x100) != 0
                lfsr = (lfsr ⊻ 0x71) & 0xff
            end
            if bit != 0
                c |= UInt64(1) << ((1 << j) - 1)
            end
        end
        rc[i] = c
    end
    return rc
end

"""
    _keccak_rho_pi() -> (Vector{Int}, Vector{Int})

Rotation offsets and the pi permutation, generated from their definitions.

Lanes are indexed `x + 5y + 1` (1-based).  Starting from `(x, y) = (1, 0)`, the
walk `(x, y) <- (y, 2x + 3y mod 5)` visits the 24 non-origin lanes, and at step
`t` the offset is the triangular number `(t+1)(t+2)/2 mod 64`.  The same walk
defines pi.
"""
function _keccak_rho_pi()
    rot = zeros(Int, 25)
    x, y = 1, 0
    for t in 0:23
        rot[x + 5y + 1] = ((t + 1) * (t + 2) ÷ 2) % 64
        x, y = y, (2x + 3y) % 5
    end
    # pi: the lane at (x, y) moves to (y, 2x + 3y).
    piperm = zeros(Int, 25)
    for yy in 0:4, xx in 0:4
        piperm[yy + 5 * ((2xx + 3yy) % 5) + 1] = xx + 5yy + 1
    end
    return (rot, piperm)
end

const _KECCAK_RC = _keccak_round_constants()
const _KECCAK_ROT, _KECCAK_PI = _keccak_rho_pi()

"""
    _keccak_f1600!(A)

The Keccak-f[1600] permutation, in place, on 25 `UInt64` lanes.
"""
function _keccak_f1600!(A::Vector{UInt64})
    C = Vector{UInt64}(undef, 5)
    B = Vector{UInt64}(undef, 25)
    @inbounds for round in 1:24
        # theta
        for x in 1:5
            C[x] = A[x] ⊻ A[x + 5] ⊻ A[x + 10] ⊻ A[x + 15] ⊻ A[x + 20]
        end
        for x in 1:5
            d = C[mod1(x - 1, 5)] ⊻ bitrotate(C[mod1(x + 1, 5)], 1)
            for y in 0:4
                A[x + 5y] ⊻= d
            end
        end
        # rho and pi
        for i in 1:25
            B[i] = bitrotate(A[_KECCAK_PI[i]], _KECCAK_ROT[_KECCAK_PI[i]])
        end
        # chi
        for y in 0:4, x in 1:5
            A[x + 5y] = B[x + 5y] ⊻ ((~B[mod1(x + 1, 5) + 5y]) & B[mod1(x + 2, 5) + 5y])
        end
        # iota
        A[1] ⊻= _KECCAK_RC[round]
    end
    return A
end

"""
    SHAKE256XOF

An incremental SHAKE256 sponge: absorb with [`absorb!`](@ref), then draw output
with [`squeeze!`](@ref).  Once squeezing has begun, further absorption is an
error -- padding has already been applied and the sponge cannot go back.
"""
mutable struct SHAKE256XOF
    "the 1600-bit state, as 25 lanes"
    state::Vector{UInt64}
    "bytes absorbed into the current block but not yet permuted"
    buf::Vector{UInt8}
    "how many bytes of `buf` are in use (absorbing) or consumed (squeezing)"
    pos::Int
    "false while absorbing, true once padded"
    squeezing::Bool
end

SHAKE256XOF() = SHAKE256XOF(zeros(UInt64, 25), zeros(UInt8, SHAKE256_RATE), 0, false)

"""
    shake256_xof(parts...) -> SHAKE256XOF

Absorb the concatenation of `parts` (byte vectors or `String`s) and return a
sponge ready to be squeezed.

The variadic form exists because FALCON always hashes a *concatenation*
(typically `salt || message`), and building that concatenation at every call
site is where an off-by-one in the domain separation likes to hide.
"""
function shake256_xof(parts...)
    x = SHAKE256XOF()
    for p in parts
        absorb!(x, _as_bytes(p))
    end
    return x
end

_as_bytes(x::AbstractVector{UInt8}) = x
_as_bytes(x::AbstractString) = codeunits(x)

# XOR one byte into the sponge state at byte offset `i` (0-based).
@inline function _xor_byte!(state::Vector{UInt64}, i::Int, b::UInt8)
    lane = (i >> 3) + 1
    shift = (i & 7) << 3
    @inbounds state[lane] ⊻= UInt64(b) << shift
    return nothing
end

@inline function _get_byte(state::Vector{UInt64}, i::Int)
    lane = (i >> 3) + 1
    shift = (i & 7) << 3
    @inbounds return UInt8((state[lane] >> shift) & 0xff)
end

"""
    absorb!(xof, data)

Absorb more input.  Errors if the sponge has already been squeezed.
"""
function absorb!(x::SHAKE256XOF, data)
    x.squeezing && throw(ArgumentError(
        "cannot absorb into a SHAKE256 sponge that has already been squeezed"))
    bytes = _as_bytes(data)
    for b in bytes
        _xor_byte!(x.state, x.pos, b)
        x.pos += 1
        if x.pos == SHAKE256_RATE
            _keccak_f1600!(x.state)
            x.pos = 0
        end
    end
    return x
end

# Apply the SHAKE padding (0x1f ... 0x80) and switch to squeezing.
function _finalize!(x::SHAKE256XOF)
    _xor_byte!(x.state, x.pos, 0x1f)
    _xor_byte!(x.state, SHAKE256_RATE - 1, 0x80)
    _keccak_f1600!(x.state)
    x.squeezing = true
    x.pos = 0
    return nothing
end

"""
    squeeze!(xof, n) -> Vector{UInt8}

Return the next `n` bytes of output, permuting the state as needed.

Unlike the first draft of this module, this is a genuine incremental squeeze:
the sponge state advances, nothing is recomputed, and squeezing `n` bytes in
any sequence of chunks gives the same stream.
"""
function squeeze!(x::SHAKE256XOF, n::Integer)
    n >= 0 || throw(ArgumentError("cannot squeeze a negative number of bytes"))
    x.squeezing || _finalize!(x)
    out = Vector{UInt8}(undef, n)
    for k in 1:n
        if x.pos == SHAKE256_RATE
            _keccak_f1600!(x.state)
            x.pos = 0
        end
        out[k] = _get_byte(x.state, x.pos)
        x.pos += 1
    end
    return out
end

"""
    shake256(data, outlen) -> Vector{UInt8}

One-shot SHAKE256 with `outlen` bytes of output.
"""
shake256(data, outlen::Integer) = squeeze!(shake256_xof(data), outlen)

# ===========================================================================
# The ChaCha20 PRNG of the reference implementation
# ===========================================================================
#
# This is *not* RFC 7539 ChaCha20.  It is ChaCha20's block function wrapped in
# the reference implementation's own key/counter layout and output interleaving,
# and the differences are exactly the kind of thing that produces a
# "everything is individually correct but the KAT still fails" afternoon.  The
# three deviations to keep in mind:
#
#   (a) The 56-byte seed becomes 14 little-endian 32-bit words s[0..13].  Words
#       s[0..9] go into the state where RFC 7539 puts the key and the nonce;
#       s[10], s[11] are XORed with the 64-bit counter; s[12], s[13] form the
#       *initial value* of that counter.  So the seed simultaneously keys the
#       cipher and initialises the counter.
#
#   (b) Output is produced eight blocks at a time and the eight blocks are
#       *interleaved word-wise*: word j of block i lands at position 8j + i of
#       a 128-word buffer.  The reference does this because its assembly
#       generates eight blocks in parallel in SIMD registers; the interleaving
#       is what falls out of the register layout, and it became normative.
#
#   (c) `randombytes(k)` does **not** refill-and-continue when the buffer runs
#       short.  It discards whatever is left and regenerates a fresh 512-byte
#       buffer.  A request that straddles the boundary therefore skips bytes.
#       Getting this wrong desynchronises the stream only *sometimes*, which
#       makes it a nasty bug; `test/vectors/chacha20_kat.jl` has cases chosen
#       to straddle the boundary on purpose.
#
# Oracle for all of the above: scripts/pyref/rng.py, itself a transliteration
# of the C reference PRNG.
#
# CONSTANT TIME: ChaCha20 is arithmetic-only and branch-free, so this is one of
# the few parts of FALCON that is constant time without any effort.  Note the
# contrast with what comes later: the *consumer* of this stream (`samplerz`)
# is where the timing leak lives, not the producer.

"""
    ChaCha20

The FALCON reference PRNG.  Create with [`chacha20`](@ref) and draw bytes with
[`randombytes!`](@ref).
"""
mutable struct ChaCha20
    "the 14 seed words s[0..13]"
    s::Vector{UInt32}
    "the 64-bit block counter"
    ctr::UInt64
    "512-byte output buffer (8 interleaved ChaCha20 blocks)"
    buf::Vector{UInt8}
    "number of bytes of `buf` already consumed"
    pos::Int
end

"ChaCha20's constants, \"expand 32-byte k\".  [Py-ref] scripts/pyref/rng.py:22"
const CHACHA_CW = UInt32[0x61707865, 0x3320646e, 0x79622d32, 0x6b206574]

"""
    chacha20(seed) -> ChaCha20

Initialise the PRNG from a `SEED_LEN`-byte (56-byte) seed.

[Py-ref] scripts/pyref/rng.py:36-42
"""
function chacha20(seed::AbstractVector{UInt8})
    length(seed) == SEED_LEN ||
        throw(ArgumentError("ChaCha20 seed must be $SEED_LEN bytes, got $(length(seed))"))
    # Assemble the words by hand rather than via `reinterpret`: the reference
    # says "little-endian", and spelling it out keeps this correct on a
    # big-endian host instead of silently depending on the platform.
    s = UInt32[_le32(seed, 4i) for i in 0:13]
    ctr = UInt64(s[13]) | (UInt64(s[14]) << 32)   # s[12], s[13] zero-indexed
    return ChaCha20(s, ctr, UInt8[], 0)
end

@inline _le32(b::AbstractVector{UInt8}, off::Int) =
    UInt32(b[off + 1]) | (UInt32(b[off + 2]) << 8) |
    (UInt32(b[off + 3]) << 16) | (UInt32(b[off + 4]) << 24)

@inline _rotl32(x::UInt32, n::Int) = (x << n) | (x >> (32 - n))

# The ChaCha20 quarter-round, operating in place on a 16-word state.
# Indices are 1-based here and 0-based in the reference; the mapping is +1.
@inline function _qround!(st::Vector{UInt32}, ai::Int, bi::Int, ci::Int, di::Int)
    a, b, c, d = st[ai], st[bi], st[ci], st[di]
    a += b; d = _rotl32(d ⊻ a, 16)
    c += d; b = _rotl32(b ⊻ c, 12)
    a += b; d = _rotl32(d ⊻ a, 8)
    c += d; b = _rotl32(b ⊻ c, 7)
    st[ai], st[bi], st[ci], st[di] = a, b, c, d
    return nothing
end

"""
    _block!(rng) -> Vector{UInt32}

One ChaCha20 block: build the state, run 20 rounds (10 double rounds), add the
original state, increment the counter, return the 16 output words.

[Py-ref] scripts/pyref/rng.py:77-100 (`update`)
"""
function _block!(rng::ChaCha20)
    st = Vector{UInt32}(undef, 16)
    st[1:4] = CHACHA_CW
    st[5:14] = rng.s[1:10]                                   # s[0..9]
    st[15] = rng.s[11] ⊻ UInt32(rng.ctr & 0xffffffff)        # s[10] ^ lo(ctr)
    st[16] = rng.s[12] ⊻ UInt32((rng.ctr >> 32) & 0xffffffff) # s[11] ^ hi(ctr)
    orig = copy(st)
    for _ in 1:10
        # column rounds
        _qround!(st, 1, 5, 9, 13)
        _qround!(st, 2, 6, 10, 14)
        _qround!(st, 3, 7, 11, 15)
        _qround!(st, 4, 8, 12, 16)
        # diagonal rounds
        _qround!(st, 1, 6, 11, 16)
        _qround!(st, 2, 7, 12, 13)
        _qround!(st, 3, 8, 9, 14)
        _qround!(st, 4, 5, 10, 15)
    end
    for i in 1:16
        st[i] += orig[i]
    end
    rng.ctr += 1
    return st
end

"""
    _refill!(rng)

Generate eight blocks and interleave them word-wise into the 512-byte buffer:
word `j` of block `i` goes to buffer word `8j + i` (zero-indexed).

[Py-ref] scripts/pyref/rng.py:102-109 (`block_update`)
"""
function _refill!(rng::ChaCha20)
    words = Vector{UInt32}(undef, 128)
    for i in 0:7
        blk = _block!(rng)
        for j in 0:15
            words[8j + i + 1] = blk[j + 1]
        end
    end
    buf = Vector{UInt8}(undef, 512)
    for (k, w) in enumerate(words)
        off = 4 * (k - 1)
        buf[off + 1] = UInt8(w & 0xff)
        buf[off + 2] = UInt8((w >> 8) & 0xff)
        buf[off + 3] = UInt8((w >> 16) & 0xff)
        buf[off + 4] = UInt8((w >> 24) & 0xff)
    end
    rng.buf = buf
    rng.pos = 0
    return nothing
end

"""
    randombytes!(rng, k) -> Vector{UInt8}

Draw `k` pseudorandom bytes.

Reproduces the reference's buffering exactly, including the discard: if fewer
than `k` bytes remain in the current 512-byte buffer, the remainder is
**thrown away** and a fresh buffer is generated.

[Py-ref] scripts/pyref/rng.py:111-122 (`randombytes`)

(The reference's `randombytes` also performs a pair of byte reversals that
cancel exactly -- `bytes.fromhex(reversed_hex)[::-1]` is the identity on the
selected slice.  They are an artefact of matching the C code's word-oriented
consumption and have no observable effect, so they are not reproduced here.
`test/vectors/chacha20_kat.jl` is what certifies that claim.)
"""
function randombytes!(rng::ChaCha20, k::Integer)
    k >= 0 || throw(ArgumentError("k must be non-negative"))
    k > 512 && throw(ArgumentError(
        "the reference PRNG cannot serve a request larger than its 512-byte " *
        "buffer in one call (asked for $k)"))
    if rng.pos + k > length(rng.buf)
        _refill!(rng)
    end
    out = rng.buf[(rng.pos + 1):(rng.pos + k)]
    rng.pos += k
    return out
end
