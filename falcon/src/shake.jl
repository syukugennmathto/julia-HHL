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

using SHA

# ===========================================================================
# SHAKE256 as an XOF
# ===========================================================================
#
# The awkwardness here is worth explaining, because it is a real limitation of
# the tooling rather than of FALCON.
#
# Keccak-based XOFs are naturally *incremental*: you absorb the message, pad
# once, and then squeeze as many output blocks as you like, the state evolving
# as you go.  `SHA.jl` does not expose that: its `shake256(data, d)` takes the
# whole message and a fixed output length `d`, and `digest!` re-applies the
# padding every time it is called, so a context cannot be squeezed twice.
#
# What saves us is the defining property of a XOF stream: for `d1 < d2`,
# `SHAKE256(m, d1)` is a *prefix* of `SHAKE256(m, d2)`.  So an incremental
# squeeze can be emulated by keeping the absorbed message around and, whenever
# more output is wanted than we have, recomputing a longer digest from scratch.
# We grow the buffer by doubling, so a sequence of small squeezes totalling `N`
# bytes costs O(N) Keccak permutations amortised, not O(N^2) -- but each
# individual regeneration re-absorbs the message, so the true bound is
# O(N + |m| log N).  For FALCON's message sizes that is irrelevant.
#
# The prefix property is an assumption about `SHA.jl`, not just about FIPS 202,
# so `test/test_shake.jl` checks it explicitly against CPython's `hashlib`
# rather than taking it on faith.
#
# CONSTANT TIME: not an issue here.  SHAKE256 is applied to public data (the
# salt and the message), and the ChaCha20 PRNG below is a stream cipher whose
# reference implementation is already constant time by construction.  The
# side-channel difficulty in FALCON lives downstream, in `samplerz` and
# `ffsampling`, not here.

"""
    SHAKE256XOF

An incremental-squeeze view of SHAKE256, built on `SHA.shake256`.

Construct with [`shake256_xof`](@ref), then draw bytes with
[`squeeze!`](@ref).  Absorption is one-shot: the whole message must be given
at construction time.
"""
mutable struct SHAKE256XOF
    "the absorbed message"
    msg::Vector{UInt8}
    "the output stream generated so far"
    buf::Vector{UInt8}
    "number of bytes of `buf` already handed out"
    pos::Int
end

"""
    shake256_xof(msg...) -> SHAKE256XOF

Absorb the concatenation of `msg` (each element a byte vector or a `String`)
and return a XOF ready to be squeezed.

The variadic form exists because FALCON always hashes a *concatenation*
(typically `salt || message`), and building that concatenation at every call
site is where an off-by-one in the domain separation likes to hide.
"""
function shake256_xof(parts...)
    msg = UInt8[]
    for p in parts
        append!(msg, _as_bytes(p))
    end
    return SHAKE256XOF(msg, UInt8[], 0)
end

_as_bytes(x::AbstractVector{UInt8}) = x
_as_bytes(x::AbstractString) = codeunits(x)

"""
    squeeze!(xof, n) -> Vector{UInt8}

Return the next `n` bytes of the XOF output stream and advance the position.
"""
function squeeze!(xof::SHAKE256XOF, n::Integer)
    n >= 0 || throw(ArgumentError("cannot squeeze a negative number of bytes"))
    need = xof.pos + n
    if need > length(xof.buf)
        # Grow by doubling, with a floor so the very first squeeze does not
        # cost a regeneration immediately afterwards.
        newlen = max(need, 2 * length(xof.buf), 64)
        xof.buf = SHA.shake256(xof.msg, UInt(newlen))
    end
    out = xof.buf[(xof.pos + 1):(xof.pos + n)]
    xof.pos += n
    return out
end

"""
    shake256(data, outlen) -> Vector{UInt8}

One-shot SHAKE256 with `outlen` bytes of output.

A thin wrapper over `SHA.shake256` whose only jobs are to accept a plain
`Int` length (`SHA.shake256` insists on `UInt`) and to accept `String` input.
"""
shake256(data, outlen::Integer) = SHA.shake256(_as_bytes(data), UInt(outlen))

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
