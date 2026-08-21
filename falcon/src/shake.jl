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
    @inbounds begin
    # The 25 lanes are pulled into locals for the whole permutation, so the
    # rounds run in registers.  The previous version kept them in the array and
    # allocated two temporaries per call; more importantly it addressed the
    # rho/pi step through `_KECCAK_PI[i]` and computed `mod1(x + 1, 5)` inside
    # chi, i.e. a table indirection and an integer division per lane per round.
    # That cost 1.11 us per permutation against the C reference's 0.49; this
    # form is generated from the same two tables, so the constants below are
    # not transcribed by hand (docs/debug_log.md #038).
        a1 = A[1]
        a2 = A[2]
        a3 = A[3]
        a4 = A[4]
        a5 = A[5]
        a6 = A[6]
        a7 = A[7]
        a8 = A[8]
        a9 = A[9]
        a10 = A[10]
        a11 = A[11]
        a12 = A[12]
        a13 = A[13]
        a14 = A[14]
        a15 = A[15]
        a16 = A[16]
        a17 = A[17]
        a18 = A[18]
        a19 = A[19]
        a20 = A[20]
        a21 = A[21]
        a22 = A[22]
        a23 = A[23]
        a24 = A[24]
        a25 = A[25]

    for round in 1:24
        # --- theta ---------------------------------------------------
        c1 = a1 ⊻ a6 ⊻ a11 ⊻ a16 ⊻ a21
        c2 = a2 ⊻ a7 ⊻ a12 ⊻ a17 ⊻ a22
        c3 = a3 ⊻ a8 ⊻ a13 ⊻ a18 ⊻ a23
        c4 = a4 ⊻ a9 ⊻ a14 ⊻ a19 ⊻ a24
        c5 = a5 ⊻ a10 ⊻ a15 ⊻ a20 ⊻ a25
        d1 = c5 ⊻ bitrotate(c2, 1)
        d2 = c1 ⊻ bitrotate(c3, 1)
        d3 = c2 ⊻ bitrotate(c4, 1)
        d4 = c3 ⊻ bitrotate(c5, 1)
        d5 = c4 ⊻ bitrotate(c1, 1)
        a1 ⊻= d1; a6 ⊻= d1; a11 ⊻= d1; a16 ⊻= d1; a21 ⊻= d1
        a2 ⊻= d2; a7 ⊻= d2; a12 ⊻= d2; a17 ⊻= d2; a22 ⊻= d2
        a3 ⊻= d3; a8 ⊻= d3; a13 ⊻= d3; a18 ⊻= d3; a23 ⊻= d3
        a4 ⊻= d4; a9 ⊻= d4; a14 ⊻= d4; a19 ⊻= d4; a24 ⊻= d4
        a5 ⊻= d5; a10 ⊻= d5; a15 ⊻= d5; a20 ⊻= d5; a25 ⊻= d5

        # --- rho and pi ---------------------------------------------
        b1 = a1
        b2 = bitrotate(a7, 44)
        b3 = bitrotate(a13, 43)
        b4 = bitrotate(a19, 21)
        b5 = bitrotate(a25, 14)
        b6 = bitrotate(a4, 28)
        b7 = bitrotate(a10, 20)
        b8 = bitrotate(a11, 3)
        b9 = bitrotate(a17, 45)
        b10 = bitrotate(a23, 61)
        b11 = bitrotate(a2, 1)
        b12 = bitrotate(a8, 6)
        b13 = bitrotate(a14, 25)
        b14 = bitrotate(a20, 8)
        b15 = bitrotate(a21, 18)
        b16 = bitrotate(a5, 27)
        b17 = bitrotate(a6, 36)
        b18 = bitrotate(a12, 10)
        b19 = bitrotate(a18, 15)
        b20 = bitrotate(a24, 56)
        b21 = bitrotate(a3, 62)
        b22 = bitrotate(a9, 55)
        b23 = bitrotate(a15, 39)
        b24 = bitrotate(a16, 41)
        b25 = bitrotate(a22, 2)

        # --- chi ----------------------------------------------------
        a1 = b1 ⊻ (~b2 & b3)
        a2 = b2 ⊻ (~b3 & b4)
        a3 = b3 ⊻ (~b4 & b5)
        a4 = b4 ⊻ (~b5 & b1)
        a5 = b5 ⊻ (~b1 & b2)
        a6 = b6 ⊻ (~b7 & b8)
        a7 = b7 ⊻ (~b8 & b9)
        a8 = b8 ⊻ (~b9 & b10)
        a9 = b9 ⊻ (~b10 & b6)
        a10 = b10 ⊻ (~b6 & b7)
        a11 = b11 ⊻ (~b12 & b13)
        a12 = b12 ⊻ (~b13 & b14)
        a13 = b13 ⊻ (~b14 & b15)
        a14 = b14 ⊻ (~b15 & b11)
        a15 = b15 ⊻ (~b11 & b12)
        a16 = b16 ⊻ (~b17 & b18)
        a17 = b17 ⊻ (~b18 & b19)
        a18 = b18 ⊻ (~b19 & b20)
        a19 = b19 ⊻ (~b20 & b16)
        a20 = b20 ⊻ (~b16 & b17)
        a21 = b21 ⊻ (~b22 & b23)
        a22 = b22 ⊻ (~b23 & b24)
        a23 = b23 ⊻ (~b24 & b25)
        a24 = b24 ⊻ (~b25 & b21)
        a25 = b25 ⊻ (~b21 & b22)

        # --- iota ---------------------------------------------------
        a1 ⊻= _KECCAK_RC[round]
    end

        A[1] = a1
        A[2] = a2
        A[3] = a3
        A[4] = a4
        A[5] = a5
        A[6] = a6
        A[7] = a7
        A[8] = a8
        A[9] = a9
        A[10] = a10
        A[11] = a11
        A[12] = a12
        A[13] = a13
        A[14] = a14
        A[15] = a15
        A[16] = a16
        A[17] = a17
        A[18] = a18
        A[19] = a19
        A[20] = a20
        A[21] = a21
        A[22] = a22
        A[23] = a23
        A[24] = a24
        A[25] = a25
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
    k = 1
    while k <= n
        if x.pos == SHAKE256_RATE
            _keccak_f1600!(x.state)
            x.pos = 0
        end
        m = min(SHAKE256_RATE - x.pos, n - k + 1)
        _copy_state_bytes!(out, k, x.state, x.pos, m)
        x.pos += m
        k += m
    end
    return out
end

"""
Copy `m` bytes of sponge state, starting at byte offset `pos`, into `out[k...]`.

On a little-endian machine the byte order of the `UInt64` lanes *is* the
sponge's byte order, so this is one `memcpy`; extracting a byte at a time cost
a lane index, a shift and a mask per byte, and squeezing is on the critical
path of `hash_to_point` (docs/debug_log.md #038).

The byte-at-a-time path is kept for big-endian machines, where the lane layout
and the sponge layout disagree and the memcpy would be wrong.  It is not dead
code that cannot be reached -- it is the correct implementation on hardware
this project has not been run on.
"""
@inline function _copy_state_bytes!(out::Vector{UInt8}, k::Int,
                                    state::Vector{UInt64}, pos::Int, m::Int)
    if Base.ENDIAN_BOM == 0x04030201                     # little-endian
        GC.@preserve out state begin
            unsafe_copyto!(pointer(out, k), Ptr{UInt8}(pointer(state)) + pos, m)
        end
    else
        @inbounds for t in 0:(m - 1)
            out[k + t] = _get_byte(state, pos + t)
        end
    end
    return nothing
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
