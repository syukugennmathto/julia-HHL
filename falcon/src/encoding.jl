# encoding.jl -- serialisation of keys and signatures.
#
# ---------------------------------------------------------------------------
# Why this module is not boring
# ---------------------------------------------------------------------------
#
# Two reasons.
#
# **It is the only part that is fully pinned by the normative implementation.**
# Everything from module 5 onwards depends on floating point, and module 8
# measured what that costs: a 1-ulp difference in a constant table changes the
# signature (docs/debug_log.md #025).  Encoding does not touch a float.  So the
# C reference's bytes are reproducible here *exactly*, and
# `test/vectors/cref_kat.jl` holds a real FALCON-512 key and signature emitted
# by the C reference for us to match.  This is where interoperability is
# actually testable.
#
# **The signature length is variable, and the specification fixes it anyway.**
# Signature coefficients are Golomb-Rice coded: sign bit, seven low bits
# verbatim, high bits in unary.  How many bytes that takes depends on the
# coefficients, so a signature can *fail to fit* in the 666 bytes FALCON-512
# allots -- in which case signing must **discard it and sample again**.  That
# retry loop is part of the scheme, not an implementation detail, and skipping
# it produces an implementation that works until it doesn't.
#
# (Note what is *not* redrawn on a retry: the salt.  Both references keep it
# and re-run only the sampler -- the C code draws the nonce once in
# `falcon_sign_start` and copies the same 40 bytes inside the loop.  An earlier
# draft of this comment said "with a fresh salt", which was wrong; see
# falcon.jl's `falcon_sign` for why the distinction matters.)
#
# ---------------------------------------------------------------------------
# The formats, from the normative C codec
# ---------------------------------------------------------------------------
#
# Public key   (897 bytes at n=512): 0x00+logn, then h at 14 bits/coefficient
# Private key (1281 bytes at n=512): 0x50+logn, then f, g at `max_fg_bits[logn]`
#                                    bits and F at `max_FG_bits[logn]` bits.
#                                    G is *not* stored -- it is recomputed from
#                                    f, g, F and the NTRU equation.
# Signature    (666 bytes at n=512): 0x30+logn, 40-byte salt, compressed s2,
#                                    zero padding to the fixed length.
#
# All three header bytes were confirmed against the C reference's own output:
# 0x09, 0x59, 0x39 at logn = 9.
#
# CONSTANT TIME: the compressor's length depends on the coefficients, hence on
# the signature, hence on the secret -- but the signature is public once
# emitted, so this is one of the few places where variable time is genuinely
# harmless.  The *retry* loop is a different matter: how many attempts a
# signature took is a function of the key, and the reference does not try to
# hide it.

# ---------------------------------------------------------------------------
# Bit packing
# ---------------------------------------------------------------------------
#
# Both key codecs pack fixed-width fields most-significant-bit first into a
# byte stream.  The C code does this with a `uint32_t` accumulator; we do the
# same, because the accumulator width is observable: it is what makes the final
# partial byte left-aligned rather than right-aligned.

"""
    _pack_bits(values, bits) -> Vector{UInt8}

Pack `values` as `bits`-bit fields, MSB first, padding the final byte with
zeros on the right.

[C-ref] codec.c:36-73 (`modq_encode`), codec.c:212-254 (`trim_i8_encode`)
"""
function _pack_bits(values::AbstractVector{<:Integer}, bits::Int)
    out = UInt8[]
    acc = UInt32(0)
    acc_len = 0
    mask = (UInt32(1) << bits) - UInt32(1)
    for v in values
        acc = (acc << bits) | (UInt32(v & Int(mask)) & mask)
        acc_len += bits
        while acc_len >= 8
            acc_len -= 8
            push!(out, UInt8((acc >> acc_len) & 0xff))
        end
    end
    if acc_len > 0
        push!(out, UInt8((acc << (8 - acc_len)) & 0xff))
    end
    return out
end

"""
    _unpack_bits(bytes, n, bits) -> Vector{Int}

Inverse of [`_pack_bits`](@ref): read `n` fields of `bits` bits, MSB first.
Returns unsigned field values; the caller applies any sign convention.
"""
function _unpack_bits(bytes::AbstractVector{UInt8}, n::Int, bits::Int)
    out = Vector{Int}(undef, n)
    acc = UInt32(0)
    acc_len = 0
    mask = (UInt32(1) << bits) - UInt32(1)
    u = 1
    pos = 1
    while u <= n
        pos <= length(bytes) || throw(ArgumentError("truncated encoding"))
        acc = (acc << 8) | UInt32(bytes[pos]); pos += 1
        acc_len += 8
        while acc_len >= bits && u <= n
            acc_len -= bits
            out[u] = Int((acc >> acc_len) & mask)
            u += 1
        end
    end
    # The C decoders reject a non-zero remainder in the last byte, so that an
    # encoding is unique.  Without this check a signature or key could be
    # mutated without detection, which is a real (if minor) malleability issue.
    leftover = acc & ((UInt32(1) << acc_len) - UInt32(1))
    leftover == 0 || throw(ArgumentError(
        "non-zero padding bits in the last byte: the encoding is not canonical"))
    return out
end

# ---------------------------------------------------------------------------
# Key coefficient widths
# ---------------------------------------------------------------------------

"""
    MAX_FG_BITS

Bits per coefficient used to store `f` and `g`, indexed by `logn + 1`.

[C-ref] codec.c:507-518 (`Zf(max_fg_bits)`)

At `logn = 9` this is 6: key generation samples `f`, `g` at `sigma_fg ~ 4.05`,
so `|coefficient| < 32` with overwhelming probability, and the C reference
simply rejects a key with a larger one.
"""
const MAX_FG_BITS = Int[0, 8, 8, 8, 8, 8, 7, 7, 6, 6, 5]

"""
    MAX_FG_BITS_F

Bits per coefficient used to store `F`, indexed by `logn + 1`: 8 everywhere.

[C-ref] codec.c:520-532 (`Zf(max_FG_bits)`)

`F` is larger than `f` and `g` -- it comes out of the NTRU solve rather than
the Gaussian -- so it gets the full byte.  `G` is not stored at all: it is
recovered from `f*G - g*F = q`.
"""
const MAX_FG_BITS_F = Int[0, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

# ---------------------------------------------------------------------------
# Public key
# ---------------------------------------------------------------------------

"""
    encode_pubkey(h, logn) -> Vector{UInt8}

Serialise the public key `h` (coefficients in `[0, q)`) at 14 bits each,
behind a `0x00 + logn` header byte.

[C-ref] falcon.c:210 (header), codec.c:36-73 (`modq_encode`)
"""
function encode_pubkey(h::AbstractVector{<:Integer}, logn::Integer)
    n = 1 << logn
    length(h) == n || throw(ArgumentError("h has length $(length(h)), expected $n"))
    all(c -> 0 <= c < Q, h) || throw(ArgumentError("public key coefficients must be in [0, q)"))
    return vcat(UInt8[0x00 + UInt8(logn)], _pack_bits(h, 14))
end

"""
    decode_pubkey(bytes) -> (logn, h)

Inverse of [`encode_pubkey`](@ref).  Rejects a wrong header, a wrong length,
non-canonical padding, and any coefficient `>= q`.

[C-ref] codec.c:76-120 (`modq_decode`)
"""
function decode_pubkey(bytes::AbstractVector{UInt8})
    isempty(bytes) && throw(ArgumentError("empty public key"))
    head = bytes[1]
    (head & 0xf0) == 0x00 || throw(ArgumentError(
        "bad public key header byte 0x$(string(head, base=16, pad=2))"))
    logn = Int(head & 0x0f)
    1 <= logn <= 10 || throw(ArgumentError("bad degree in header: logn = $logn"))
    n = 1 << logn
    want = 1 + ((n * 14 + 7) >> 3)
    length(bytes) == want || throw(ArgumentError(
        "public key is $(length(bytes)) bytes, expected $want for logn = $logn"))
    h = _unpack_bits(@view(bytes[2:end]), n, 14)
    all(c -> c < Q, h) || throw(ArgumentError("public key coefficient >= q"))
    return (logn, h)
end

# ---------------------------------------------------------------------------
# Private key
# ---------------------------------------------------------------------------

"""
    encode_privkey(f, g, F, logn) -> Vector{UInt8}

Serialise the private key behind a `0x50 + logn` header byte: `f` and `g` at
`MAX_FG_BITS[logn+1]` bits, then `F` at `MAX_FG_BITS_F[logn+1]` bits.

`G` is deliberately not stored; `recover_G` reconstructs it.

[C-ref] falcon.c:175 (header), codec.c:212-254 (`trim_i8_encode`)
"""
function encode_privkey(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer},
                        F::AbstractVector{<:Integer}, logn::Integer)
    n = 1 << logn
    (length(f) == length(g) == length(F) == n) ||
        throw(ArgumentError("f, g, F must all have length $n"))
    bfg = MAX_FG_BITS[logn + 1]
    bF = MAX_FG_BITS_F[logn + 1]
    _check_trim_range(f, bfg, "f"); _check_trim_range(g, bfg, "g")
    _check_trim_range(F, bF, "F")
    return vcat(UInt8[0x50 + UInt8(logn)],
                _pack_bits(f, bfg), _pack_bits(g, bfg), _pack_bits(F, bF))
end

function _check_trim_range(x, bits::Int, name::AbstractString)
    maxv = (1 << (bits - 1)) - 1
    all(c -> -maxv <= c <= maxv, x) || throw(ArgumentError(
        "$name has a coefficient outside [-$maxv, $maxv], which does not fit " *
        "in $bits bits; the reference rejects such a key"))
    return nothing
end

"""
    decode_privkey(bytes) -> (logn, f, g, F)

Inverse of [`encode_privkey`](@ref).

The value `-2^(bits-1)` is **forbidden** by the C decoder even though it fits,
so that the encoding is a bijection; we reject it too.

[C-ref] codec.c:256-310 (`trim_i8_decode`)
"""
function decode_privkey(bytes::AbstractVector{UInt8})
    isempty(bytes) && throw(ArgumentError("empty private key"))
    head = bytes[1]
    (head & 0xf0) == 0x50 || throw(ArgumentError(
        "bad private key header byte 0x$(string(head, base=16, pad=2))"))
    logn = Int(head & 0x0f)
    1 <= logn <= 10 || throw(ArgumentError("bad degree in header: logn = $logn"))
    n = 1 << logn
    bfg = MAX_FG_BITS[logn + 1]
    bF = MAX_FG_BITS_F[logn + 1]
    lfg = (n * bfg + 7) >> 3
    lF = (n * bF + 7) >> 3
    want = 1 + 2 * lfg + lF
    length(bytes) == want || throw(ArgumentError(
        "private key is $(length(bytes)) bytes, expected $want for logn = $logn"))
    o = 2
    f = _signed_fields(@view(bytes[o:(o + lfg - 1)]), n, bfg); o += lfg
    g = _signed_fields(@view(bytes[o:(o + lfg - 1)]), n, bfg); o += lfg
    F = _signed_fields(@view(bytes[o:(o + lF - 1)]), n, bF)
    return (logn, f, g, F)
end

function _signed_fields(bytes, n::Int, bits::Int)
    raw = _unpack_bits(bytes, n, bits)
    half = 1 << (bits - 1)
    out = Vector{Int}(undef, n)
    for i in 1:n
        w = raw[i]
        v = w >= half ? w - (1 << bits) : w
        v == -half && throw(ArgumentError(
            "the value -$half is forbidden in a $bits-bit field (non-canonical)"))
        out[i] = v
    end
    return out
end

"""
    recover_G(f, g, F; q = Q) -> Vector{BigInt}

Recover `G` from `f*G - g*F = q`, i.e. `G = (q + g*F) / f` computed in `R_q`
and then lifted.

The private key format omits `G` because it is redundant.  Recovering it needs
a division in `R_q`, so it fails exactly when `f` is not invertible -- which
key generation already excluded.

[C-ref] inner.h:618 (`Zf(complete_private)`)
"""
function recover_G(f::AbstractVector{<:Integer}, g::AbstractVector{<:Integer},
                   F::AbstractVector{<:Integer}; q::Integer = Q)
    n = length(f)
    fq = Int[mod(c, q) for c in f]
    gq = Int[mod(c, q) for c in g]
    Fq = Int[mod(c, q) for c in F]
    is_invertible_zq(fq) || throw(ArgumentError("f is not invertible mod q"))
    rhs = polyaddq([mod(q, q); zeros(Int, n - 1)], polymulq(gq, Fq), q)
    Gq = polydivq(rhs, fq)
    # Lift to the centred representatives: G is small, like F.
    return BigInt.(centered(Gq, q))
end

# ---------------------------------------------------------------------------
# Signature compression (Golomb-Rice)
# ---------------------------------------------------------------------------

"""
    compress_sig(v, slen) -> Union{Vector{UInt8},Nothing}

Golomb-Rice encode the signature coefficients `v` into exactly `slen` bytes,
or return `nothing` if they do not fit.

Per coefficient: one sign bit, then the seven low bits of `|v|` verbatim, then
`|v| >> 7` zero bits followed by a terminating one.  Small coefficients cost
nine bits; large ones cost more, which is why the total length is data
dependent.

[Py-ref] scripts/pyref/encoding.py:6-33 (`compress`)

Returning `nothing` is not an error -- it is the signal that signing must
sample again (keeping the same salt).  See the module header.
"""
function compress_sig(v::AbstractVector{<:Integer}, slen::Integer)
    bits = BitVector()
    for coef in v
        a = abs(Int(coef))
        push!(bits, coef < 0)                       # sign
        for k in 6:-1:0                             # seven low bits, MSB first
            push!(bits, ((a >> k) & 1) == 1)
        end
        for _ in 1:(a >> 7)                         # high bits, unary
            push!(bits, false)
        end
        push!(bits, true)
    end
    length(bits) > 8 * slen && return nothing       # does not fit: retry signing
    while length(bits) < 8 * slen
        push!(bits, false)
    end
    out = zeros(UInt8, slen)
    for i in 1:length(bits)
        if bits[i]
            byte = (i - 1) >> 3 + 1
            out[byte] |= UInt8(0x80 >> ((i - 1) & 7))
        end
    end
    return out
end

"""
    decompress_sig(x, slen, n) -> Union{Vector{Int},Nothing}

Inverse of [`compress_sig`](@ref); returns `nothing` if `x` is not a valid
encoding of exactly `n` coefficients.

[Py-ref] scripts/pyref/encoding.py:36-79 (`decompress`)

Two validity rules beyond "the bits run out", both there to make the encoding
canonical:

  * `-0` is forbidden (sign bit set with all magnitude bits zero);
  * after the `n`-th coefficient every remaining bit must be zero.

Without them a signature could be altered without becoming invalid.
"""
function decompress_sig(x::AbstractVector{UInt8}, slen::Integer, n::Integer)
    length(x) > slen && return nothing
    total = 8 * length(x)
    bit(i) = ((x[(i - 1) >> 3 + 1] >> (7 - ((i - 1) & 7))) & 1) == 1

    v = Int[]
    i = 1
    while length(v) < n
        i + 7 <= total || return nothing
        neg = bit(i)
        low = 0
        for k in 1:7
            low = (low << 1) | (bit(i + k) ? 1 : 0)
        end
        i += 8
        high = 0
        while true
            i <= total || return nothing
            bit(i) && break
            high += 1
            i += 1
            high > 2040 && return nothing        # runaway unary run
        end
        i += 1                                   # consume the terminating 1
        coef = low + (high << 7)
        (coef == 0 && neg) && return nothing     # -0 is not a valid encoding
        push!(v, neg ? -coef : coef)
    end
    # every remaining bit must be zero
    while i <= total
        bit(i) && return nothing
        i += 1
    end
    return v
end

# ---------------------------------------------------------------------------
# Signature container
# ---------------------------------------------------------------------------

"""
    encode_signature(salt, s2, logn, sig_bytes) -> Union{Vector{UInt8},Nothing}

Assemble a PADDED-format signature: `0x30 + logn`, the 40-byte salt, then the
compressed `s2`, zero-padded to exactly `sig_bytes`.

Returns `nothing` when `s2` does not compress into the available room, which
is the signal to sample again (keeping the same salt).

[C-ref] falcon.c:464 (header byte), falcon.h:335 (padded size)
"""
function encode_signature(salt::AbstractVector{UInt8}, s2::AbstractVector{<:Integer},
                          logn::Integer, sig_bytes::Integer)
    length(salt) == SALT_LEN ||
        throw(ArgumentError("salt must be $SALT_LEN bytes, got $(length(salt))"))
    room = sig_bytes - HEAD_LEN - SALT_LEN
    room > 0 || throw(ArgumentError("sig_bytes = $sig_bytes leaves no room"))
    body = compress_sig(s2, room)
    body === nothing && return nothing
    return vcat(UInt8[0x30 + UInt8(logn)], collect(salt), body)
end

"""
    decode_signature(bytes) -> (logn, salt, s2)

Inverse of [`encode_signature`](@ref).  Throws on a malformed signature rather
than returning a wrong one.
"""
function decode_signature(bytes::AbstractVector{UInt8})
    length(bytes) > HEAD_LEN + SALT_LEN ||
        throw(ArgumentError("signature is too short: $(length(bytes)) bytes"))
    head = bytes[1]
    (head & 0xf0) == 0x30 || throw(ArgumentError(
        "bad signature header byte 0x$(string(head, base=16, pad=2))"))
    logn = Int(head & 0x0f)
    1 <= logn <= 10 || throw(ArgumentError("bad degree in header: logn = $logn"))
    n = 1 << logn
    salt = collect(bytes[2:(1 + SALT_LEN)])
    body = @view bytes[(2 + SALT_LEN):end]
    s2 = decompress_sig(collect(body), length(body), n)
    s2 === nothing && throw(ArgumentError("malformed compressed signature"))
    return (logn, salt, s2)
end
