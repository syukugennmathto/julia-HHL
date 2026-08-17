# test_encoding.jl -- module 9.
#
# The one module whose bytes are pinned by the *normative* implementation.
# Everything from module 5 onwards depends on floating point, and module 8
# measured what that costs (docs/debug_log.md #025); encoding touches no float,
# so the C reference's bytes are reproducible here exactly.
#
# test/vectors/cref_kat.jl holds a real FALCON-512 key and signature emitted by
# the C reference from a fixed seed, and verified by it. Matching those bytes
# is genuine interoperability rather than agreement with a transliteration.

# Golden vectors are included by runtests.jl.

@testset "encoding" begin

    @testset "the C reference's public key" begin
        logn, h = decode_pubkey(CREF_PUBKEY)
        @test logn == CREF_LOGN == 9
        @test length(h) == 512
        @test all(c -> 0 <= c < Q, h)
        # and we reproduce its bytes exactly
        @test encode_pubkey(h, logn) == CREF_PUBKEY
        @test length(CREF_PUBKEY) == FALCON_512.pubkey_bytes
        @test CREF_PUBKEY[1] == 0x00 + UInt8(logn)
    end

    @testset "the C reference's private key" begin
        logn, f, g, F = decode_privkey(CREF_PRIVKEY)
        @test logn == 9
        @test length(f) == length(g) == length(F) == 512
        @test encode_privkey(f, g, F, logn) == CREF_PRIVKEY
        @test length(CREF_PRIVKEY) == FALCON_512.privkey_bytes
        @test CREF_PRIVKEY[1] == 0x50 + UInt8(logn)

        # the widths the format allots really are enough for this key
        @test maximum(abs, f) <= (1 << (MAX_FG_BITS[logn + 1] - 1)) - 1
        @test maximum(abs, g) <= (1 << (MAX_FG_BITS[logn + 1] - 1)) - 1
        @test maximum(abs, F) <= (1 << (MAX_FG_BITS_F[logn + 1] - 1)) - 1
    end

    @testset "G is recoverable, and the key is consistent" begin
        # The private key format omits G. Recovering it and checking the NTRU
        # equation is a complete end-to-end validation of a real key produced by
        # the normative implementation -- exactly, in BigInt.
        _, f, g, F = decode_privkey(CREF_PRIVKEY)
        G = recover_G(f, g, F)
        @test length(G) == 512
        @test ntru_equation_holds(BigInt.(f), BigInt.(g), BigInt.(F), G)

        # and the public key really is g/f mod q
        _, h = decode_pubkey(CREF_PUBKEY)
        fq = Int[mod(c, Q) for c in f]
        gq = Int[mod(c, Q) for c in g]
        @test polydivq(gq, fq) == h
        @test polymulq(h, fq) == gq

        # the Gram-Schmidt bound key generation enforced
        @test gs_norm_ok(Float64.(f), Float64.(g))
    end

    @testset "the C reference's signature" begin
        logn, salt, s2 = decode_signature(CREF_SIGNATURE)
        @test logn == 9
        @test length(salt) == SALT_LEN == 40
        @test length(s2) == 512
        @test encode_signature(salt, s2, logn, length(CREF_SIGNATURE)) == CREF_SIGNATURE
        @test length(CREF_SIGNATURE) == FALCON_512.sig_bytes
        @test CREF_SIGNATURE[1] == 0x30 + UInt8(logn)
    end

    @testset "compression against the reference" begin
        for (v, slen, want) in COMPRESS_KAT
            got = compress_sig(v, slen)
            if want === nothing
                @test got === nothing        # does not fit: signing must retry
            else
                @test got == want
                @test length(got) == slen
                # and it decodes back
                @test decompress_sig(got, slen, length(v)) == v
            end
        end
        # at least one recorded case must be a non-fit, or the retry path is
        # untested
        @test any(t -> t[3] === nothing, COMPRESS_KAT)
    end

    @testset "compression round-trips" begin
        for (v, slen, want) in COMPRESS_KAT
            want === nothing && continue
            @test decompress_sig(compress_sig(v, slen), slen, length(v)) == v
        end
        # and on random signature-shaped data
        rng = chacha20(collect(UInt8, 0:55))
        for _ in 1:50
            v = [Int(b) - 128 for b in randombytes!(rng, 64)]
            enc = compress_sig(v, 200)
            @test enc !== nothing
            @test decompress_sig(enc, 200, 64) == v
        end
    end

    @testset "invalid encodings are rejected, not crashed on" begin
        for (x, slen, n) in DECOMPRESS_INVALID
            @test decompress_sig(x, slen, n) === nothing
        end
        # The all-zero encoding in particular. The Python reference *raises
        # IndexError* here rather than rejecting it -- its trailing-zero strip
        # runs off the front of the string, outside the try/except meant to
        # catch that (docs/debug_log.md #027). A decoder fed attacker-supplied
        # bytes must not crash.
        @test decompress_sig(zeros(UInt8, 8), 8, 4) === nothing
        @test decompress_sig(UInt8[], 0, 1) === nothing
        # too long for the declared length
        @test decompress_sig(zeros(UInt8, 10), 4, 1) === nothing
    end

    @testset "malformed containers are rejected" begin
        bad_head = copy(CREF_PUBKEY); bad_head[1] = 0x99
        @test_throws ArgumentError decode_pubkey(bad_head)
        @test_throws ArgumentError decode_pubkey(CREF_PUBKEY[1:end-1])
        @test_throws ArgumentError decode_pubkey(UInt8[])

        bad_sk = copy(CREF_PRIVKEY); bad_sk[1] = 0x00
        @test_throws ArgumentError decode_privkey(bad_sk)
        @test_throws ArgumentError decode_privkey(CREF_PRIVKEY[1:end-1])

        bad_sig = copy(CREF_SIGNATURE); bad_sig[1] = 0x00
        @test_throws ArgumentError decode_signature(bad_sig)
        @test_throws ArgumentError decode_signature(UInt8[0x39])
    end

    @testset "encodings are canonical" begin
        # Non-zero padding bits in the final byte must be rejected, or a key
        # could be mutated without becoming invalid.
        #
        # NOTE which degree: at logn = 9 there ARE no padding bits -- 512*14 =
        # 7168 bits is exactly 896 bytes, and likewise for f, g (512*6) and F
        # (512*8).  The first version of this test corrupted the last byte of
        # the FALCON-512 public key and expected a rejection; it got a
        # different but perfectly valid key instead.  Padding exists only when
        # n*bits is not a multiple of 8, e.g. logn = 1: 2*14 = 28 bits in 4
        # bytes, so four padding bits.  (docs/debug_log.md #028.)
        @test (512 * 14) % 8 == 0                      # why the above matters
        small = encode_pubkey([1, 2], 1)
        @test length(small) == 1 + 4
        @test decode_pubkey(small) == (1, [1, 2])
        mutated = copy(small)
        mutated[end] |= 0x01                            # set a padding bit
        @test_throws ArgumentError decode_pubkey(mutated)

        # -0 is not a valid coefficient encoding
        @test decompress_sig(UInt8[0b10000000, 0b10000000], 2, 1) === nothing
        # ... while +0 is
        @test decompress_sig(compress_sig([0], 2), 2, 1) == [0]
    end

    @testset "round-tripping our own keys" begin
        # A key we generated ourselves must survive the format too, which is a
        # different check from decoding the C reference's: it exercises the
        # encoder's range validation on real Gaussian coefficients.
        rng = chacha20(collect(UInt8, 1:56))
        src = bytesource(rng)
        f, g, F, G = ntru_gen(512, src)
        sk = encode_privkey(Int.(f), Int.(g), Int.(F), 9)
        @test length(sk) == FALCON_512.privkey_bytes
        lg, f2, g2, F2 = decode_privkey(sk)
        @test (lg, f2, g2, F2) == (9, Int.(f), Int.(g), Int.(F))
        @test recover_G(f2, g2, F2) == G

        fq = Int[mod(c, Q) for c in f]; gq = Int[mod(c, Q) for c in g]
        h = polydivq(gq, fq)
        pk = encode_pubkey(h, 9)
        @test length(pk) == FALCON_512.pubkey_bytes
        @test decode_pubkey(pk) == (9, h)
    end

    @testset "range checking on encode" begin
        # A coefficient too large for its field must be refused, not truncated.
        f = zeros(Int, 512); f[1] = 1 << 20
        @test_throws ArgumentError encode_privkey(f, zeros(Int, 512), zeros(Int, 512), 9)
        h = zeros(Int, 512); h[1] = Q
        @test_throws ArgumentError encode_pubkey(h, 9)
        @test_throws ArgumentError encode_signature(zeros(UInt8, 39), zeros(Int, 512), 9, 666)
    end
end
