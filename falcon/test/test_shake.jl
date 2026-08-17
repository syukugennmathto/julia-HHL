# test_shake.jl -- module 2.
#
# Oracles:
#   * SHAKE256: CPython's hashlib (an independent FIPS 202 implementation),
#     via the vectors in test/vectors/shake256_kat.jl.  Our Keccak is our own,
#     so these vectors are the only thing standing between a mistyped rotation
#     offset and a hash that is wrong but looks perfectly random.
#   * ChaCha20 PRNG: the Python reference implementation, via
#     test/vectors/chacha20_kat.jl.

# Golden vectors are included by runtests.jl.

@testset "shake" begin

    @testset "SHAKE256 one-shot against hashlib" begin
        for (msg, dlen, expected) in SHAKE256_KAT
            @test shake256(msg, dlen) == expected
        end
    end

    @testset "SHAKE256 prefix property" begin
        # SHAKE256(m, d1) must be a prefix of SHAKE256(m, d2) for d1 < d2.
        # This started life as a check on SHA.jl, whose shake256 re-pads on
        # every digest! call.  That dependency is gone -- SHA v0.7.0, the
        # version Julia actually bundles, has no SHAKE at all (debug_log #014)
        # -- so we now own the sponge and this tests our own squeeze.
        for (msg, d1, d2, expected_long) in SHAKE256_PREFIX
            @test shake256(msg, d2) == expected_long
            @test shake256(msg, d1) == expected_long[1:d1]
        end
    end

    @testset "XOF incremental squeeze matches one-shot" begin
        for (msg, _, _) in SHAKE256_KAT
            total = 300
            reference = shake256(msg, total)

            # one byte at a time
            x = shake256_xof(msg)
            got = UInt8[]
            for _ in 1:total
                append!(got, squeeze!(x, 1))
            end
            @test got == reference

            # irregular chunk sizes, including a zero-length squeeze
            x = shake256_xof(msg)
            got = UInt8[]
            for k in (0, 1, 7, 0, 64, 100, 128)
                append!(got, squeeze!(x, k))
            end
            @test got == reference[1:length(got)]
            @test length(got) == 300
        end
    end

    @testset "XOF absorbs a concatenation" begin
        a = UInt8[0x01, 0x02, 0x03]
        b = UInt8[0xfe, 0xff]
        @test squeeze!(shake256_xof(a, b), 64) == shake256(vcat(a, b), 64)
        @test squeeze!(shake256_xof("abc"), 32) == shake256(b"abc", 32)
        # Concatenation is *not* domain separation: these must collide, and
        # knowing that they do is why the salt has a fixed length.
        @test squeeze!(shake256_xof(UInt8[0x01], UInt8[0x02]), 16) ==
              squeeze!(shake256_xof(UInt8[0x01, 0x02]), 16)
    end

    @testset "XOF argument checking" begin
        x = shake256_xof(UInt8[])
        @test_throws ArgumentError squeeze!(x, -1)
        # absorbing after squeezing is an error: the padding is already applied
        y = shake256_xof(UInt8[0x01])
        squeeze!(y, 1)
        @test_throws ArgumentError absorb!(y, UInt8[0x02])
    end

    @testset "absorbing in pieces equals absorbing at once" begin
        msg = collect(UInt8, 0:250)
        for cut in (0, 1, 135, 136, 137, 200, 251)
            x = SHAKE256XOF()
            absorb!(x, msg[1:cut])
            absorb!(x, msg[(cut + 1):end])
            @test squeeze!(x, 64) == shake256(msg, 64)
        end
    end

    @testset "ChaCha20 PRNG against the reference" begin
        for (seed, sizes, expected) in CHACHA20_KAT
            rng = chacha20(seed)
            for (k, want) in zip(sizes, expected)
                @test randombytes!(rng, k) == want
            end
        end
    end

    @testset "ChaCha20 buffer boundary behaviour" begin
        # The reference *discards* the tail of the buffer when a request does
        # not fit.  So consuming 511 bytes and then 8 must NOT return bytes
        # 512..519 of the concatenated stream: it returns the first 8 bytes of
        # a freshly generated buffer.  This is the single most surprising
        # property of the reference PRNG, so it gets its own test rather than
        # relying on the KAT to notice.
        seed = collect(UInt8, 0:55)

        a = chacha20(seed)
        first511 = randombytes!(a, 511)
        after = randombytes!(a, 8)          # forces a refill; byte 512 is discarded

        # b consumes the first buffer exactly (256 + 256 = 512), so its next
        # request is served from the head of the *second* buffer.
        b = chacha20(seed)
        buffer1 = vcat(randombytes!(b, 256), randombytes!(b, 256))
        head_of_buffer2 = randombytes!(b, 8)

        @test length(buffer1) == 512
        @test buffer1[1:511] == first511
        @test after == head_of_buffer2      # a skipped straight to the new buffer
        @test after != buffer1[505:512]     # ... rather than continuing the old one
    end

    @testset "ChaCha20 argument checking" begin
        @test_throws ArgumentError chacha20(UInt8[])
        @test_throws ArgumentError chacha20(zeros(UInt8, 55))
        @test_throws ArgumentError chacha20(zeros(UInt8, 57))
        rng = chacha20(zeros(UInt8, SEED_LEN))
        @test_throws ArgumentError randombytes!(rng, -1)
        @test_throws ArgumentError randombytes!(rng, 513)
        @test randombytes!(rng, 0) == UInt8[]
    end

    @testset "ChaCha20 is deterministic and seed-sensitive" begin
        s1 = collect(UInt8, 0:55)
        s2 = copy(s1); s2[end] ⊻= 0x01
        @test randombytes!(chacha20(s1), 64) == randombytes!(chacha20(s1), 64)
        @test randombytes!(chacha20(s1), 64) != randombytes!(chacha20(s2), 64)
        # A one-bit change in the seed word that feeds the counter must also
        # change the stream (it XORs into the state, so it had better).
        s3 = copy(s1); s3[41] ⊻= 0x01     # byte 40 = word s[10], zero-indexed
        @test randombytes!(chacha20(s1), 64) != randombytes!(chacha20(s3), 64)
    end
end
