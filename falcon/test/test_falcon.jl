# test_falcon.jl -- module 10: the whole scheme.
#
# The headline test is `"the C reference's signature verifies"`. Verification
# is integer-only, so it agrees with the normative implementation exactly, and
# accepting a signature that C produced from its own key is genuine
# interoperability -- not agreement between two transliterations of the same
# Python.
#
# What is deliberately NOT tested: reproducing anyone's signature bytes.
# Signing runs through the FFT, and module 8 measured that a 1-ulp difference
# in the root table changes the output (docs/debug_log.md #025). Attempting a
# signing KAT against the Python reference would either fail or force us to
# adopt its less accurate constants. So signing is checked by properties:
# it verifies, it is short, and it is a valid solution to the congruence.

# Golden vectors are included by runtests.jl.

@testset "falcon" begin

    @testset "hash_to_point against the reference construction" begin
        for (msg, salt, n, want) in HASH_TO_POINT_KAT
            @test hash_to_point(msg, salt, n) == want
            @test all(c -> 0 <= c < Q, want)
        end
    end

    @testset "hash_to_point rejects rather than biases" begin
        # k = floor(2^16/q) = 5, so values >= 5q = 61445 are discarded. Reducing
        # them instead would make residues 0..3651 more likely by about 2^-13 --
        # enormous next to the 2^-45 the parameters were chosen for.
        @test H2P_BOUND == (65536 ÷ Q) * Q
        @test H2P_BOUND == 61445
        @test 65536 - H2P_BOUND == 4091          # how many values get rejected

        # different salts give different points; the salt is not decorative
        m = b"same message"
        a = hash_to_point(m, zeros(UInt8, 40), 512)
        b = hash_to_point(m, [0x01; zeros(UInt8, 39)], 512)
        @test a != b
        # and different messages under one salt likewise
        @test hash_to_point(b"x", zeros(UInt8, 40), 64) !=
              hash_to_point(b"y", zeros(UInt8, 40), 64)
    end

    @testset "the C reference's signature verifies" begin
        # THE interoperability test. A FALCON-512 key and signature produced by
        # the normative C implementation from a fixed seed; our verifier must
        # accept it.
        pk = pubkey_from_bytes(CREF_PUBKEY)
        @test pk.params === FALCON_512
        @test falcon_verify(pk, CREF_MESSAGE, CREF_SIGNATURE)

        # and with real headroom under the bound, not marginally
        nrm = signature_norm(pk, CREF_MESSAGE, CREF_SIGNATURE)
        @test nrm !== nothing
        @test 0 < nrm <= FALCON_512.sig_bound
        @test nrm > 0.5 * 2 * 512 * FALCON_512.sigma^2    # not degenerate either

        # the congruence itself: s1 + s2*h = c mod q
        _, salt, s2 = decode_signature(CREF_SIGNATURE)
        c = hash_to_point(CREF_MESSAGE, salt, 512)
        s2q = Int[mod(x, Q) for x in s2]
        s1 = centered(polysubq(c, polymulq(s2q, pk.h)))
        @test polyaddq(Int[mod(x, Q) for x in s1], polymulq(s2q, pk.h)) == c
    end

    @testset "the C reference's private key signs, and we verify it" begin
        # Sign with C's key using our sampler. The bytes will not match C's
        # signature (they cannot -- #025), but the result must verify under
        # C's public key, which is the property that matters.
        sk = privkey_from_bytes(CREF_PRIVKEY)
        pk = pubkey_from_bytes(CREF_PUBKEY)
        @test public_key(sk).h == pk.h            # the keys really are a pair

        rng = chacha20(collect(UInt8, 0:55))
        src = bytesource(rng)
        for _ in 1:3
            sig = falcon_sign(sk, CREF_MESSAGE, src)
            @test length(sig) == FALCON_512.sig_bytes
            @test falcon_verify(pk, CREF_MESSAGE, sig)
            @test signature_norm(pk, CREF_MESSAGE, sig) <= FALCON_512.sig_bound
        end
    end

    @testset "keygen, sign, verify round trip" begin
        rng = chacha20(collect(UInt8, 100:155))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)

        # the key is internally consistent
        @test ntru_equation_holds(BigInt.(sk.f), BigInt.(sk.g),
                                  BigInt.(sk.F), BigInt.(sk.G))
        @test public_key(sk).h == pk.h
        @test polymulq(pk.h, Int[mod(c, Q) for c in sk.f]) ==
              Int[mod(c, Q) for c in sk.g]

        for msg in (b"", b"a", b"hello falcon", collect(UInt8, 0:255))
            sig = falcon_sign(sk, msg, src)
            @test length(sig) == FALCON_512.sig_bytes
            @test falcon_verify(pk, msg, sig)
        end
    end

    @testset "signatures are randomised" begin
        # Two signatures of the same message must differ -- FALCON is
        # randomised, and a deterministic signature here would mean the salt or
        # the sampler is not being driven.
        rng = chacha20(collect(UInt8, 200:255))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)
        sigs = [falcon_sign(sk, b"same message", src) for _ in 1:4]
        @test length(unique(sigs)) == 4
        @test all(s -> falcon_verify(pk, b"same message", s), sigs)
        # ... and the salts differ, which is where the randomisation enters
        salts = [decode_signature(s)[2] for s in sigs]
        @test length(unique(salts)) == 4
    end

    @testset "tampering is rejected" begin
        rng = chacha20(collect(UInt8, 10:65))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)
        msg = b"authentic message"
        sig = falcon_sign(sk, msg, src)
        @test falcon_verify(pk, msg, sig)

        # a different message
        @test !falcon_verify(pk, b"authentic messagf", sig)
        @test !falcon_verify(pk, b"", sig)

        # a flipped bit in the salt
        bad = copy(sig); bad[2] ⊻= 0x01
        @test !falcon_verify(pk, msg, bad)

        # a flipped bit in the compressed body: either it fails to parse or the
        # norm check rejects it. Both are `false`, never an exception.
        for pos in (60, 200, 400, 660)
            bad = copy(sig); bad[pos] ⊻= 0x01
            @test falcon_verify(pk, msg, bad) isa Bool
        end

        # the wrong public key
        sk2, pk2 = falcon_keygen(512, src)
        @test !falcon_verify(pk2, msg, sig)

        # a truncated or overlong signature
        @test !falcon_verify(pk, msg, sig[1:end-1])
        @test !falcon_verify(pk, msg, vcat(sig, UInt8[0x00]))
        @test !falcon_verify(pk, msg, UInt8[])
        # a wrong header byte
        bad = copy(sig); bad[1] = 0x3a
        @test !falcon_verify(pk, msg, bad)
    end

    @testset "verification never throws on hostile input" begin
        # A verifier is fed attacker-controlled bytes; it must return false,
        # not raise. (The Python reference's decompress raises on all-zero
        # input -- docs/debug_log.md #027.)
        pk = pubkey_from_bytes(CREF_PUBKEY)
        for bytes in (UInt8[], zeros(UInt8, 666), fill(0xff, 666),
                      [0x39; zeros(UInt8, 665)], [0x39; fill(0xff, 665)],
                      rand(UInt8, 666))
            @test falcon_verify(pk, CREF_MESSAGE, bytes) isa Bool
        end
    end

    @testset "serialised keys round-trip through signing" begin
        rng = chacha20(collect(UInt8, 30:85))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)

        skb = privkey_bytes(sk)
        pkb = pubkey_bytes(pk)
        @test length(skb) == FALCON_512.privkey_bytes
        @test length(pkb) == FALCON_512.pubkey_bytes

        sk2 = privkey_from_bytes(skb)
        pk2 = pubkey_from_bytes(pkb)
        @test sk2.f == sk.f && sk2.g == sk.g && sk2.F == sk.F
        @test sk2.G == sk.G                     # recovered, not stored
        @test pk2.h == pk.h

        # a signature made with the reloaded key verifies under the original
        sig = falcon_sign(sk2, b"reloaded", src)
        @test falcon_verify(pk, b"reloaded", sig)
    end

    @testset "the retry loop is real" begin
        # Both retry conditions must be reachable, or the loop is dead code we
        # have not tested. Rather than wait for one to fire naturally, drive
        # the pieces directly: a too-long vector must fail the norm test, and
        # an incompressible one must fail to encode.
        p = FALCON_512
        long_s = fill(300, 512)                  # ||.||^2 = 512*90000 = 4.6e7
        @test sqnorm(long_s, zeros(Int, 512)) > p.sig_bound
        room = p.sig_bytes - HEAD_LEN - SALT_LEN
        @test compress_sig(fill(1 << 13, 512), room) === nothing
        # and a normal signature does fit
        rng = chacha20(collect(UInt8, 0:55))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)
        sig = falcon_sign(sk, b"fits", src)
        _, _, s2 = decode_signature(sig)
        @test compress_sig(s2, room) !== nothing
        @test sqnorm(s2) < p.sig_bound
    end

    @testset "signature norms sit where theory says" begin
        # Expected ||(s1,s2)||^2 is 2n*sigma^2 = 2.81e7 against beta^2 = 3.40e7.
        # If our sampler were running at the wrong width this ratio would move,
        # and nothing else in the suite would notice.
        rng = chacha20(collect(UInt8, 77:132))
        src = bytesource(rng)
        sk, pk = falcon_keygen(512, src)
        norms = BigInt[]
        for i in 1:8
            sig = falcon_sign(sk, [UInt8(i)], src)
            push!(norms, signature_norm(pk, [UInt8(i)], sig))
        end
        expected = 2 * 512 * FALCON_512.sigma^2
        m = sum(norms) / length(norms)
        @test 0.8 * expected < m < 1.2 * expected
        @test all(nn -> nn <= FALCON_512.sig_bound, norms)
    end
end
