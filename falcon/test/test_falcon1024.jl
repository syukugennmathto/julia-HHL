# test_falcon1024.jl -- the *other* parameter set, end to end.
#
# Everything else in the suite runs at n = 512 (plus toy degrees).  FALCON-1024
# has had its parameters in params.jl since module 1 and every algorithm is
# written for a general power of two, so "it should just work" -- which is
# exactly the kind of claim that turns out to be false at one specific place.
# The places it could plausibly break are all *table* sizes and *width* limits
# rather than algorithms:
#
#   - `ntt.jl` needs a primitive 2n-th root of unity mod q.  At n = 1024 that
#     means order 2048, and q - 1 = 12288 = 6 * 2048, so it exists -- but only
#     just.  There is no FALCON-2048, and this is why.
#   - `encoding.jl`'s `MAX_FG_BITS[logn + 1]` is 5 at logn = 10 against 8 at
#     logn = 9: the private key's f and g are packed into *fewer* bits at the
#     larger degree, because sigma_fg shrinks.  An off-by-one in that table
#     shows up here and nowhere else.
#   - `ffsampling.jl`'s tree is one level deeper.
#   - `samplerz.jl`'s sigma_min is different (params.jl), so the rejection
#     rate differs.
#
# COST: key generation at n = 1024 is slow in this implementation -- `ntru_solve`
# is the bottleneck and it carries BigInt (docs/debug_log.md #032).  One key is
# generated here, once, and reused for every check below.  If this file starts
# dominating the suite's runtime, the fix is to make `ntru_solve` faster, not to
# delete the test.

@testset "falcon-1024" begin

    rng = chacha20(collect(UInt8, 0x40:0x77))
    src = bytesource(rng)
    p = FALCON_1024

    sk, pk = falcon_keygen(1024, src)

    @testset "the key is internally consistent" begin
        @test sk.params === FALCON_1024
        @test pk.params === FALCON_1024
        @test length(sk.f) == length(sk.g) == length(sk.F) == length(sk.G) == 1024
        # f*G - g*F = q, exactly, over the integers.  No oracle needed.
        @test ntru_equation_holds(BigInt.(sk.f), BigInt.(sk.g),
                                  BigInt.(sk.F), BigInt.(sk.G))
        # h = g/f mod q, and h*f = g mod q as a consequence
        @test public_key(sk).h == pk.h
        @test polymulq(pk.h, Int[mod(c, Q) for c in sk.f]) ==
              Int[mod(c, Q) for c in sk.g]
        # the basis is short enough for the scheme's own criterion
        @test gs_norm_ok(sk.f, sk.g; q = p.q)
    end

    @testset "sign and verify" begin
        for msg in (b"", b"a", b"hello falcon-1024", collect(UInt8, 0:255))
            sig = falcon_sign(sk, msg, src)
            @test length(sig) == p.sig_bytes
            @test falcon_verify(pk, msg, sig)
            @test signature_norm(pk, msg, sig) <= p.sig_bound
            # ... and the same signature must NOT verify at the other degree,
            # because the encoded header carries logn
            @test decode_signature(sig)[1] == 10
        end
    end

    @testset "a 512 verifier rejects a 1024 signature and vice versa" begin
        # The header byte is 0x30 + logn, so this is really testing that
        # `falcon_verify` checks it rather than trusting the length.
        sk5, pk5 = falcon_keygen(512, bytesource(chacha20(collect(UInt8, 0x80:0xb7))))
        sig10 = falcon_sign(sk, b"cross", src)
        sig5 = falcon_sign(sk5, b"cross", src)
        @test !falcon_verify(pk5, b"cross", sig10)
        @test !falcon_verify(pk, b"cross", sig5)
    end

    @testset "keys round-trip through their serialised form" begin
        # This is where MAX_FG_BITS[11] = 5 gets exercised: at logn = 10 the
        # private key packs f and g into 5 bits each, not 8.
        skb = privkey_bytes(sk)
        pkb = pubkey_bytes(pk)
        @test length(skb) == p.privkey_bytes
        @test length(pkb) == p.pubkey_bytes
        @test skb[1] == 0x50 + 0x0a            # header: 0101 nnnn, logn = 10
        @test pkb[1] == 0x00 + 0x0a

        sk2 = privkey_from_bytes(skb)
        pk2 = pubkey_from_bytes(pkb)
        @test sk2.f == sk.f && sk2.g == sk.g && sk2.F == sk.F && sk2.G == sk.G
        @test pk2.h == pk.h

        # and the reconstructed key still signs into the original public key
        sig = falcon_sign(sk2, b"round trip", src)
        @test falcon_verify(pk, b"round trip", sig)
    end

    @testset "signature norms sit where theory says" begin
        # Same check as at n = 512, with n doubled: E||(s1,s2)||^2 = 2n*sigma^2.
        # sigma is *larger* at 1024 (params.jl), so this is not just the 512
        # number times two, and asserting it separately is the point.
        norms = [signature_norm(pk, "message $i", falcon_sign(sk, "message $i", src))
                 for i in 1:8]
        expected = 2 * 1024 * p.sigma^2
        @test 0.8 * expected < sum(norms) / length(norms) < 1.2 * expected
        @test all(nn -> nn <= p.sig_bound, norms)
    end
end
