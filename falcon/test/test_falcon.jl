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

    @testset "hash_to_point is unchanged by the block squeeze" begin
        # Regression test for docs/debug_log.md #038.  hash_to_point used to
        # draw two bytes at a time; it now squeezes a block and indexes into it,
        # and reduces by conditional subtraction rather than by `%`.  Both are
        # meant to be invisible.  The oracle is the original formulation, kept
        # here verbatim, because "same bytes read differently" is exactly the
        # kind of claim that is easy to get subtly wrong at a block boundary.
        function h2p_twobytes(message, salt, n; q = Q)
            k = (1 << 16) ÷ q
            bound = k * q
            xof = shake256_xof(salt, message)
            out = Int[]
            while length(out) < n
                b = squeeze!(xof, 2)
                elt = (Int(b[1]) << 8) + Int(b[2])
                elt < bound && push!(out, elt % q)
            end
            return out
        end
        rng = MersenneTwister(20260822)
        for n in (2, 8, 512, 1024), mlen in (0, 1, 17, 200)
            m = rand(rng, UInt8, mlen)
            s = rand(rng, UInt8, SALT_LEN)
            @test hash_to_point(m, s, n) == h2p_twobytes(m, s, n)
        end
        # every coefficient is a valid residue, and the rejection really rejects
        pt = hash_to_point(rand(rng, UInt8, 32), rand(rng, UInt8, SALT_LEN), 1024)
        @test all(c -> 0 <= c < Q, pt)
        @test length(unique(pt)) > 500          # not degenerate
    end

    @testset "signing reproduces the C reference byte for byte" begin
        # THIS BRANCH's headline result (docs/debug_log.md #050).
        #
        # Since #025 this project could verify the C reference's signatures but
        # not *produce* them: given the same key, message and randomness, the
        # two implementations returned different (valid) signatures.  #048
        # bisected the floating-point path and found three places where the
        # specification's formulas and the reference's spelling are
        # algebraically identical and differ in rounding.  With those three
        # respelled -- complex division, LDL*'s D11, and the bottom two levels
        # of ffSampling -- plus the reference's reciprocal leaf convention and
        # its `fpr_inv_sigma` table, signing agrees exactly.
        #
        # The vectors are the C reference's own output, recorded by
        # scripts/gen_cref_sign_kat.jl through a shim that seeds the sampler's
        # ChaCha20 state directly (scripts/cref_shim.c), so the byte stream is
        # pinned as well as the key and the message.  Regenerating them needs a
        # C compiler; checking them does not.
        p = FALCON_512
        @test length(CREF_SIGN_KAT) == 8
        for (f, g, F, G, salt, state, pt, want_s2, want_bytes) in CREF_SIGN_KAT
            sk = expand_privkey(f, g, F, G, p)

            # 1. the expanded key: same basis, same tree
            @test length(leaf_sigmas(sk.tree)) == p.n

            # 2. the sampled short vector, from the same 56-byte PRNG state
            rng = chacha20(state)
            s1, s2 = Falcon.sample_preimage(sk, pt, k -> randombytes!(rng, k))
            @test Int.(s2) == want_s2

            # 3. the compressed bytes.  C's encoder returns the natural length;
            #    the padded format zero-fills to sig_bytes - 41, which is what
            #    `compress_sig` produces.
            ours = compress_sig(Int.(s2), p.sig_bytes - 41)
            @test ours !== nothing
            @test ours[1:length(want_bytes)] == want_bytes
            @test all(iszero, ours[(length(want_bytes) + 1):end])

            # 4. and it is a signature that verifies, which the equality above
            #    does not by itself establish.
            h = polydivq(Int[mod(c, p.q) for c in g], Int[mod(c, p.q) for c in f])
            pk = FalconPublicKey(p, h)
            sig = encode_signature(salt, Int.(s2), p.logn, p.sig_bytes)
            @test sig !== nothing
            @test sqnorm(s1, s2) <= p.sig_bound
        end
    end

    @testset "the two spellings differ in bits, and agree on these 8 vectors" begin
        # The finding whose interpretation this project got wrong three times
        # and then corrected (docs/debug_log.md #048, #050, #053, #054, #056).
        #
        # #048 identified three places where the specification's formulas and
        # the C reference's spelling are algebraically identical and round
        # differently.  The tempting conclusion -- that these differences never
        # reach the signature -- was asserted here on the strength of 480
        # signatures showing zero divergence.  That was under-powered: the
        # divergence rate is on the order of 1e-5 per signature (A1, the
        # hand-unrolled bottom levels, at 5e-5; A2, complex division and D11,
        # at 7.5e-6 -- scripts/divergence_rate.jl), and 480 signatures cannot
        # distinguish 1e-5 from 0.  The mechanism is NOT a sub-ulp BerExp
        # margin as once claimed here; it is `s = floor(mu)` straddling an
        # integer centre (ePrint 2024/1709, Lemma 1), which is a difference of
        # exactly 1, not a rare coin flip.
        #
        # What is true, and all this test now asserts, is narrower: the
        # arithmetic genuinely differs bit-for-bit, and on these 8 specific KAT
        # vectors the two spellings happen to agree (none is one of the rare
        # divergent cases).  The divergence itself is measured, at scale, in
        # scripts/divergence_rate.jl, not here.
        p = FALCON_512

        # (a) the arithmetic really is different -- otherwise (b) is vacuous
        let a = ComplexF64(1.0, 3.0), b = ComplexF64(7.0, 11.0)
            @test Falcon._cdiv_cref(a, b) != a / b
        end

        # (b) on these 8 vectors specifically, both spellings give C's output.
        #     This is NOT a claim that they always agree -- see #056.
        for (f, g, F, G, salt, state, pt, want_s2, _) in CREF_SIGN_KAT
            sk = expand_privkey(f, g, F, G, p)
            r1 = chacha20(state); r2 = chacha20(state)
            _, a = Falcon.sample_preimage(sk, pt, k -> randombytes!(r1, k))
            _, b = with_spec_ffsampling() do
                Falcon.sample_preimage(sk, pt, k -> randombytes!(r2, k))
            end
            @test Int.(a) == Int.(b)
            @test Int.(a) == want_s2
        end
    end

    @testset "blind which-call solve: rank-2 lift recovers f without labels" begin
        # docs/debug_log.md #075, docs/paper.md sec 6.6.  The first-two event
        # channel with UNKNOWN which-call label: each event gives one of two rows
        # a=form_row(c,k1), b=form_row(c,k2), one satisfying row.f=0.  Since
        # b = x^{n/2} a, the disjunction is (a.f)(a.f')=0 with f'=-x^{n/2}f, a
        # linear measurement on the rank-2 S=sym(f f'^T).  Checked at n=16 with a
        # fixed seed: S is pinned (kernel dim 1) and its column space contains f.
        qq = 12289
        nn = 16; k1 = nn ÷ 2 - 1; k2 = nn - 1; Dd = nn*(nn+1) ÷ 2
        rng = MersenneTwister(20260823)
        fsec = Int[rand(rng, -5:5) for _ in 1:nn]; fsec[1] = fsec[1] == 0 ? 1 : fsec[1]
        shp(v, m) = Int[(i = j - m; i >= 0 ? v[i+1] : -v[i+nn+1]) for j in 0:nn-1]
        fp = mod.(-shp(fsec, nn ÷ 2), qq)
        frw(c, k) = Int[(i = k - j; mod(i >= 0 ? c[i+1] : -c[i+nn+1], qq)) for j in 0:nn-1]
        sm(a) = (o = Int[]; for i in 1:nn, j in i:nn; push!(o, i==j ? mod(a[i]*a[j], qq) : mod(2*a[i]*a[j], qq)); end; o)
        function rrefq(M0)
            M = [mod(x, qq) for x in M0]; rows, cols = size(M); piv = Int[]; r = 1
            for col in 1:cols
                pr = findfirst(i -> M[i,col] % qq != 0, r:rows); pr === nothing && continue; pr += r - 1
                M[r,:], M[pr,:] = M[pr,:], M[r,:]; iv = invmod(M[r,col], qq); M[r,:] = mod.(M[r,:] .* iv, qq)
                for i in 1:rows; i == r && continue; fc = M[i,col]; fc == 0 && continue; M[i,:] = mod.(M[i,:] .- fc .* M[r,:], qq); end
                push!(piv, col); r += 1; r > rows && break
            end
            M, piv
        end
        meas = Vector{Int}[]
        for _ in 1:Int(round(1.6*Dd))
            k = rand(rng, (k1, k2)); c = Int[rand(rng, 0:qq-1) for _ in 1:nn]
            row = frw(c, k); sv = mod(sum(row[j]*fsec[j] for j in 1:nn), qq)
            j0 = findfirst(j -> gcd(fsec[j], qq) == 1, 1:nn); p0 = k - (j0-1)
            if p0 >= 0; c[p0+1] = mod(c[p0+1] - sv*invmod(fsec[j0], qq), qq)
            else; c[p0+nn+1] = mod(c[p0+nn+1] + sv*invmod(fsec[j0], qq), qq) end
            push!(meas, sm(frw(c, k1)))          # solver only ever forms a_i (row at k1)
        end
        A = reduce(vcat, [reshape(r, 1, :) for r in meas])
        R, piv = rrefq(A); free = setdiff(1:Dd, piv)
        @test length(free) == 1                  # S pinned up to scale (kernel dim 1)
        Svec = zeros(Int, Dd); Svec[free[1]] = 1
        for (ri, cc) in enumerate(piv); Svec[cc] = mod(-R[ri, free[1]], qq); end
        S = zeros(Int, nn, nn); t = 1
        for i in 1:nn, j in i:nn; S[i,j] = Svec[t]; S[j,i] = Svec[t]; t += 1; end
        # colspace(S) contains f  <=>  rank[S | f] == rank[S]
        _, ps = rrefq(S); _, pa = rrefq(hcat(S, reshape(mod.(fsec, qq), nn, 1)))
        @test length(ps) == 2                     # rank-2
        @test length(pa) == length(ps)            # f is in the column space
    end

    @testset "the FMA arm, and the parity that closes the first two calls" begin
        # docs/debug_log.md #070, docs/paper.md sec 6.6.
        #
        # (a) FMA_FFT is OFF by default.  Every byte-exact KAT above depends on
        #     that, so this is the guard against a future edit flipping it.
        @test Falcon.FMA_FFT[] == false

        # (b) with_fma restores the flag even when the body throws.
        @test_throws ErrorException with_fma(() -> error("boom"))
        @test Falcon.FMA_FFT[] == false

        # (c) the arm is not a no-op, and -- the reason it is `fma` and not
        #     `muladd` -- it is not a no-op in EVERY build configuration.  The
        #     first version of this arm used `muladd`, which is only permitted
        #     to fuse: LLVM contracted it inside merge_fft's loop, declined
        #     inside mul_fft's comprehension, and declined in both under
        #     --check-bounds=yes, so this very assertion passed under `julia
        #     script.jl` and failed under `Pkg.test()`.  That is a faithful
        #     picture of the hazard sec 6.6 is about and a useless experimental
        #     arm.  `fma` rounds once by specification.
        let a = ComplexF64(0.1, 1.0), b = ComplexF64(0.2, 0.020000000000000004)
            @test Falcon._cmul_fma(a, b) != a * b
        end
        n = 512
        v = Float64[Float64(mod(7*i*i + 3i, 12289)) for i in 1:n]
        w = Float64[Float64(mod(5*i*i + 11i, 12289)) for i in 1:n]
        @test polymul_fft(v, w) != with_fma(() -> polymul_fft(v, w))

        # (d) the parity theorem of sec 6.6.  With `round` in place of `floor`
        #     (ePrint 2024/1709 Algorithm 4) the sensitive centres are the
        #     half-integers, and a first-two centre is N/q with N an integer.
        #     N/q = m + 1/2 needs 2N = q(2m+1): even on the left, odd on the
        #     right.  So for odd q it never happens -- part 1 of that
        #     countermeasure closes the first two sampler calls unconditionally.
        #     Checked exhaustively over the numerator's residues mod 2q.
        q = FALCON_512.q
        @test isodd(q)
        @test !any(N -> mod(2N, 2q) == mod(q, 2q), 0:(2q - 1))
    end

    @testset "key recovery algebra (ePrint 2024/1709 sec 5.1)" begin
        # The identity behind scripts/key_recovery.jl, which recovers the
        # private key from a single A2 discrepant pair (docs/paper.md sec 6.1).
        # A last-two-call divergence gives, over R = Z[x]/(x^n+1),
        #     Δs0 = Δz0 · g,   Δz0 = a + b x^{n/2},
        # and (a + b x^{n/2})(a - b x^{n/2}) = a^2 + b^2 because x^n = -1, so
        #     g = Δs0 · (a - b x^{n/2}) / (a^2 + b^2)
        # with no general ring inversion.  This tests that recovery on synthetic
        # data, so the full 382-pair search in the script is regression-guarded
        # without paying for the 70-key reproduction.
        ringmul(u, v) = begin              # multiply in Z[x]/(x^n+1)
            n = length(u); w = zeros(Int, n)
            for i in 0:n-1, j in 0:n-1
                k = i + j; c = u[i+1] * v[j+1]
                w[mod(k, n) + 1] += k < n ? c : -c
            end
            w
        end
        mulsparse(v, a, b) = begin          # v * (a - b x^{n/2}) in the ring
            n = length(v); h = n ÷ 2; out = a .* v
            for k in 0:n-1
                src = k - h
                out[k+1] -= b * (src >= 0 ? v[src+1] : -v[src+n+1])
            end
            out
        end
        rng = MersenneTwister(0xFA1C0)
        for n in (8, 16, 32), _ in 1:20
            g = rand(rng, -12:12, n); f = rand(rng, -12:12, n)
            a = rand(rng, -19:19); b = rand(rng, -19:19)
            (a == 0 && b == 0) && continue
            dz0 = zeros(Int, n); dz0[1] = a; dz0[n ÷ 2 + 1] = b
            ds0 = ringmul(dz0, g); ds1 = ringmul(dz0, .-f)
            d = a * a + b * b
            grec = mulsparse(ds0, a, b) .÷ d
            frec = .-(mulsparse(ds1, a, b) .÷ d)
            @test all(iszero, mulsparse(ds0, a, b) .% d)   # exact division
            @test grec == g
            @test frec == f
        end
    end
end
