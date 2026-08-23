#!/usr/bin/env julia
#
# key_recovery.jl -- recover the FALCON private key from a single pair of
# signatures that diverge at the last two sampler calls, exactly as ePrint
# 2024/1709 section 5.1 describes, using the A2 divergence this project found
# (complex division + D11, which that paper does not consider).
#
#     julia --project=falcon falcon/scripts/key_recovery.jl [key] [sig]
#
# This turns the paper's phrase "at the key-recovery position" into "we
# recovered the key".  The divergent pair is the A2 case at n=512 located by
# scripts/first_divergence.jl (default key 70, signature 534).
#
# ===========================================================================
# THE MATH (ePrint 2024/1709 section 5.1)
# ===========================================================================
#
# Two signatures s, s' on the SAME syndrome u differ only when ffSampling's
# last two integer samples differ, i.e. in z(2n-2) and z(2n-1), which are the
# degree-0 and degree-n/2 coefficients of z0.  Writing
#
#     dz0 = a + b*x^{n/2},   (a, b) = (z(2n-2)-z'(2n-2), z(2n-1)-z'(2n-1))
#
# the signature difference is, over R = Z[x]/(x^n+1),
#
#     ds0 = dz0 * g,   ds1 = dz0 * (-f)
#
# (up to this implementation's sign convention for the two components, which we
# resolve by trying both).  dz0 is sparse, and
#
#     (a + b*x^{n/2})(a - b*x^{n/2}) = a^2 - b^2 x^n = a^2 + b^2   (x^n = -1),
#
# so dz0^{-1} = (a - b*x^{n/2})/(a^2 + b^2): no general ring inversion is
# needed.  The attacker does not know (a, b), but a,b in {-18,...,19} because
# SamplerZ's outputs are bounded, so a search over 40*40 pairs recovers the
# key: for each (a,b), g_cand = ds0 * (a - b x^{n/2}) / (a^2+b^2); the right
# pair is the one giving integer, short coefficients that reproduce the public
# key.  The attacker sees only (s, s') and h -- not (a,b), not z.
#
# CONSTANT TIME: irrelevant here; this is the attacker's side.

using Falcon
using Printf

const F = Falcon

keyidx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 70
sigidx = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 534

p = FALCON_512
n = p.n
q = p.q

# --- reproduce the A2 divergent pair, exactly as first_divergence.jl does ---
function respelled_tree(sk)
    t = with_spec_spelling(cdiv = true, ldl = true) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

r = chacha20(collect(UInt8, 0x00:0x37))
local sk
for _ in 1:keyidx
    global sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
end
sk2 = respelled_tree(sk)
pt = hash_to_point(collect(codeunits("power analysis message")),
                   shake256(codeunits("power salt"), SALT_LEN), 512; q = q)
st = shake256(codeunits("A2/$keyidx/$sigidx"), 56)
r1 = chacha20(st); r2 = chacha20(st)
sA0, sA1 = F.sample_preimage(sk,  pt, x -> randombytes!(r1, x))   # C spelling
sB0, sB1 = F.sample_preimage(sk2, pt, x -> randombytes!(r2, x))   # spec spelling

ds0 = Int.(sA0) .- Int.(sB0)
ds1 = Int.(sA1) .- Int.(sB1)
ndiff = count(!=(0), ds0) + count(!=(0), ds1)
@printf("divergent pair: key %d sig %d, %d of %d signature coefficients differ\n",
        keyidx, sigidx, count(!=(0), ds0) + count(!=(0), ds1), 2n)
if ndiff == 0
    println("no divergence on this pair -- pick one from scripts/divergence_rate.jl")
    exit(1)
end

# --- the true secret, for verification only (the attacker does not have it) --
truef = Int.(sk.f); trueg = Int.(sk.g)

# multiply an integer vector by (a - b*x^{n/2}) in Z[x]/(x^n+1)
function mul_sparse(v::Vector{Int}, a::Int, b::Int)
    n = length(v)
    h = n ÷ 2
    out = a .* v
    # (-b x^{n/2}) * v : (x^{n/2} v)_k = v_{k-h} for k>=h, -v_{k+h} for k<h
    for k in 0:(n - 1)
        src = k - h
        val = src >= 0 ? v[src + 1] : -v[src + n + 1]
        out[k + 1] -= b * val
    end
    return out
end

# exact integer division of every coefficient by d, or `nothing`
function exact_div(v::Vector{Int}, d::Int)
    d == 0 && return nothing
    out = similar(v)
    for i in eachindex(v)
        q_, rem = divrem(v[i], d)
        rem == 0 || return nothing
        out[i] = q_
    end
    return out
end

# a plausible FALCON key coefficient is small
smallkey(v) = maximum(abs, v) <= 256

"reduce a poly mod q to canonical residues"
redq(v) = Int[mod(c, q) for c in v]

# recover: search (a,b), try both sign conventions and both components
function recover()
    for a in -19:19, b in -19:19
        (a == 0 && b == 0) && continue
        d = a * a + b * b
        # g from ds0, f from ds1 (with sign options); dz0^{-1} = (a - b x^{h})/d
        g0 = exact_div(mul_sparse(ds0, a, b), d)
        f0 = exact_div(mul_sparse(ds1, a, b), d)
        g0 === nothing && continue
        f0 === nothing && continue
        for (gc, fc) in ((g0, -f0), (-g0, f0), (g0, f0), (-g0, -f0))
            (smallkey(gc) && smallkey(fc)) || continue
            # public-key check: g ≡ f*h (mod q), independent of the true secret
            h = F.polydivq(redq(trueg), redq(truef))   # the real public key
            lhs = redq(gc)
            rhs = redq(F.polymulq(redq(fc), h))
            lhs == rhs || continue
            return (a, b, gc, fc)
        end
    end
    return nothing
end

res = recover()
if res === nothing
    println("recovery FAILED (unexpected for a genuine last-two-call divergence)")
    exit(1)
end
a, b, gc, fc = res

# (f, g) and (-f, -g) are the same key: g/f = (-g)/(-f), and both are the
# NTRU lattice's shortest vectors.  Accept either as full recovery.
exact = (gc == trueg && fc == truef)
negated = (gc == .-trueg && fc == .-truef)

@printf("\nrecovered from (s, s') and the public key h alone:\n")
@printf("  searched (a,b) in [-19,19]^2 (< 2^11 pairs, per 2024/1709 sec 5.1)\n")
@printf("  hit at (a, b) = (%d, %d)\n", a, b)
@printf("  candidate reproduces the public key h (mod q): yes  [checked in recover()]\n")
@printf("  |g_cand| max = %d, |f_cand| max = %d   (true |g| = %d, |f| = %d)\n",
        maximum(abs, gc), maximum(abs, fc), maximum(abs, trueg), maximum(abs, truef))
if exact
    println("  matches the stored secret (f, g) exactly.")
elseif negated
    println("  matches the stored secret up to the key's negation symmetry: (f,g) = -(f_true,g_true).")
    println("  This is the same signing key -- g/f = (-g)/(-f).")
else
    println("  a short (f, g) reproducing h -- the NTRU lattice's shortest vector, i.e. the key.")
end
if exact || negated
    println("\n*** FULL KEY RECOVERY from a single A2 discrepant pair. ***")
    println("    Complex division + D11 (which ePrint 2024/1709 does not consider),")
    println("    at the last two sampler calls, is sufficient for full key recovery.")
end
