#!/usr/bin/env julia
#
# cref_contract_recover.jl -- recover the private key from a pair of signatures
# produced by two conforming BUILDS of the FALCON reference implementation.
#
#     julia --project=falcon falcon/scripts/cref_contract_recover.jl dump_a.txt dump_b.txt
#
# The two dumps come from scripts/cref_contract_diff.c, built once with
# `clang -O2 -march=native -ffp-contract=fast` and once with the same clang at
# its default.  Both are conforming C; C99 6.5p8 licenses the first to contract
# `a*b + c` into one fused multiply-add and the FALCON specification says
# nothing about FP_CONTRACT.  Same source, same key, same message, same PRNG
# tape -- and, on about one signature in ten thousand, different bytes.
#
# This script is the attacker: it sees the two signatures, the public key and
# the message.  It runs the recovery of ePrint 2024/1709 section 5.1 --
#
#     ds1 = dz0 * g,  ds2 = dz0 * (-f),  dz0 = a + b x^{n/2},
#     dz0^{-1} = (a - b x^{n/2}) / (a^2 + b^2),
#
# searching the small (a,b) and accepting the candidate that reproduces the
# public key.  s1 is not transmitted, but s1 = c - s2*h (mod q) and ||ds1|| is
# far below q/2, so the attacker reconstructs it from public data.
#
# CONSTANT TIME: irrelevant; this is the attacker's side.

using Falcon
using Printf

const F = Falcon
const p = FALCON_512
const n = p.n
const q = p.q

function readdump(path)
    key = Dict{String,Vector{UInt8}}()
    msg = Dict{Int,Vector{UInt8}}(); sig = Dict{Int,Vector{UInt8}}()
    for line in eachline(path)
        parts = split(line)
        length(parts) == 3 || continue
        tag, j, hx = parts[1], parse(Int, parts[2]), parts[3]
        b = [parse(UInt8, hx[i:i+1], base = 16) for i in 1:2:length(hx)]
        tag == "privkey" && (key["sk"] = b)
        tag == "pubkey"  && (key["pk"] = b)
        tag == "message" && (msg[j] = b)
        tag == "signature" && (sig[j] = b)
    end
    return key, msg, sig
end

redq(v) = Int[mod(c, q) for c in v]
centre(v) = Int[(y = mod(c, q); y > q ÷ 2 ? y - q : y) for c in v]

function mul_sparse(v::Vector{Int}, a::Int, b::Int)
    m = length(v); h = m ÷ 2
    out = a .* v
    for k in 0:(m - 1)
        src = k - h
        out[k + 1] -= b * (src >= 0 ? v[src + 1] : -v[src + m + 1])
    end
    return out
end

function exact_div(v::Vector{Int}, d::Int)
    d == 0 && return nothing
    out = similar(v)
    for i in eachindex(v)
        qq, rr = divrem(v[i], d)
        rr == 0 || return nothing
        out[i] = qq
    end
    return out
end

smallkey(v) = maximum(abs, v) <= 256

"Section 5.1, from the two signature differences and the public key alone."
function recover(ds1::Vector{Int}, ds2::Vector{Int}, h::Vector{Int})
    for a in -19:19, b in -19:19
        (a == 0 && b == 0) && continue
        d = a * a + b * b
        g0 = exact_div(mul_sparse(ds1, a, b), d); g0 === nothing && continue
        f0 = exact_div(mul_sparse(ds2, a, b), d); f0 === nothing && continue
        for (gc, fc) in ((g0, -f0), (-g0, f0), (g0, f0), (-g0, -f0))
            (smallkey(gc) && smallkey(fc)) || continue
            redq(gc) == redq(F.polymulq(redq(fc), h)) || continue
            return (a, b, gc, fc)
        end
    end
    return nothing
end

function main()
    fa, fb = ARGS[1], ARGS[2]
    ka, msgs, sa = readdump(fa)
    kb, _,   sb = readdump(fb)
    @assert ka["pk"] == kb["pk"] "the two builds disagree on the public key"
    _, h = decode_pubkey(ka["pk"])
    # the true secret, used only to say whether the recovery is right
    _, tf, tg, _tF = decode_privkey(ka["sk"])
    truef = Int.(tf); trueg = Int.(tg)

    println("# cref_contract_recover.jl -- two builds of the C reference, n=", n)
    println("# ", fa, "  vs  ", fb)
    println()
    nrec = 0; ntried = 0; nforge = 0
    forge = !("--no-forge" in ARGS)
    for j in sort(collect(keys(sa)))
        haskey(sb, j) || continue
        sa[j] == sb[j] && continue
        ntried += 1
        _, salta, s2a = decode_signature(sa[j])
        _, saltb, s2b = decode_signature(sb[j])
        @assert salta == saltb "the two builds drew different salts"
        c = hash_to_point(msgs[j], salta, n; q = q)
        # s1 = c - s2*h  (mod q), centre-lifted: ||s1|| << q/2
        s1a = centre(c .- F.polymulq(redq(Int.(s2a)), h))
        s1b = centre(c .- F.polymulq(redq(Int.(s2b)), h))
        ds1 = s1a .- s1b
        ds2 = Int.(s2a) .- Int.(s2b)
        nd = count(!=(0), ds1) + count(!=(0), ds2)
        @printf("message %-6d : %4d of %d signature coefficients differ  ", j, nd, 2n)
        r = recover(ds1, ds2, h)
        if r === nothing
            println("-> no (a,b): a first-two divergence")
            continue
        end
        a, b, gc, fc = r
        nrec += 1
        # Is it a USABLE key?  The public-key relation is already imposed inside
        # recover, and the norm bound is what key generation enforces.  But the
        # only check that settles the question is to USE it: complete (f,g) to a
        # full trapdoor basis by solving the NTRU equation, sign a fresh message
        # the victim never signed, and offer it to the ORIGINAL public key.
        nrm = sum(x -> x * x, gc) + sum(x -> x * x, fc)
        bound = 1.17^2 * q
        sym = symmetry(fc, gc, truef, trueg)
        @printf("-> *** KEY RECOVERED *** (a,b)=(%2d,%2d)  ||(g,-f)||^2 = %d %s bound %d  %s\n",
                a, b, nrm, nrm <= bound ? "<=" : ">", round(Int, bound), sym)
        if forge
            ok = forge_with(fc, gc, h)
            @printf("     forged a signature on a NEW message with the recovered key: %s\n",
                    ok ? "the ORIGINAL public key ACCEPTS it" : "rejected")
            ok && (nforge += 1)
        end
    end
    println()
    @printf("# %d divergent pairs; %d yielded the private key by section 5.1", ntried, nrec)
    forge && @printf(", %d of which\n#   signed a new message that the original public key accepted", nforge)
    println(";")
    @printf("# %d did not -- their dz is not 2-sparse in z0, so the divergence was not\n", ntried - nrec)
    println("# at the last two sampler calls.  Section 6.2 measured zero interior integer")
    println("# centres in 1.02e8 draws, so those are first-two divergences: the positions")
    println("# ePrint 2024/1709 section 5 sets aside, and the ones section 6.5 turns into")
    println("# one F_q equation each.")
end

"""
Complete the recovered `(f,g)` to a trapdoor basis, sign a message the victim
never signed, and verify it against the victim's public key.

This is the only check that settles "is it a working key".  The public-key
relation says the pair is *an* NTRU pair for `h`; the norm bound says it is
short enough; neither says the sampler will run on it.  Solving `fG - gF = q`
is the same Pornin-Prest descent key generation uses, and it needs nothing the
attacker does not already have.
"""
function forge_with(fc::Vector{Int}, gc::Vector{Int}, h::Vector{Int})
    try
        bF, bG = F.ntru_solve(BigInt.(fc), BigInt.(gc); q = q)
        sk = expand_privkey(fc, gc, Int.(bF), Int.(bG), p)
        pk = FalconPublicKey(p, h)
        msg = collect(codeunits("a message the victim never signed"))
        rb = chacha20(shake256(codeunits("forge"), 56))
        sig = falcon_sign(sk, msg, x -> randombytes!(rb, x))
        return falcon_verify(pk, msg, sig)
    catch e
        @printf("     (completing the basis failed: %s)\n", sprint(showerror, e))
        return false
    end
end

"""Which of the NTRU symmetries +-x^k (f,g) the recovered pair is, or "" if none."""
function symmetry(fc, gc, truef, trueg)
    m = length(fc)
    shift(v, k) = Int[(src = i - k; src >= 0 ? v[src + 1] : -v[src + m + 1])
                      for i in 0:(m - 1)]
    for k in 0:(m - 1), sgn in (1, -1)
        shiftf = shift(truef, k); shiftg = shift(trueg, k)
        if fc == sgn .* shiftf && gc == sgn .* shiftg
            return sgn > 0 ? "= x^$k (f,g)" : "= -x^$k (f,g)"
        end
    end
    return "(an equivalent key, not a symmetry of the stored one)"
end

main()
