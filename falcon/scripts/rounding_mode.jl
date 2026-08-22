#!/usr/bin/env julia
#
# rounding_mode.jl -- the floating-point rounding mode is process-global mutable
# state that FALCON never sets, and it changes the signature.
#
#     julia --project=falcon falcon/scripts/rounding_mode.jl [keys] [sigs]
#
# Everything measured elsewhere in this project perturbs the computation by
# rewriting SOURCE (a different but algebraically equal spelling).  This
# perturbs it without touching the source at all: the IEEE-754 rounding
# direction lives in the x87 control word and in MXCSR, it is per-thread process
# state, and nothing in the FALCON specification, the C reference, or this
# implementation ever sets it.  Any library in the address space may change it
# -- and some do.
#
# For a standard that intends bit-exact known-answer tests, this matters
# directly: the bytes of a FALCON signature are not a function of (key, message,
# randomness) alone.  They are a function of (key, message, randomness, MXCSR).
#
# Control (verified before the measurement): with the same inputs, a plain
# floating-point reduction gives 7.564912371767929 under FE_TONEAREST,
# 7.564912371768772 under FE_UPWARD and 7.564912371767021 under FE_DOWNWARD, so
# the mode is live in compiled Julia code.
#
# The key is generated once, under FE_TONEAREST, and both runs use it, so the
# expanded key and its tree are identical; only the SIGNING arithmetic differs.
# The PRNG state is identical too.  Any divergence is therefore attributable to
# the rounding mode alone.

using Falcon
using Printf

const F = Falcon

const FE_TONEAREST  = Cint(0)
const FE_DOWNWARD   = Cint(0x400)
const FE_UPWARD     = Cint(0x800)
const FE_TOWARDZERO = Cint(0xc00)
setround(m) = ccall(:fesetround, Cint, (Cint,), m)

Falcon.eval(quote
    const MU = Float64[]
    const ON = Ref(false)
    const CALLN = Ref(0)
    const FLIP_AT = Ref(0)      # 0 = never; otherwise flip the mode at this call
    const FLIP_MODE = Ref(Cint(0))
    function samplerz_isigma(mu::Float64, isigma::Float64, sigmin::Float64, rb)
        CALLN[] += 1
        if FLIP_AT[] != 0 && CALLN[] == FLIP_AT[]
            ccall(:fesetround, Cint, (Cint,), FLIP_MODE[])
        end
        ON[] && push!(MU, mu)
        s = Int(floor(mu)); r = mu - s
        dss = 0.5*(isigma*isigma); ccs = isigma*sigmin
        while true
            z0 = basesampler(rb); b = Int(rb(1)[1]) & 1; z = b + (2b-1)*z0
            x = ((z-r)^2)*dss - (z0^2)*INV_2SIGMA2
            berexp(x, ccs, rb) && return z + s
        end
    end
end)

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    p = FALCON_512; n = p.n; L = 2n
    pt = hash_to_point(collect(codeunits("power analysis message")),
                       shake256(codeunits("power salt"), SALT_LEN), n; q = p.q)

    println("# rounding_mode.jl  n=", n, "  ", nkeys, " keys x ", nsig, " signatures per mode")
    println("# key and PRNG identical; only the IEEE-754 rounding direction differs")
    println()

    setround(FE_TONEAREST)
    r = chacha20(collect(UInt8, 0x00:0x37))
    keys = [falcon_keygen(n, k -> randombytes!(r, k))[1] for _ in 1:nkeys]

    for (name, mode) in (("control(nearest)", FE_TONEAREST),
                         ("FE_UPWARD", FE_UPWARD), ("FE_DOWNWARD", FE_DOWNWARD),
                         ("FE_TOWARDZERO", FE_TOWARDZERO))
        ndiff = 0; ntot = 0
        firstpos = Int[]
        for (ki, sk) in enumerate(keys)
            for j in 1:nsig
                st = shake256(codeunits("rm/$ki/$j"), 56)
                setround(FE_TONEAREST)
                r1 = chacha20(st); empty!(F.MU); F.ON[] = true; F.CALLN[] = 0
                _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
                A = copy(F.MU); F.ON[] = false
                setround(mode)
                r2 = chacha20(st); empty!(F.MU); F.ON[] = true; F.CALLN[] = 0
                _, b = F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
                B = copy(F.MU); F.ON[] = false
                setround(FE_TONEAREST)
                ntot += 1
                if Int.(a) != Int.(b)
                    ndiff += 1
                    # where did the sampler first see a different centre that
                    # actually changed the floor?
                    for i in 1:min(length(A), length(B))
                        if floor(A[i]) != floor(B[i]); push!(firstpos, i); break; end
                    end
                end
            end
        end
        @printf("%-14s : %6d of %6d signatures differ   rate %.4g\n", name, ndiff, ntot, ndiff/ntot)
        if !isempty(firstpos)
            f2 = count(<=(2), firstpos); l2 = count(>=(L-1), firstpos)
            @printf("%-14s   first floor-flip position: %d at calls 1-2, %d at calls 2n-1..2n, %d elsewhere (of %d)\n",
                    "", f2, l2, length(firstpos)-f2-l2, length(firstpos))
        end
    end
    setround(FE_TONEAREST)
    println()
    println("# Reference for scale: the spelling differences measured elsewhere in this")
    println("# project diverge at 1.9e-5 (A1) and 7.5e-6 (A2).")
    println()

    # ---- a LATE flip: the mode changes only near the end of the traversal ----
    println("# A late flip: the adversary changes the rounding mode only after call k,")
    println("# so only the tail of the tree traversal is perturbed.  If the resulting")
    println("# difference is confined to the last two sampled integers, it is the")
    println("# STRUCTURED difference that ePrint 2024/1709 section 5.1 turns into a key.")
    for k in (L-8, L-4, L-2, L-1)
        conf = 0; div = 0; tot = 0
        for (ki, sk) in enumerate(keys)
            for j in 1:min(nsig, 60)
                st = shake256(codeunits("rmlate/$ki/$j"), 56)
                setround(FE_TONEAREST)
                F.FLIP_AT[] = 0; F.CALLN[] = 0
                r1 = chacha20(st)
                _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
                setround(FE_TONEAREST)
                F.FLIP_AT[] = k; F.FLIP_MODE[] = FE_UPWARD; F.CALLN[] = 0
                r2 = chacha20(st)
                _, b = F.sample_preimage(sk, pt, x -> randombytes!(r2, x))
                setround(FE_TONEAREST); F.FLIP_AT[] = 0
                tot += 1
                d = count(Int.(a) .!= Int.(b))
                if d > 0
                    div += 1
                    # structured == the difference is a 2-sparse polynomial
                    nz = findall(!=(0), Int.(a) .- Int.(b))
                    length(nz) <= 2 && (conf += 1)
                end
            end
        end
        @printf("  flip at call %4d of %d : %4d of %4d diverge, %4d of them 2-sparse (exploitable)\n",
                k, L, div, tot, conf)
    end
    setround(FE_TONEAREST)
end

main()
