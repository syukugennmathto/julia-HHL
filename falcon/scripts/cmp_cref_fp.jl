#!/usr/bin/env julia
#
# Bisect the floating-point path against the C reference, bit for bit.
#
#     cd falcon/scripts/cref
#     cc -O2 -fPIC -shared -o /tmp/libfalcon.so \
#        codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c sign.c vrfy.c -lm
#     cc -O2 -fPIC -shared -o /tmp/libfalconshim.so ../cref_shim.c \
#        codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c vrfy.c -lm
#     julia --project=falcon falcon/scripts/cmp_cref_fp.jl /tmp/libfalcon.so /tmp/libfalconshim.so
#
# WHY BIT-EXACT AND NOT "CLOSE"
# -----------------------------
# FALCON's signature is a *rounded* function of these values: `samplerz`
# compares a fixed-point exponential against random bytes, so a difference of
# one ulp in its centre can flip an accept/reject and change the signature.
# "Agrees to 1e-12" therefore says nothing about whether the bytes will match.
# Every comparison below is on the raw 64-bit patterns.
#
# WHAT THIS FOUND (docs/debug_log.md #048)
# ----------------------------------------
# Bit-exact with C, at every degree from 8 to 1024, with no change needed:
#     fft  ifft  split_fft  merge_fft  add_fft  sub_fft  mul_fft  adj_fft
#     mul_fft(a, adj_fft(b))   mul_fft(a, adj_fft(a))
#     the B0 matrix of expand_privkey
# Bit-exact only after being respelled the way C spells them:
#     div_fft      (fft.jl: `_cdiv_cref` -- C forms 1/|b|^2 and multiplies,
#                   Julia and Python use Smith's algorithm)
#     ldl_fft      (ffsampling.jl -- C computes D11 = G11 - mu*adj(G01),
#                   the specification computes G11 - L10*adj(L10)*G00)
#     ffSampling   (ffsampling.jl: `_ffsampling_c4` -- C hand-unrolls the
#                   bottom two levels with the twiddles pre-folded into
#                   1/sqrt(2) and 1/sqrt(8); NOT wired in, see its docstring)
# Different by construction, and not yet reconciled:
#     the tree leaves.  C stores sqrt(x)*(1/sigma), we store sigma/sqrt(x) --
#     reciprocals of each other, and it is a storage convention, not an
#     arithmetic difference.  C's sampler takes `isigma` and forms
#     `dss = 0.5*isigma^2`, `ccs = sigma_min*isigma`, where the specification
#     forms `1/(2*sigma^2)` and `sigma_min/sigma`.  Those constants differ in
#     the last bit 55% and 32% of the time respectively -- and over 20000
#     draws on identical byte streams the sampled integer never differed.

using Falcon, Printf, Random, Libdl
const F = Falcon

length(ARGS) >= 1 || error("usage: cmp_cref_fp.jl <libfalcon.so> [<libfalconshim.so>]")
const LIB  = ARGS[1]
const SHIM = length(ARGS) >= 2 ? ARGS[2] : ""
const H = Libdl.dlopen(LIB)
sym(s) = Libdl.dlsym(H, s)

exact(a, b) = length(a) == length(b) &&
              all(reinterpret(UInt64, collect(Float64, a)) .== reinterpret(UInt64, collect(Float64, b)))
cfft(v, logn)  = (w = copy(v); ccall(sym(:falcon_inner_FFT),  Cvoid, (Ptr{Float64}, Cuint), w, logn); w)
cifft(v, logn) = (w = copy(v); ccall(sym(:falcon_inner_iFFT), Cvoid, (Ptr{Float64}, Cuint), w, logn); w)
c1(s, a, logn)    = (w = copy(a); ccall(sym(s), Cvoid, (Ptr{Float64}, Cuint), w, logn); w)
c2(s, a, b, logn) = (w = copy(a); ccall(sym(s), Cvoid, (Ptr{Float64}, Ptr{Float64}, Cuint), w, b, logn); w)

rng = MersenneTwister(20260821)
randf(n) = Float64[(rand(rng) - 0.5) * exp10(rand(rng) * 6 - 2) for _ in 1:n]

println("# transforms  (arbitrary doubles, not small integers)")
println("op                    ", join(lpad.("logn=" .* string.(3:10), 8)))
function row(name, f)
    @printf("%-20s  %s\n", name, join(lpad.([f(logn) ? "ok" : "DIFF" for logn in 3:10], 8)))
end
row("fft",       logn -> (n = 1 << logn; f = randf(n); exact(F.to_c_fft(fft(f)), cfft(f, logn))))
row("ifft",      logn -> (n = 1 << logn; Fc = cfft(randf(n), logn);
                          exact(real.(ifft(F.from_c_fft(Fc))), cifft(Fc, logn))))
row("split_fft", logn -> begin
        n = 1 << logn; A = cfft(randf(n), logn)
        j0, j1 = F.split_fft(F.from_c_fft(A))
        c0 = zeros(Float64, n >> 1); c1v = zeros(Float64, n >> 1)
        ccall(sym(:falcon_inner_poly_split_fft), Cvoid,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cuint), c0, c1v, A, logn)
        exact(F.to_c_fft(j0), c0) && exact(F.to_c_fft(j1), c1v)
    end)
row("merge_fft", logn -> begin
        n = 1 << logn; A = cfft(randf(n), logn)
        c0 = zeros(Float64, n >> 1); c1v = zeros(Float64, n >> 1)
        ccall(sym(:falcon_inner_poly_split_fft), Cvoid,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cuint), c0, c1v, A, logn)
        cm = zeros(Float64, n)
        ccall(sym(:falcon_inner_poly_merge_fft), Cvoid,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cuint), cm, c0, c1v, logn)
        exact(F.to_c_fft(merge_fft(F.from_c_fft(c0), F.from_c_fft(c1v))), cm)
    end)

println()
println("# FFT-domain arithmetic")
println("op                    ", join(lpad.("logn=" .* string.(3:10), 8)))
for (name, jf, cs) in (
        ("add_fft",    (a,b) -> add_fft(a,b),           :falcon_inner_poly_add),
        ("sub_fft",    (a,b) -> sub_fft(a,b),           :falcon_inner_poly_sub),
        ("mul_fft",    (a,b) -> mul_fft(a,b),           :falcon_inner_poly_mul_fft),
        ("div_fft",    (a,b) -> div_fft(a,b),           :falcon_inner_poly_div_fft),
        ("mul*adj",    (a,b) -> mul_fft(a, adj_fft(b)), :falcon_inner_poly_muladj_fft))
    res = map(3:10) do logn
        n = 1 << logn; A = cfft(randf(n), logn); B = cfft(randf(n), logn)
        exact(F.to_c_fft(jf(F.from_c_fft(A), F.from_c_fft(B))), c2(cs, A, B, logn)) ? "ok" : "DIFF"
    end
    @printf("%-20s  %s\n", name, join(lpad.(res, 8)))
end
for (name, jf, cs) in (
        ("adj_fft",    a -> adj_fft(a),             :falcon_inner_poly_adj_fft),
        ("mul selfadj", a -> mul_fft(a, adj_fft(a)), :falcon_inner_poly_mulselfadj_fft))
    res = map(3:10) do logn
        n = 1 << logn; A = cfft(randf(n), logn)
        exact(F.to_c_fft(jf(F.from_c_fft(A))), c1(cs, A, logn)) ? "ok" : "DIFF"
    end
    @printf("%-20s  %s\n", name, join(lpad.(res, 8)))
end

println()
println("# LDL*  (the specification's spelling vs the C reference's)")
res = map(3:10) do logn
    n = 1 << logn
    ff = F.from_c_fft(cfft(randf(n), logn)); gf = F.from_c_fft(cfft(randf(n), logn))
    Ff = F.from_c_fft(cfft(randf(n), logn)); Gf = F.from_c_fft(cfft(randf(n), logn))
    g00 = add_fft(mul_fft(ff, adj_fft(ff)), mul_fft(gf, adj_fft(gf)))
    g01 = add_fft(mul_fft(ff, adj_fft(Ff)), mul_fft(gf, adj_fft(Gf)))
    g11 = add_fft(mul_fft(Ff, adj_fft(Ff)), mul_fft(Gf, adj_fft(Gf)))
    G = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    G[1,1] = g00; G[2,1] = adj_fft(g01); G[1,2] = g01; G[2,2] = g11
    l10, d00, d11 = F.ldl_fft(G)
    a = F.to_c_fft(g00); b = F.to_c_fft(g01); c = F.to_c_fft(g11)
    ccall(sym(:falcon_inner_poly_LDL_fft), Cvoid,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cuint), a, b, c, logn)
    (exact(F.to_c_fft(l10), b) && exact(F.to_c_fft(d11), c)) ? "ok" : "DIFF"
end
@printf("%-20s  %s\n", "ldl_fft", join(lpad.(res, 8)))

println()
println("# expand_privkey: the B0 matrix and the ffLDL tree")
flat!(t::F.FFLDLNode, out) = (append!(out, F.to_c_fft(t.l10)); flat!(t.left, out); flat!(t.right, out))
flat!(t::F.FFLDLLeaf, out) = push!(out, t.sigma)
for logn in (9, 10)
    n = 1 << logn
    r = chacha20(collect(UInt8, 0x00:0x37))
    sk, _ = falcon_keygen(n, k -> randombytes!(r, k))
    ek = zeros(Float64, (8*logn + 40) * n ÷ 8)
    ccall(sym(:falcon_inner_expand_privkey), Cvoid,
          (Ptr{Float64}, Ptr{Int8}, Ptr{Int8}, Ptr{Int8}, Ptr{Int8}, Cuint, Ptr{UInt8}),
          ek, Int8.(sk.f), Int8.(sk.g), Int8.(sk.F), Int8.(sk.G), logn, zeros(UInt8, 48n))
    # C's layout is b00 = fft(g), b01 = fft(-f), b10 = fft(G), b11 = fft(-F),
    # which is our B0_fft read row-major.
    order = [sk.B0_fft[1,1], sk.B0_fft[1,2], sk.B0_fft[2,1], sk.B0_fft[2,2]]
    b = all(k -> exact(F.to_c_fft(order[k]), ek[((k-1)*n + 1):(k*n)]), 1:4)
    ours = Float64[]; flat!(sk.tree, ours); want = ek[(4n + 1):end]
    nd = count(i -> reinterpret(UInt64, ours[i]) != reinterpret(UInt64, want[i]), 1:length(ours))
    @printf("logn=%-3d B0 exact=%-6s  tree: %d of %d differ (= the %d leaves; storage convention)\n",
            logn, b, nd, length(ours), n)
end
