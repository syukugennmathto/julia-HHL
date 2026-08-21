using Falcon, Printf
import Falcon: bitsize, karamul, field_norm, lift, galois_conjugate, xgcd_floor, NTRUSolveFailure

# A copy of babai_reduce that reports, per level, how many iterations it ran and
# *why* it stopped.  The two exits are very different findings:
#   "Size < size"  -- the reduction finished
#   "k rounded to 0" -- it gave up with (F, G) still enormous
function babai_watch(f, g, F, G, n)
    F = copy(F); G = copy(G)
    size_ = max(53, maximum(bitsize, f), maximum(bitsize, g))
    fa = Float64[Float64(c >> (size_ - 53)) for c in f]
    ga = Float64[Float64(c >> (size_ - 53)) for c in g]
    fa_fft = fft(fa); ga_fft = fft(ga)
    den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)), mul_fft(ga_fft, adj_fft(ga_fft)))
    iters = 0; reason = "?"; kmax_first = 0.0
    while true
        Size = max(53, maximum(bitsize, F), maximum(bitsize, G))
        if Size < size_; reason = "Size < size (finished)"; break; end
        Fa = Float64[Float64(c >> (Size - 53)) for c in F]
        Ga = Float64[Float64(c >> (Size - 53)) for c in G]
        num_fft = add_fft(mul_fft(fft(Fa), adj_fft(fa_fft)), mul_fft(fft(Ga), adj_fft(ga_fft)))
        k = ifft(div_fft(num_fft, den_fft))
        iters == 0 && (kmax_first = maximum(abs, real.(k)))
        ki = BigInt[BigInt(round(elt)) for elt in k]
        if all(iszero, ki); reason = "k rounded to 0 (gave up)"; break; end
        iters += 1
        fk = karamul(f, ki); gk = karamul(g, ki)
        shift = Size - size_
        for i in 1:n
            F[i] -= fk[i] << shift; G[i] -= gk[i] << shift
        end
    end
    return (F, G, iters, reason, size_, kmax_first)
end

function solve_watch(f::Vector{BigInt}, g::Vector{BigInt}, q, rows)
    n = length(f)
    if n == 1
        d, u, v = xgcd_floor(f[1], g[1])
        d == 1 || throw(NTRUSolveFailure("gcd != 1"))
        return (BigInt[-q * v], BigInt[q * u])
    end
    Fp, Gp = solve_watch(field_norm(f), field_norm(g), q, rows)
    F = karamul(lift(Fp), galois_conjugate(g))
    G = karamul(lift(Gp), galois_conjugate(f))
    raw = max(maximum(bitsize, F), maximum(bitsize, G))
    Fo, Go, iters, reason, size_, kmax = babai_watch(f, g, F, G, n)
    out = max(maximum(bitsize, Fo), maximum(bitsize, Go))
    push!(rows, (n, maximum(bitsize, f), size_, raw, out, iters, kmax, reason))
    return (Fo, Go)
end

function main()
    src = ReplayBytes(shake256(codeunits("why/stall"), 1 << 23))
    local f, g
    while true
        f = gen_poly(512, src); g = gen_poly(512, src)
        gs_norm(Float64.(f), Float64.(g)) > gram_schmidt_quality()^2 * 12289 && continue
        is_invertible_zq(Int.(f)) || continue
        try; ntru_solve(f, g); break; catch e; e isa NTRUSolveFailure || rethrow(); end
    end
    rows = Any[]
    F, G = solve_watch(BigInt.(f), BigInt.(g), 12289, rows)
    @assert Falcon.ntru_equation_holds(BigInt.(f), BigInt.(g), F, G)
    @printf("%6s %9s %7s %9s %9s %7s %11s  %s\n",
            "n", "bits(f,g)", "size", "raw(F,G)", "out(F,G)", "iters", "max|k| 1st", "why it stopped")
    println("-"^92)
    for r in sort(rows, by = x -> -x[1])
        @printf("%6d %9d %7d %9d %9d %7d %11.3g  %s\n", r...)
    end
end
main()
