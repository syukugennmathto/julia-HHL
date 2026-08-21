using Falcon, Printf

# Same LCG as same.c, so both sides multiply the same polynomials and the
# checksum is comparable.  If the checksums differ the comparison is void.
function inputs(n)
    s = UInt32(12345)
    f = Vector{Int}(undef, n); g = Vector{Int}(undef, n)
    for i in 1:n
        s = s * UInt32(1103515245) + UInt32(12345); f[i] = Int((s >> 8) % 12289)
        s = s * UInt32(1103515245) + UInt32(12345); g[i] = Int((s >> 8) % 12289)
    end
    return f, g
end

function main()
    for (n, iters) in ((128, 2000), (256, 1000), (512, 300), (1024, 100))
    f, g = inputs(n)
    o = Falcon.polymulq(f, g)
    ts = Float64[]
    for _ in 1:iters
        t0 = time_ns(); Falcon.polymulq(f, g); push!(ts, (time_ns() - t0) / 1e6)
    end
    sort!(ts)
    @printf("Jl  polymulq n=%4d: median %.4f ms   (checksum %d)\n",
            n, ts[iters ÷ 2 + 1], o[1] + o[n])
    end
end
main()
