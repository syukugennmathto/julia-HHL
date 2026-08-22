#!/usr/bin/env julia
#
# event_solve.jl -- recover f from integer-centre EVENTS alone, by linear
# algebra over F_q: no difference vectors, no lattice reduction.
#
#     julia --project=falcon falcon/scripts/event_solve.jl [n] [trials]
#
# ePrint 2024/1709 section 5 recovers the key from the difference VECTOR of two
# signatures that diverge at the LAST two sampler calls, and dismisses a
# divergence at the FIRST two: "a short lattice vector, but it is not expected
# to be short enough to make key recovery feasible."
#
# Section 6.4 of this paper shows the first-two centre is exactly
# (c*f)[n/2-1]/q, a PUBLIC LINEAR FORM in the secret.  So an integer centre at
# call 1 is exactly the condition
#
#     (c*f)[n/2-1] = 0   (mod q),
#
# one F_q-linear equation on f with coefficients the adversary can compute from
# the public message.  This script shows what that buys: given n-1 independent
# events, f is the kernel of the resulting system, and Gaussian elimination
# recovers it exactly.  The difference vectors are never used, and no lattice
# reduction is involved -- which is precisely the step section 5 assumed was
# necessary and infeasible.
#
# What this script does and does not establish.  It establishes the SOLVE: that
# events determine the key, cheaply and exactly.  It does not supply the oracle
# that tells the adversary which messages produced an event; section 6.4 states
# that requirement plainly and scripts/window.jl characterises when a given
# implementation difference provides one.  The syndromes here are drawn
# uniformly rather than hashed, which is what hash_to_point produces and what
# the measured event rate confirms (2.08e-4 observed against 1.63e-4 predicted,
# scripts/first_two_probe.jl); using real hashes would only cost SHAKE time.

using Falcon
using Printf
using Random

const F = Falcon

"Row of the linear form: coefficient of f_j in (c*f)[k], negacyclic."
function form_row(c::Vector{Int}, k::Int, n::Int, q::Int)
    row = Vector{Int}(undef, n)
    for j in 0:(n-1)
        idx = k - j
        row[j+1] = idx >= 0 ? c[idx+1] : mod(-c[idx+n+1], q)
    end
    return Int[mod(x, q) for x in row]
end

"Gaussian elimination mod q; returns (rank, row-reduced rows, pivot columns)."
function rref_modq(rows::Vector{Vector{Int}}, n::Int, q::Int)
    A = [copy(r) for r in rows]
    piv = Int[]
    r = 1
    for col in 1:n
        p = findfirst(i -> A[i][col] % q != 0, r:length(A))
        p === nothing && continue
        p += r - 1
        A[r], A[p] = A[p], A[r]
        inv = invmod(A[r][col], q)
        A[r] = Int[mod(x * inv, q) for x in A[r]]
        for i in 1:length(A)
            i == r && continue
            f = A[i][col]
            f == 0 && continue
            A[i] = Int[mod(A[i][t] - f * A[r][t], q) for t in 1:n]
        end
        push!(piv, col); r += 1
        r > length(A) && break
    end
    return r - 1, A, piv
end

"The unique (up to scale) kernel vector of a rank n-1 system."
function kernel_vector(A, piv, n, q)
    free = setdiff(1:n, piv)
    length(free) == 1 || return nothing
    fc = free[1]
    v = zeros(Int, n); v[fc] = 1
    for (r, c) in enumerate(piv)
        v[c] = mod(-A[r][fc], q)
    end
    return v
end

centre(x, q) = (y = mod(x, q); y > q ÷ 2 ? y - q : y)

function main()
    n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 512
    p = n == 512 ? FALCON_512 : FALCON_1024
    q = p.q
    rng = MersenneTwister(20260822)

    r = chacha20(collect(UInt8, 0x00:0x37))
    sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
    f = Int.(sk.f)
    @printf("# event_solve.jl  n=%d  q=%d\n", n, q)
    @printf("# target f: max |f_i| = %d\n\n", maximum(abs, f))

    k = n ÷ 2 - 1                     # the coefficient the first call reads
    rows = Vector{Int}[]
    tried = 0; events = 0
    rank = 0; A = nothing; piv = nothing
    while rank < n - 1
        tried += 1
        c = rand(rng, 0:(q-1), n)                    # a syndrome
        row = form_row(c, k, n, q)
        val = mod(sum(row[j] * f[j] for j in 1:n), q)
        val == 0 || continue                          # not an event
        events += 1
        push!(rows, row)
        if events % 128 == 0 || events >= n - 1
            rank, A, piv = rref_modq(rows, n, q)
        end
        events > n + 40 && break
    end
    @printf("collected %d events from %d syndromes (rate %.3g, predicted 1/q = %.3g)\n",
            events, tried, events / tried, 1 / q)
    @printf("system rank %d of %d unknowns\n", rank, n)

    v = kernel_vector(A, piv, n, q)
    if v === nothing
        println("kernel is not one-dimensional; collect more events")
        return
    end
    # scale the kernel vector so its entries are small: f is short
    best = nothing; bestmax = typemax(Int)
    for s in 1:(q-1)
        w = Int[centre(s * v[i], q) for i in 1:n]
        m = maximum(abs, w)
        if m < bestmax; bestmax = m; best = w; end
        bestmax <= maximum(abs, f) && break
    end
    ok = best == f || best == .-f
    @printf("\nrecovered vector: max |.| = %d   (true f has max |.| = %d)\n", bestmax, maximum(abs, f))
    @printf("matches the secret f: %s%s\n", ok, best == .-f ? "  (up to sign)" : "")
    if ok
        println("\n*** KEY RECOVERED FROM EVENTS ALONE ***")
        println("    No signature difference vectors were used.  No lattice reduction.")
        println("    Only the boolean 'this message produced an integer centre at call 1',")
        println("    plus the public syndrome, for ", events, " messages.")
        println("    ePrint 2024/1709 section 5 assumed a first-two divergence could not")
        println("    yield the key; over F_q it does, and in ", n, " x ", n, " Gaussian elimination.")
    end
end

main()
