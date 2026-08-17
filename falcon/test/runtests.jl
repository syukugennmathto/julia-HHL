using Test
using Falcon

# --- small helpers shared by the test files ---------------------------------

"Trial division primality test; only used on q, so speed is irrelevant."
function isprime_trial(n::Integer)
    n < 2 && return false
    n % 2 == 0 && return n == 2
    d = 3
    while d * d <= n
        n % d == 0 && return false
        d += 2
    end
    return true
end

"Evaluate a polynomial at a complex point by Horner, for property checks."
function _evalpoly_c(f::AbstractVector, z::Complex)
    acc = zero(z)
    for c in reverse(f)
        acc = acc * z + ComplexF64(c)
    end
    return acc
end

@testset "Falcon.jl" begin
    include("test_params.jl")   # module 1
    include("test_shake.jl")    # module 2
    include("test_poly.jl")     # module 3
end
