using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "julia_env"); io=devnull)

using Yao
using LinearAlgebra
using BenchmarkTools
using Printf
using JSON

# -----------------------------------------------------------------------------
# Ch.5: 基本操作
# -----------------------------------------------------------------------------
function ch5_basics()
    println("=== Yao.jl 基本 ===")
    h_circuit = chain(1, put(1 => H))
    reg = zero_state(1)
    reg |> h_circuit
    println("  H |0> => statevec = ", statevec(reg))
    println("  X gate matrix = ", mat(X))
    println("  H gate matrix = ", mat(H))
end

# -----------------------------------------------------------------------------
# HHL 回路 (Yao.jl 版 — 本のスニペットに準拠)
# Yao の qubit 番号は 1..N で、`put(k => op)` は k 番目に op を適用。
# 配線: 1=ancilla, 2,3=clock, 4=b
# -----------------------------------------------------------------------------
function build_hhl_circuit()
    n_total = 4
    C = 2.0 / 3.0
    lambda1 = 2.0 / 3.0
    lambda2 = 4.0 / 3.0
    theta1 = 2.0 * asin(C / lambda1)
    theta2 = 2.0 * asin(C / lambda2)

    circuit = chain(n_total,
        put(2 => H),
        put(3 => H),
        put(2 => shift(3π / 2)),
        put(3 => shift(π / 2)),
        control(2, 4 => X),
        put(4 => shift(-3π / 4)),
        control(2, 4 => X),
        put(4 => shift(3π / 4)),
        # 逆QFT
        put(3 => H),
        control(2, 3 => shift(-π / 2)),
        put(2 => H),
        swap(2, 3),
        # 制御回転
        control(2, 1 => Ry(theta1)),
        control(3, 1 => Ry(theta2)),
        # 逆QPE
        swap(2, 3),
        put(2 => H),
        control(2, 3 => shift(π / 2)),
        put(3 => H),
        control(3, 4 => X),
        put(4 => shift(-π / 4)),
        control(3, 4 => X),
        put(4 => shift(π / 4)),
        control(2, 4 => X),
        put(4 => shift(-3π / 4)),
        control(2, 4 => X),
        put(4 => shift(3π / 4)),
        put(2 => H),
        put(3 => H),
    )
    return circuit
end

function run_hhl_once()
    circuit = build_hhl_circuit()
    reg = zero_state(4)
    reg |> circuit
    statevec(reg)
end

"""
Yao のビット順: statevec の index i (1-origin) について、(i-1) を base2 展開すると
bit_1 (LSB) = qubit 1 (ancilla), bit_2 = qubit 2, ..., bit_N (MSB) = qubit N (b).
よって「アンシラ=1, クロック=00」の成分のみが |x> に比例する。
-> i-1 の下位 1 ビットが 1, 次2ビットが 00, 最上位1ビットが b の値。
"""
function extract_solution(sv)
    ancilla_1_clock_0 = ComplexF64[]
    for i in 1:length(sv)
        bits = i - 1
        ancilla = bits & 0b1
        clock = (bits >> 1) & 0b11
        if ancilla == 1 && clock == 0
            push!(ancilla_1_clock_0, sv[i])
        end
    end
    nrm = norm(ancilla_1_clock_0)
    if nrm < 1e-12
        return ancilla_1_clock_0
    end
    ancilla_1_clock_0 / nrm
end

function measure_ancilla_success(nshots::Int=4096)
    circuit = build_hhl_circuit()
    reg = zero_state(4)
    reg |> circuit
    sv = statevec(reg)
    # アンシラ = qubit 1 = bit 1 (LSB). P(anc=1) = sum_i |sv[i]|^2 where (i-1)&1==1.
    p = 0.0
    for i in 1:length(sv)
        if (i - 1) & 0b1 == 1
            p += abs2(sv[i])
        end
    end
    p
end

function main()
    ch5_basics()

    println("\n=== HHL 2x2 (Yao.jl) ===")
    sv = run_hhl_once()
    sol = extract_solution(sv)
    println("  量子解 |x> ∝ ", abs.(sol))

    A = [1.0 -1/3; -1/3 1.0]
    b = [1.0, 0.0]
    x_cl = A \ b
    x_cl_norm = x_cl / norm(x_cl)
    println("  古典解 x = ", x_cl, " (正規化後 ", x_cl_norm, ")")

    fidelity = abs(dot(sol, ComplexF64.(x_cl_norm)))^2
    println("  Fidelity |<x_q|x_cl>|^2 = ", fidelity)

    p_succ = measure_ancilla_success()
    println("  P(ancilla=1) = ", p_succ)

    # --- ベンチマーク (build + run)
    println("\n=== Benchmark (Julia / Yao.jl) ===")
    bench = @benchmark run_hhl_once() samples=50 evals=1
    println(bench)

    minimum_ms = minimum(bench.times) / 1e6
    median_ms = median(bench.times) / 1e6
    mean_ms = mean(bench.times) / 1e6
    maxval_ms = maximum(bench.times) / 1e6
    std_ms = std(bench.times) / 1e6
    @printf "  min=%.3fms, median=%.3fms, mean=%.3fms, max=%.3fms, std=%.3fms\n" minimum_ms median_ms mean_ms maxval_ms std_ms

    result = Dict(
        "ancilla_success_prob" => p_succ,
        "fidelity" => fidelity,
        "classical_solution" => collect(x_cl),
        "quantum_solution_abs" => abs.(sol),
        "benchmark_build_run" => Dict(
            "samples" => length(bench.times),
            "min_ms" => minimum_ms,
            "mean_ms" => mean_ms,
            "median_ms" => median_ms,
            "max_ms" => maxval_ms,
            "stdev_ms" => std_ms,
        ),
    )

    out_path = joinpath(@__DIR__, "..", "results", "julia_hhl.json")
    open(out_path, "w") do io
        JSON.print(io, result, 2)
    end
    println("Saved -> ", out_path)
end

main()
