using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "julia_env"); io=devnull)

using Yao
using LinearAlgebra
using BenchmarkTools
using JSON
using Printf
using Statistics

function build_hhl_circuit()
    C = 2.0 / 3.0
    lambda1 = 2.0 / 3.0
    lambda2 = 4.0 / 3.0
    theta1 = 2.0 * asin(C / lambda1)
    theta2 = 2.0 * asin(C / lambda2)
    chain(4,
        put(2 => H), put(3 => H),
        put(2 => shift(3π / 2)), put(3 => shift(π / 2)),
        control(2, 4 => X), put(4 => shift(-3π / 4)),
        control(2, 4 => X), put(4 => shift(3π / 4)),
        put(3 => H), control(2, 3 => shift(-π / 2)),
        put(2 => H), swap(2, 3),
        control(2, 1 => Ry(theta1)), control(3, 1 => Ry(theta2)),
        swap(2, 3),
        put(2 => H), control(2, 3 => shift(π / 2)), put(3 => H),
        control(3, 4 => X), put(4 => shift(-π / 4)),
        control(3, 4 => X), put(4 => shift(π / 4)),
        control(2, 4 => X), put(4 => shift(-3π / 4)),
        control(2, 4 => X), put(4 => shift(3π / 4)),
        put(2 => H), put(3 => H),
    )
end

function run_prebuilt(circuit)
    reg = zero_state(4)
    reg |> circuit
    return statevec(reg)
end

function end_to_end()
    circuit = build_hhl_circuit()
    reg = zero_state(4)
    reg |> circuit
    return statevec(reg)
end

function simulate_measurements(circuit, nshots::Int)
    reg = zero_state(4)
    reg |> circuit
    return measure(reg; nshots=nshots)
end

function main()
    println("=== Ch.5 HHL 細粒度ベンチマーク (Julia / Yao.jl) ===")

    prebuilt = build_hhl_circuit()

    # warmup (JITコンパイルを起こす)
    run_prebuilt(prebuilt)
    end_to_end()
    simulate_measurements(prebuilt, 10)

    function summarize(bench)
        Dict(
            "samples" => length(bench.times),
            "min_ms" => minimum(bench.times) / 1e6,
            "median_ms" => median(bench.times) / 1e6,
            "mean_ms" => mean(bench.times) / 1e6,
            "max_ms" => maximum(bench.times) / 1e6,
            "stdev_ms" => std(bench.times) / 1e6,
        )
    end

    println("  -- build only --")
    b_build = @benchmark build_hhl_circuit() samples=200 evals=1
    display(b_build)
    println()

    println("  -- run only (prebuilt) --")
    b_run = @benchmark run_prebuilt($prebuilt) samples=200 evals=1
    display(b_run)
    println()

    println("  -- end-to-end (build + run) --")
    b_e2e = @benchmark end_to_end() samples=200 evals=1
    display(b_e2e)
    println()

    println("  -- simulate + measure 4096 shots --")
    b_measure = @benchmark simulate_measurements($prebuilt, 4096) samples=100 evals=1
    display(b_measure)
    println()

    out = Dict(
        "build_only" => summarize(b_build),
        "run_only" => summarize(b_run),
        "end_to_end" => summarize(b_e2e),
        "measure_4096" => summarize(b_measure),
    )

    @printf "\nSummary (ms):\n"
    for (k, v) in out
        @printf "  %-22s min=%.4f med=%.4f mean=%.4f\n" k v["min_ms"] v["median_ms"] v["mean_ms"]
    end

    out_path = joinpath(@__DIR__, "..", "results", "julia_bench.json")
    open(out_path, "w") do io
        JSON.print(io, out, 2)
    end
    println("Saved -> ", out_path)
end

main()
