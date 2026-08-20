#!/usr/bin/env julia
#
# Put the C reference's timings and ours side by side.
#
#     julia falcon/scripts/bench_compare.jl c_emu.txt c_native.txt julia.txt
#
# Each argument is the output of `cref_bench` or of `bench.jl` -- they share a
# format on purpose:
#
#     # <a comment line naming the build>
#     # op n median_ms mean_ms min_ms iters
#     keygen 512 14.421030 16.608956 10.855031 20
#
# The first `#` line of each file is used as its column label, so a run is
# always labelled with which build produced it.  Reporting a FALCON timing
# without saying whether the floating point was emulated is meaningless: the
# two differ by more than 10x on signing (docs/debug_log.md #031, #032).
#
# The operation names are already aligned by the two drivers: our `falcon_sign`
# takes an expanded key, so it lines up with C's `sign_tree`, not `sign_dyn`.

using Printf

const ORDER = ["keygen", "expand_privkey", "sign_dyn", "sign_tree", "verify"]

function read_bench(path::AbstractString)
    label = basename(path)
    rows = Dict{Tuple{String,Int},NamedTuple}()
    for (i, line) in enumerate(eachline(path))
        if startswith(line, "#")
            i == 1 && (label = strip(lstrip(line, ['#', ' '])))
            continue
        end
        isempty(strip(line)) && continue
        f = split(line)
        length(f) == 6 || continue
        rows[(String(f[1]), parse(Int, f[2]))] =
            (median = parse(Float64, f[3]), mean = parse(Float64, f[4]),
             min = parse(Float64, f[5]), iters = parse(Int, f[6]))
    end
    return label, rows
end

fmt(x) = x >= 1000 ? @sprintf("%.0f", x) :
         x >= 10   ? @sprintf("%.1f", x) :
                     @sprintf("%.3f", x)

function main(paths)
    isempty(paths) && (println("usage: bench_compare.jl FILE..."); return)
    runs = [read_bench(p) for p in paths]

    ns = sort(unique(n for (_, rows) in runs for (_, n) in keys(rows)))
    for n in ns
        println("\n## n = ", n, "  (median ms per call)\n")
        print("| op |")
        for (label, _) in runs
            print(" ", label, " |")
        end
        # A ratio column only makes sense against a single baseline; use the
        # last file, which by the usage above is ours.
        println(" ratio (last / first) |")
        print("|:---|")
        for _ in runs
            print("---:|")
        end
        println("---:|")

        for op in ORDER
            any(haskey(rows, (op, n)) for (_, rows) in runs) || continue
            print("| `", op, "` |")
            for (_, rows) in runs
                r = get(rows, (op, n), nothing)
                print(r === nothing ? " -- |" : " " * fmt(r.median) * " |")
            end
            # The ratio compares the first file against the last, and only when
            # *both* have the row.  Falling back to "the last file that has it"
            # would silently compare two C builds to each other under a heading
            # that says otherwise -- `sign_dyn` has no Julia counterpart, and a
            # 0.1x there would read as a Julia win.
            a = get(runs[1][2], (op, n), nothing)
            b = get(runs[end][2], (op, n), nothing)
            if a !== nothing && b !== nothing && a.median > 0
                @printf(" %.1fx |\n", b.median / a.median)
            else
                println(" n/a |")
            end
        end
    end
    println()
end

main(ARGS)
