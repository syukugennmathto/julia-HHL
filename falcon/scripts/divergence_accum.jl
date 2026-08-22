#!/usr/bin/env julia
#
# divergence_accum.jl -- one INDEPENDENT chunk of a divergence measurement,
# appended to a results file, so a long measurement can be built up across
# many short foreground runs that each survive a session suspension.
#
#     julia --project=falcon falcon/scripts/divergence_accum.jl \
#           <arm> <logn> <nkeys> <nsig> <chunk_id> <outfile>
#
#   arm    : A1 | A2 | A2a | A2b
#   logn   : 9 (n=512) or 10 (n=1024)
#   chunk_id: any integer; it seeds BOTH the keys and the per-signature PRNG,
#             so distinct chunk_ids are independent samples that may be summed.
#
# Appends one line to <outfile>:
#     arm logn chunk nkeys nsig events ncoef ncoefdiff  [ (ki,j,d) ... ]
# and prints the running total across all lines already in the file.
#
# Why chunks instead of one long run: background jobs in this environment are
# killed by SIGTERM when the session suspends between turns (docs/debug_log.md
# #057).  A foreground run completes within its turn, and persisting each chunk
# means a kill costs at most the chunk in flight, not the whole measurement.
# Summing independent Poisson chunks is exact: total events over total n.

using Falcon
using Printf

const F = Falcon

arm    = ARGS[1]
logn   = parse(Int, ARGS[2])
nkeys  = parse(Int, ARGS[3])
nsig   = parse(Int, ARGS[4])
chunk  = parse(Int, ARGS[5])
out    = ARGS[6]

n = 1 << logn
params = logn == 9 ? FALCON_512 : FALCON_1024

pt = hash_to_point(collect(codeunits("power analysis message")),
                   shake256(codeunits("power salt"), SALT_LEN), n; q = params.q)

# Distinct key stream per chunk: fold chunk id into the keygen seed.
keyseed = shake256(codeunits("accum-keys/$logn/$chunk"), 56)

function respelled_tree(sk; cdiv::Bool, ldl::Bool)
    t = with_spec_spelling(cdiv = cdiv, ldl = ldl) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

build, sign2 =
    arm == "A1"  ? (identity,                                    with_spec_ffsampling) :
    arm == "A2"  ? (sk -> respelled_tree(sk, cdiv=true,  ldl=true),  f -> f()) :
    arm == "A2a" ? (sk -> respelled_tree(sk, cdiv=true,  ldl=false), f -> f()) :
    arm == "A2b" ? (sk -> respelled_tree(sk, cdiv=false, ldl=true),  f -> f()) :
    error("unknown arm $arm")

function measure(build, sign2, keyseed)
    r = chacha20(keyseed)
    events = 0; ncoef = 0; ncoefdiff = 0
    divergent = Tuple{Int,Int,Int}[]
    for ki in 1:nkeys
        sk = falcon_keygen(n, k -> randombytes!(r, k))[1]
        sk2 = build(sk)
        for j in 1:nsig
            st = shake256(codeunits("accum/$arm/$logn/$chunk/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            _, a = F.sample_preimage(sk,  pt, x -> randombytes!(r1, x))
            _, b = sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
            d = count(Int.(a) .!= Int.(b))
            ncoef += length(a); ncoefdiff += d
            d > 0 && (events += 1; push!(divergent, (ki, j, d)))
        end
    end
    return events, ncoef, ncoefdiff, divergent
end

events, ncoef, ncoefdiff, divergent = measure(build, sign2, keyseed)

open(out, "a") do io
    @printf(io, "%s %d %d %d %d %d %d %d  %s\n",
            arm, logn, chunk, nkeys, nsig, events, ncoef, ncoefdiff,
            isempty(divergent) ? "-" : join(string.(divergent), " "))
end

# running total across the whole file
function running_total(out, arm, logn)
    tot_ev = 0; tot_n = 0
    for line in eachline(out)
        f = split(line)
        length(f) >= 8 || continue
        (f[1] == arm && parse(Int, f[2]) == logn) || continue
        tot_n  += parse(Int, f[4]) * parse(Int, f[5])
        tot_ev += parse(Int, f[6])
    end
    return tot_ev, tot_n
end

tot_ev, tot_n = running_total(out, arm, logn)
@printf("chunk done: %s n=%d chunk=%d -> %d events in %d sigs\n",
        arm, n, chunk, events, nkeys * nsig)
@printf("RUNNING TOTAL %s n=%d : %d events in %d signatures  (rate %.3g)\n",
        arm, n, tot_ev, tot_n, tot_ev / max(tot_n, 1))
