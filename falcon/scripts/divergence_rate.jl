#!/usr/bin/env julia
#
# divergence_rate.jl -- do the specification's underdetermined choices reach
# the signature, and how often?
#
#     julia --project=falcon falcon/scripts/divergence_rate.jl [keys] [sigs_per_key]
#
# ===========================================================================
# THE FOUR ARMS, AND WHY THEY ARE SEPARATE
# ===========================================================================
#
# Two implementations of Falcon can both follow the specification and still
# compute different intermediate values, because the specification does not
# determine them.  This project found four such places, in three groups.  They
# are NOT equally interesting, and the whole point of running them separately
# is that lumping them together hid which was which (docs/debug_log.md #055).
#
#   A1. THE HAND-UNROLLED BOTTOM LEVELS of ffSampling.  The C reference
#       special-cases `logn == 2` and folds the split/merge twiddles into
#       `1/sqrt(2)` and `1/sqrt(8)`, one multiplication where the generic
#       recursion does two.
#
#       **This is already published.**  ePrint 2024/1709 (Lin, Tibouchi, Yu,
#       Zhang, EUROCRYPT 2025) section 6.1 identifies exactly this difference,
#       as the gap between the reference's own `sign_dyn` and `sign_tree`
#       modes, and writes out the same two formulas we do (their p.22).
#       Section 7.2 gives a countermeasure.  Running A1 reproduces a known
#       result; it does not establish a new one.
#
#   A2. COMPLEX DIVISION, and THE `D11` ENTRY OF LDL*.  The C reference forms
#       the reciprocal explicitly (`m = 1/(b_re^2+b_im^2)`) where the
#       specification gives no formula and its readers use their language's
#       complex division (Smith's algorithm, in both Python and Julia); and it
#       computes `D11 = G11 - mu*adj(G01)` where Algorithm 8 and the Python
#       reference compute `G11 - L10*adj(L10)*G00`.
#
#       These do NOT appear in 2024/1709.  They are differences between the
#       specification and the C reference rather than between two entry points
#       of the C reference, so they are present in BOTH of C's signing modes
#       and are untouched by that paper's countermeasure.  They live in the
#       tree: they perturb `l10` and the leaf widths, not the unrolled block.
#
#   B.  A CONSTANT THAT CANNOT BE DERIVED.  `fpr_inv_sigma[logn]` is not the
#       correctly rounded reciprocal of the sigma Table 3.3 publishes
#       (scripts/inv_sigma_audit.jl).  This perturbs the sampler's WIDTH.
#
#       **Also already settled**, and by a theorem rather than a measurement:
#       Lemma 2 of 2024/1709 bounds the probability of an inconsistent
#       execution by 160*eps for |sigma - sigma'| <= eps.  At eps ~ 2^-52 that
#       is ~1e-11 per signature.  No feasible experiment can see it, so the
#       zero this arm reports is not evidence of anything (#055).
#
# A1 and A2 both reach the sampler's CENTRE; B reaches its WIDTH.  That
# distinction is Lemma 1 against Lemma 2 of the same paper.
#
# ===========================================================================
# THE MECHANISM
# ===========================================================================
#
# `SamplerZ` begins `s = floor(mu)`.  `floor` is discontinuous, so a
# discrepancy of 1e-13 in `mu` becomes a difference of 1 in `s` whenever the
# two evaluations land on opposite sides of an integer.  That changes `z`,
# changes how many bytes `BerExp` consumes, and desynchronises the PRNG, so
# the rest of the signature is unrelated.
#
# Integer centres are not uniformly likely across the traversal.  Heuristic 1
# of 2024/1709 puts the probability at 1/q for the first two calls and
# 1/||(g,-f)||^2 for the last two, and negligible elsewhere.  We measured that
# directly and it holds (scripts/heuristic1_check.jl):
#
#     calls 1, 2        17 near-integer centres in 200000   (heuristic 8.1e-5)
#     calls 2n-1, 2n    14 near-integer centres in 200000   (heuristic 6.1e-5)
#     calls 3 .. 2n-2    0 in 102000000
#
# A divergence in the LAST two calls is the dangerous one: 2024/1709 section 5
# recovers the whole private key from a single such pair.  A divergence in the
# first two yields only a short lattice vector, which is not enough.
#
# ===========================================================================
# WHAT TO READ OFF THE OUTPUT
# ===========================================================================
#
# The question this script now exists to answer is A2, and only A2:
#
#   * if A2 diverges, then an implementer who reads only the standard
#     disagrees with the reference in a way that 2024/1709's countermeasure
#     does not fix, and the standard's text has to pin the two formulas;
#   * if A2 does not diverge, then everything measured here is 2024/1709's,
#     and what this project has to say about it is a reproduction.
#
# A1 is kept as a positive control -- it should diverge, at the rate that
# paper's Table 2 reports (about 3e-5; we measure 5e-5) -- and B as a
# negative control whose expected value is zero for a known reason.

using Falcon
using Printf

const F = Falcon

point() = hash_to_point(collect(codeunits("power analysis message")),
                        shake256(codeunits("power salt"), SALT_LEN), 512;
                        q = FALCON_512.q)

"""
Sign `nsig` times per key with the two spellings on identical randomness, and
count how often the results differ.

`build` is called with a private key and must return the private key to sign
with in the second run; `sign2` wraps the second run's signing call.  Between
them these cover both places a respelling can live -- the tree (A2, B) and the
sampling itself (A1).
"""
function rate(nkeys::Int, nsig::Int, tag::String; build = identity, sign2 = f -> f())
    pt = point()
    r = chacha20(collect(UInt8, 0x00:0x37))
    n = 0; ndiff = 0; ncoef = 0; ncoefdiff = 0
    divergent = Tuple{Int,Int,Int}[]
    for ki in 1:nkeys
        sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
        sk2 = build(sk)
        for j in 1:nsig
            st = shake256(codeunits("$tag/$ki/$j"), 56)
            r1 = chacha20(st); r2 = chacha20(st)
            _, a = F.sample_preimage(sk, pt, x -> randombytes!(r1, x))
            _, b = sign2(() -> F.sample_preimage(sk2, pt, x -> randombytes!(r2, x)))
            n += 1
            d = count(Int.(a) .!= Int.(b))
            ncoef += length(a); ncoefdiff += d
            d > 0 && (ndiff += 1; push!(divergent, (ki, j, d)))
        end
    end
    return (n, ndiff, ncoef, ncoefdiff, divergent)
end

"""
Rebuild the ffLDL tree from the same basis with `cdiv` and/or `ldl` spelled the
specification's way, leaving the KEY ITSELF untouched.

The key must not move: `div_fft` is also used inside key generation (the Babai
reduction in ntrugen.jl), so toggling it there would produce a different
(f, g, F, G) and the comparison would be meaningless.
"""
function respelled_tree(sk; cdiv::Bool, ldl::Bool)
    t = with_spec_spelling(cdiv = cdiv, ldl = ldl) do
        F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    end
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

"Replace every leaf width by one built from `1/sigma` instead of the table."
function recip_tree(sk)
    recip = 1 / sk.params.sigma
    t = F.normalize_tree!(F.ffldl_fft(F.gram_fft(sk.B0_fft)), sk.params.sigma)
    patch!(x) = x isa F.FFLDLNode ? (patch!(x.left); patch!(x.right)) :
                (x.isigma = sqrt(real(x.value[1])) * recip)
    patch!(t)
    return F.FalconPrivateKey(sk.params, sk.f, sk.g, sk.F, sk.G, sk.B0_fft, t)
end

leafisig(t::F.FFLDLLeaf) = [t.isigma]
leafisig(t::F.FFLDLNode) = vcat(leafisig(t.left), leafisig(t.right))

"Every Float64 in the tree: the `l10` vectors of the nodes and the leaf widths."
treedoubles(t::F.FFLDLLeaf) = [t.isigma]
treedoubles(t::F.FFLDLNode) =
    vcat(reinterpret(Float64, t.l10), treedoubles(t.left), treedoubles(t.right))

"""
How many of the tree's doubles the respelling actually moved.

**A zero here would mean the arm is a no-op** and its zero divergence rate
says nothing -- which is the failure mode this check exists to catch.  The
harness in this project has silently measured nothing before
(docs/debug_log.md #043, and the dudect controls in scripts/dudect.jl).
"""
function treediff(sk, sk2)
    a = treedoubles(sk.tree); b = treedoubles(sk2.tree)
    length(a) == length(b) || error("trees differ in shape: $(length(a)) vs $(length(b))")
    return (count(i -> reinterpret(UInt64, a[i]) != reinterpret(UInt64, b[i]),
                  eachindex(a)), length(a))
end

function report(title::String, res, note::String = "")
    n, nd, nc, ncd, div = res
    println("## ", title)
    isempty(note) || println("   ", note)
    @printf("  signatures differing   : %d of %d   (rate %.3g)\n", nd, n, nd / n)
    @printf("  coefficients differing : %d of %d\n", ncd, nc)
    if nd == 0
        @printf("  zero events -> 95%% upper bound on the rate: %.3g\n", 3 / n)
    else
        @printf("  when one differs, %.0f of 512 coefficients do\n", ncd / nd)
        println("  divergent cases (key, signature, coefficients):")
        for d in div
            @printf("    %s\n", string(d))
        end
    end
    println()
end

function main()
    nkeys = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    nsig  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1000
    println("# divergence_rate.jl  n=512  ", nkeys, " keys x ", nsig, " signatures per arm")
    println()

    # Before measuring anything: confirm each tree-side arm actually perturbs
    # the tree.  A no-op arm would report a zero divergence rate that looks
    # exactly like a real negative result.
    let r0 = chacha20(collect(UInt8, 0x00:0x37)),
        sk0 = falcon_keygen(512, k -> randombytes!(r0, k))[1]
        println("# sanity: doubles moved in the ffLDL tree by each perturbation")
        for (name, mk) in (("A2  cdiv+ldl", sk -> respelled_tree(sk, cdiv = true,  ldl = true)),
                           ("A2a cdiv    ", sk -> respelled_tree(sk, cdiv = true,  ldl = false)),
                           ("A2b ldl     ", sk -> respelled_tree(sk, cdiv = false, ldl = true)),
                           ("B   1/sigma ", recip_tree))
            d, tot = treediff(sk0, mk(sk0))
            @printf("  %s : %5d of %5d%s\n", name, d, tot,
                    d == 0 ? "   <-- NO-OP, this arm measures nothing" : "")
        end
        println()
    end

    report("A2 -- complex division and D11, the tree only  (THE OPEN QUESTION)",
           rate(nkeys, nsig, "A2"; build = sk -> respelled_tree(sk, cdiv = true, ldl = true)),
           "not in ePrint 2024/1709; present in both of C's signing modes")

    report("A1 -- the hand-unrolled bottom levels  (positive control)",
           rate(nkeys, nsig, "A1"; sign2 = with_spec_ffsampling),
           "= ePrint 2024/1709 section 6.1, sign_dyn vs sign_tree; Table 2 reports ~3e-5")

    # split A2, so a divergence can be attributed
    report("A2a -- complex division alone",
           rate(nkeys, nsig, "A2a"; build = sk -> respelled_tree(sk, cdiv = true, ldl = false)))
    report("A2b -- the D11 spelling alone",
           rate(nkeys, nsig, "A2b"; build = sk -> respelled_tree(sk, cdiv = false, ldl = true)))

    # B, with the leaf count that shows the perturbation really is applied
    pt = point()
    r = chacha20(collect(UInt8, 0x00:0x37))
    sk = falcon_keygen(512, k -> randombytes!(r, k))[1]
    nleaf = count(zip(leafisig(sk.tree), leafisig(recip_tree(sk).tree))) do (x, y)
        reinterpret(UInt64, x) != reinterpret(UInt64, y)
    end
    report("B -- fpr_inv_sigma vs 1/sigma  (negative control)",
           rate(nkeys, nsig, "B"; build = recip_tree),
           "$nleaf of 512 tree leaves differ; Lemma 2 of 2024/1709 predicts ~1e-11, " *
           "so zero here is expected and carries no information")
end

main()
