#!/usr/bin/env julia
#
# inv_sigma_audit.jl -- is the C reference's `fpr_inv_sigma[]` the reciprocal
# of the sigma the specification publishes?
#
#     julia --project=falcon falcon/scripts/inv_sigma_audit.jl
#
# Short answer: no, and not of any sigma derivable from the specification
# either.  It is one ulp above at logn = 9 and two above at logn = 10, and the
# 35-digit decimal in fpr.h is not the expansion of 1/sigma to 35 digits -- it
# agrees to about sixteen significant figures and then diverges.
#
# WHY THIS MATTERS, AND HOW MUCH
# ------------------------------
# The reference multiplies every leaf of the ffLDL tree by this constant
# (`ffLDL_binary_normalize` in sign.c: "we actually store in the tree leaf the
# inverse of the value mandated by the specification").  An implementation that
# computes `1/sigma` from Table 3.3 instead of transcribing the table therefore
# builds a DIFFERENT expanded private key -- 444 of 512 leaves differ.
#
# It does not build different signatures.  Measured over 100000 signatures at
# n = 512 with identical key, message and PRNG state, not one of 51200000
# coefficients changed (95% upper bound on the rate, 3e-5).  So the constant
# matters for a KAT on intermediate values or on the expanded-key format, and
# not for a KAT on signatures.
#
# The contrast is the point.  This constant perturbs the sampler's WIDTH, which
# reaches `BerExp`'s comparison smoothly.  The three respellings of
# docs/debug_log.md #048 perturb its CENTRE, which passes through `floor(mu)`,
# and those DO change the signature -- 5 of the same 100000
# (scripts/divergence_rate.jl, docs/debug_log.md #054).  Whether an
# underdetermined choice reaches the signature depends on which quantity it
# reaches, not on how large it is.
#
# Everything below is computed at 400 bits and compared on raw bit patterns.

using Printf
setprecision(BigFloat, 400)

# ---------------------------------------------------------------------------
# Transcribed from scripts/cref/fpr.h, the FALCON_FPNATIVE branch.
# ---------------------------------------------------------------------------

"[C-ref] fpr.h, `fpr_inv_sigma[]`, indexed by logn (entry 1 is the unused slot)"
const INV_SIGMA_DECIMAL = (
    "0.0",
    "0.0069054793295940891952143765991630516",
    "0.0068102267767177975961393730687908629",
    "0.0067188101910722710707826117910434131",
    "0.0065883354370073665545865037227681924",
    "0.0064651781207602900738053897763485516",
    "0.0063486788828078995327741182928037856",
    "0.0062382586529084374473367528433697537",
    "0.0061334065020930261548984001431770281",
    "0.0060336696681577241031668062510953022",
    "0.0059386453095331159950250124336477482",
)

"[C-ref] fpr.h, `fpr_sigma_min[]`.  The header comment at inner.h:749 calls
this `1/sigma_min`; the values are `sigma_min` itself."
const SIGMA_MIN_DECIMAL = (
    "0.0",
    "1.1165085072329102588881898380334015",
    "1.1321247692325272405718031785357108",
    "1.1475285353733668684571123112513188",
    "1.1702540788534828939713084716509250",
    "1.1925466358390344011122170489094133",
    "1.2144300507766139921088487776957699",
    "1.2359260567719808790104525941706723",
    "1.2570545284063214162779743112075080",
    "1.2778336969128335860256340575729042",
    "1.2982803343442918539708792538826807",
)

"[Spec] Table 3.3, as printed (twelve significant figures)"
const SIGMA_PRINTED = Dict(9 => "165.736617183", 10 => "168.388571447")

"The Float64 literal the Python reference carries.  [Py-ref] falcon.py:131,139"
const SIGMA_F64 = Dict(9 => 165.7366171829776, 10 => 168.38857144654395)

bits(x::Float64) = string(reinterpret(UInt64, x); base = 16, pad = 16)
ulps(a::Float64, b::Float64) = Int(reinterpret(Int64, a) - reinterpret(Int64, b))

"""
Equation (2.13) of the specification, evaluated at `prec` bits.

    sigma = (1/pi) * sqrt( log(4n(1 + 1/eps)) / 2 ) * 1.17 * sqrt(q)

with `eps <= 1/sqrt(Q_s * lambda)`, `Q_s = 2^64`, `lambda = 128` (NIST Level I)
or `256` (Level V).  Note the `4n` inside the logarithm: this is the smoothing
parameter of `Z^{2n}`, not of `Z`.  The one-dimensional
`eta_eps(Z) = (1/pi) sqrt( (1/2) ln(2(1 + 1/eps)) )`, which is what reproduces
`sigma_min`, is a different formula with a different epsilon -- see
`smoothing_eta` in src/params.jl and docs/debug_log.md #046.  Conflating the
two is easy and gives a value 10% off.
"""
function spec_sigma(logn::Int)
    lambda = BigFloat(logn == 9 ? 128 : 256)
    n = BigFloat(2)^logn
    eps = 1 / sqrt(BigFloat(2)^64 * lambda)
    return (1 / BigFloat(pi)) * sqrt(log(4 * n * (1 + 1 / eps)) / 2) *
           BigFloat("1.17") * sqrt(BigFloat(12289))
end

function main()
    println("# What sigma does each fpr_inv_sigma entry imply?")
    println()
    @printf("%5s  %-18s  %-22s  %s\n", "logn", "table (bits)", "1 / table", "")
    for logn in 1:10
        t = parse(Float64, INV_SIGMA_DECIMAL[logn + 1])
        tb = parse(BigFloat, INV_SIGMA_DECIMAL[logn + 1])
        note = logn == 9 ? "Falcon-512" : logn == 10 ? "Falcon-1024" : ""
        @printf("%5d  %-18s  %-22.13f  %s\n", logn, bits(t), Float64(1 / tb), note)
    end

    println()
    println("# The two standardised sets, against every derivation we can construct")
    for logn in (9, 10)
        t = parse(Float64, INV_SIGMA_DECIMAL[logn + 1])
        se = spec_sigma(logn)
        candidates = (
            ("correctly rounded 1 / sigma_printed  ",
             Float64(1 / parse(BigFloat, SIGMA_PRINTED[logn]))),
            ("correctly rounded 1 / sigma_Float64  ",
             Float64(1 / BigFloat(SIGMA_F64[logn]))),
            ("Float64 division 1.0 / sigma_Float64 ", 1 / SIGMA_F64[logn]),
            ("correctly rounded 1 / sigma_(2.13)   ", Float64(1 / se)),
            ("(2.13) in Float64, then 1.0 / it     ", begin
                 epsf = 1 / sqrt(2.0^64 * (logn == 9 ? 128.0 : 256.0))
                 sf = (1 / pi) * sqrt(log(4 * 2.0^logn * (1 + 1 / epsf)) / 2) *
                      1.17 * sqrt(12289.0)
                 1 / sf
             end),
        )
        println()
        @printf("logn = %d      fpr_inv_sigma = %s\n", logn, bits(t))
        for (name, v) in candidates
            @printf("  %s %s  ulp %+d\n", name, bits(v), ulps(t, v))
        end
        # is the 35-digit literal the expansion of 1/sigma, or of a Float64?
        exact_inv = @sprintf("%.36f", 1 / se)
        litstr = INV_SIGMA_DECIMAL[logn + 1]
        agree = 0
        for i in 1:min(length(exact_inv), length(litstr))
            exact_inv[i] == litstr[i] || break
            agree += 1
        end
        @printf("  1/sigma at 400 bits    %s\n", exact_inv)
        @printf("  the literal in fpr.h   %s\n", litstr)
        @printf("  -> they agree to %d characters, then diverge\n", agree)
    end

    println()
    println("# inner.h:749 says fpr_sigma_min[] holds 1/sigma_min.  It holds sigma_min.")
    for logn in (9, 10)
        v = parse(Float64, SIGMA_MIN_DECIMAL[logn + 1])
        @printf("  logn=%2d   stored %.15f   Table 3.3 sigma_min %s   (1/stored = %.6f)\n",
                logn, v, logn == 9 ? "1.277833697" : "1.298280334", 1 / v)
    end
end

main()
