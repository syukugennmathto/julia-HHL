"""
    Falcon

A specification-conformant Julia implementation of FALCON / FN-DSA
(NIST FIPS 206, draft).

Scope and non-goals
-------------------
The target is FALCON-512 key generation, signing and verification passing the
Known Answer Tests.  **Constant-time execution and side-channel resistance are
explicitly out of scope.**  Where a production implementation would have to do
something different -- and in FALCON that is a long list, because its
floating-point Gaussian sampler is the reason FN-DSA was the last of the four
NIST PQC selections to reach standardisation -- the code carries a comment
saying so.  Those comments are the point, not an afterthought.

Module map (built in this order, each cross-checked against a reference
implementation before the next is started):

| module         | contents                                             |
|:---------------|:-----------------------------------------------------|
| `params.jl`    | parameter sets, with per-constant provenance         |
| `shake.jl`     | SHAKE256 XOF, ChaCha20 PRNG of the reference         |
| `poly.jl`      | schoolbook arithmetic in `Z[x]/(x^n+1)` and mod q    |
| `ntt.jl`       | NTT over `Z_q`, `q = 12289 = 3*2^12 + 1`             |
| `fft.jl`       | FFT over `C`, `splitfft`/`mergefft`, FFT-domain ops  |
| `ntrugen.jl`   | solving `f*G - g*F = q` by field-norm descent        |
| `samplerz.jl`  | discrete Gaussian sampler over `Z`                   |
| `ffsampling.jl`| LDL* decomposition, Falcon tree, fast Fourier sampling|
| `encoding.jl`  | key/signature serialisation, Golomb-Rice compression |
| `falcon.jl`    | `keygen` / `sign` / `verify`                         |

Modules below the horizontal rule are not yet written; see docs/debug_log.md
for the running state of the project.
"""
module Falcon

# --- module 1 --------------------------------------------------------------
export FalconParams, FALCON_512, FALCON_1024, params, Q, HEAD_LEN, SALT_LEN, SEED_LEN
export SIGMA_FG_BASE, smoothing_eta, falcon_eps, gram_schmidt_quality

# --- module 2 --------------------------------------------------------------
export SHAKE256XOF, shake256_xof, squeeze!, shake256
export ChaCha20, chacha20, randombytes!

# --- module 3 --------------------------------------------------------------
export polyadd, polysub, polyneg, polymul, polyadj, sqnorm
export polyaddq, polysubq, polymulq, centered
export polysplit, polymerge

include("params.jl")
include("shake.jl")
include("poly.jl")

end # module
