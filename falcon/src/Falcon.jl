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

All ten modules are implemented and their tests pass.  See docs/debug_log.md
for the project's running state, and docs/math/ for the mathematical notes.
"""
module Falcon

# --- module 1 --------------------------------------------------------------
export FalconParams, FALCON_512, FALCON_1024, params, Q, HEAD_LEN, SALT_LEN, SEED_LEN
export SIGMA_FG_BASE, smoothing_eta, falcon_eps, gram_schmidt_quality

# --- module 2 --------------------------------------------------------------
export SHAKE256XOF, shake256_xof, absorb!, squeeze!, shake256
export ChaCha20, chacha20, randombytes!

# --- module 3 --------------------------------------------------------------
export polyadd, polysub, polyneg, polymul, polyadj, sqnorm
export polyaddq, polysubq, polymulq, centered
export polysplit, polymerge

# --- module 4 --------------------------------------------------------------
export ntt, intt, ntt_roots, ntt_add, ntt_sub, ntt_mul, ntt_div
export polymulq_ntt, polydivq, is_invertible_zq
export ntt_ip!, intt_ip!, ntt_zetas, intt_zetas, polymulq_fast, polymulq_fast!

# --- module 5 --------------------------------------------------------------
export fft, ifft, fft_roots, unit_root, split_fft, merge_fft
export add_fft, sub_fft, neg_fft, mul_fft, div_fft, adj_fft
export polymul_fft, polydiv_fft, polyadj_fft
export set_fft_roots!, reset_fft_roots!, with_fft_roots

# --- module 6 --------------------------------------------------------------
export karamul, galois_conjugate, field_norm, lift, bitsize
export ntru_solve, NTRUSolveFailure, ntru_equation_residual, ntru_equation_holds
export gs_norm, gs_norm_ok, gen_poly, ntru_gen, SIGMA_FG_MIN
export mkgauss, mkgauss_u64, gen_poly_cdt, GAUSS_1024_12289   # branch variant, see README

# --- module 7 --------------------------------------------------------------
export samplerz, basesampler, approxexp, berexp
export RCDT, EXP_COEFFS, bytesource, ReplayBytes

# --- module 8 --------------------------------------------------------------
export FalconTree, FFLDLNode, FFLDLLeaf
export gram_fft, ldl_fft, ffldl_fft, normalize_tree!, falcon_tree
export ffnp_fft, ffsampling_fft, leaf_sigmas, node_l10s, nleaves, treedepth

# --- module 9 --------------------------------------------------------------
export encode_pubkey, decode_pubkey, encode_privkey, decode_privkey, recover_G
export compress_sig, decompress_sig, encode_signature, decode_signature
export MAX_FG_BITS, MAX_FG_BITS_F

# --- module 10 -------------------------------------------------------------
export FalconPublicKey, FalconPrivateKey
export hash_to_point, falcon_keygen, falcon_sign, falcon_verify
export expand_privkey, public_key, sample_preimage, signature_norm
export privkey_from_bytes, pubkey_from_bytes, pubkey_bytes, privkey_bytes

include("params.jl")
include("shake.jl")
include("poly.jl")
include("ntt.jl")
include("fft.jl")
include("ntrugen.jl")
include("samplerz.jl")
include("ffsampling.jl")
include("encoding.jl")
include("falcon.jl")

end # module
