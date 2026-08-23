#!/bin/sh
# Build the C reference under several compilers and optimisation levels and
# time all of them, so the comparison is against a RANGE of C rather than one
# build.  See docs/benchmarks.md, "Methods".
#
#     sh falcon/scripts/bench_matrix.sh [outdir] [iters_keygen] [iters_sign] [iters_verify]
#
# The official config.h leaves FALCON_FPEMU and FALCON_FPNATIVE both commented
# out, which selects native double (docs/debug_log.md #045).  Passing
# -DFALCON_FPEMU=1 selects the emulated build -- the one Algorand ships, and
# therefore the one most deployed FALCON code runs.
set -e
OUT=${1:-/tmp/benchmatrix}
KG=${2:-40}; SG=${3:-200}; VF=${4:-400}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT" && cd "$OUT"
cp "$HERE"/cref/*.c "$HERE"/cref/*.h . && cp "$HERE"/cref_bench.c .
CFILES="cref_bench.c codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c sign.c vrfy.c"

for CC in gcc clang; do
  command -v "$CC" >/dev/null 2>&1 || continue
  for OPT in "-O2" "-O3" "-O3 -march=native"; do
    tag="$CC$(echo "$OPT" | tr -d ' -')"
    $CC $OPT -o "b_nat_$tag" $CFILES -lm
    $CC $OPT -DFALCON_FPEMU=1 -o "b_emu_$tag" $CFILES -lm
  done
done

for b in b_*; do
  printf '### %s\n' "$b"
  "./$b" 9 "$KG" "$SG" "$VF"
done
