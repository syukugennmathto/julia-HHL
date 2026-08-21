#!/bin/sh
# bench_all.sh -- the full three-way, multi-build performance matrix in one
# run, so the paper's performance section rests on a range of C builds and
# both Julia and Python rather than one number each.
#
#     sh falcon/scripts/bench_all.sh [outdir]
#
# Produces, in OUT:
#   c_<cc><opt>_<emu|nat>_<n>.txt   C reference, per compiler/opt/backend/degree
#   jl_O<k>_<n>.txt                 this project's Julia, per -O level/degree
#   py_<n>.txt                      the Python reference
#   MANIFEST                        what ran, with versions
#
# Every file is in the shared "op n median mean min p25 p75 iters" format.
#
# One machine only.  We cannot vary hardware here, so the spread we CAN show
# is across compilers (gcc, clang), optimisation levels (-O2, -O3,
# -O3 -march=native) and floating-point back ends (emulated vs native) -- and
# the run-to-run distribution via the quartile columns.  The single-machine
# limitation is stated in docs/benchmarks.md and the paper; it is not hidden.
set -e
OUT=${1:-/tmp/bench_all}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
mkdir -p "$OUT"

# iteration counts: generous for the fast back ends, small for slow ones
KG_C=40;  SG_C=400; VF_C=800
KG_JL=8;  SG_JL=200; VF_JL=400
KG_PY=3;  SG_PY=20;  VF_PY=60

{
  echo "date_utc: (Date.now unavailable in-harness; stamp externally)"
  echo "uname: $(uname -a)"
  echo "nproc: $(nproc 2>/dev/null || echo '?')"
  echo "gcc: $(gcc --version 2>/dev/null | head -1)"
  echo "clang: $(clang --version 2>/dev/null | head -1)"
  echo "julia: $(julia --version 2>/dev/null)"
  echo "python: $(python3 --version 2>/dev/null)"
  grep -m1 'model name' /proc/cpuinfo 2>/dev/null || true
} > "$OUT/MANIFEST"

# ---- C: build the matrix once, reuse for both degrees ----------------------
CB="$OUT/cbuild"
mkdir -p "$CB"
cp "$HERE"/cref/*.c "$HERE"/cref/*.h "$CB"/
cp "$HERE"/cref_bench.c "$CB"/
CFILES="cref_bench.c codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c sign.c vrfy.c"
( cd "$CB"
  for CC in gcc clang; do
    command -v "$CC" >/dev/null 2>&1 || continue
    for OPT in "-O2" "-O3" "-O3 -march=native"; do
      tag="$CC$(echo "$OPT" | tr -d ' -')"
      $CC $OPT              -o "b_nat_$tag" $CFILES -lm
      $CC $OPT -DFALCON_FPEMU=1 -o "b_emu_$tag" $CFILES -lm
    done
  done )

for LOGN in 9 10; do
  N=$((1 << LOGN))
  for b in "$CB"/b_*; do
    [ -x "$b" ] || continue
    base=$(basename "$b")            # e.g. b_emu_gccO3marchnative
    kind=${base#b_}                  # emu_gccO3marchnative
    echo "# $base  logn=$LOGN"
    "$b" "$LOGN" "$KG_C" "$SG_C" "$VF_C" > "$OUT/c_${kind}_${N}.txt" 2>&1 || \
      echo "  (failed: $base)"
  done

  # ---- Julia at three optimisation levels -------------------------------
  for O in 1 2 3; do
    julia -O$O --project="$ROOT/falcon" "$HERE/bench.jl" "$LOGN" "$KG_JL" "$SG_JL" "$VF_JL" \
      > "$OUT/jl_O${O}_${N}.txt" 2>&1 || echo "  (julia -O$O logn=$LOGN failed)"
  done

  # ---- Python reference -------------------------------------------------
  PYTHONPATH="$HERE/pyref" python3 "$HERE/pyref_bench.py" "$LOGN" "$KG_PY" "$SG_PY" "$VF_PY" \
    > "$OUT/py_${N}.txt" 2>&1 || echo "  (python logn=$LOGN failed)"
done

echo "done -> $OUT"
