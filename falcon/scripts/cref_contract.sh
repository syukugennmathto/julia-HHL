#!/bin/sh
# cref_contract.sh -- do two conforming builds of the FALCON reference
# implementation produce the same signatures?
#
#     sh falcon/scripts/cref_contract.sh [nsig]        (default 100000)
#
# Step 0 probes which compiler settings contract `a*b + c` into a fused
# multiply-add, on the reference's own `fpr` wrapper shape as well as on plain
# doubles.  Step 1 builds the reference three ways and counts FMA instructions.
# Step 2 signs `nsig` messages under one key with a per-message deterministic
# tape and hashes each signature.  Step 3 diffs the streams.  Step 4 hands every
# divergent pair to the key recovery of ePrint 2024/1709 section 5.1.
#
# Everything but the floating point is identical across builds by construction:
# same source, same key, same salt, same PRNG tape.
set -e
N=${1:-100000}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-/tmp/cref_contract}
mkdir -p "$OUT"
SRC="$HERE/cref/codec.c $HERE/cref/common.c $HERE/cref/falcon.c $HERE/cref/fft.c \
$HERE/cref/fpr.c $HERE/cref/keygen.c $HERE/cref/rng.c $HERE/cref/shake.c \
$HERE/cref/sign.c $HERE/cref/vrfy.c"

echo "## step 0 -- which settings contract, on which expression shape"
cat > "$OUT/probe.c" <<'PROBE'
typedef struct { double v; } fpr;                 /* the reference's wrapper */
static inline fpr FPR(double v){ fpr x; x.v=v; return x; }
static inline fpr fpr_mul(fpr a, fpr b){ return FPR(a.v*b.v); }
static inline fpr fpr_sub(fpr a, fpr b){ return FPR(a.v-b.v); }
double plain(double a,double b,double c,double d){ return a*b - c*d; }
double wrapped(double a,double b,double c,double d){
  return fpr_sub(fpr_mul(FPR(a),FPR(b)), fpr_mul(FPR(c),FPR(d))).v;
}
PROBE
printf '%-34s %7s %7s\n' "compiler and flags" "plain" "wrapped"
for CC in gcc clang; do
  for F in "-O2" "-O2 -ffp-contract=on" "-O2 -ffp-contract=fast" "-O2 -ffp-contract=off" "-Ofast"; do
    $CC $F -march=native -c "$OUT/probe.c" -o "$OUT/probe.o" 2>/dev/null || continue
    P=$(objdump -d "$OUT/probe.o" | awk '/<plain>:/,/^$/'   | grep -c 'fmadd\|fmsub' || true)
    W=$(objdump -d "$OUT/probe.o" | awk '/<wrapped>:/,/^$/' | grep -c 'fmadd\|fmsub' || true)
    printf '%-34s %7s %7s\n' "$CC $F" "$P" "$W"
  done
done

echo
echo "## step 1 -- three builds of the reference, native double"
for B in "cl_fast:clang -O2 -march=native -ffp-contract=fast" \
         "cl_def:clang -O2 -march=native" \
         "gc_def:gcc -O2 -march=native"; do
  NAME=${B%%:*}; CMD=${B#*:}
  $CMD -DFALCON_FPNATIVE=1 -I "$HERE/cref" -o "$OUT/$NAME" \
       "$HERE/cref_contract_diff.c" $SRC -lm
  echo "$NAME: $(objdump -d "$OUT/$NAME" | grep -c 'fmadd\|fmsub\|fnmadd\|fnmsub' || true) fma instructions   [$CMD]"
done

echo
echo "## step 2 -- $N signatures per build, one key, deterministic tape"
for NAME in cl_fast cl_def gc_def; do
  "$OUT/$NAME" "$N" > "$OUT/$NAME.bin" 2> "$OUT/$NAME.err"
  echo "$NAME  $(head -1 "$OUT/$NAME.err")"
done

echo
echo "## step 3 -- do the streams agree?"
if cmp -s "$OUT/cl_def.bin" "$OUT/gc_def.bin"; then
  echo "gcc default == clang default : identical over $N signatures"
else
  echo "gcc default != clang default"
fi
IDX=$(python3 - "$OUT/cl_fast.bin" "$OUT/cl_def.bin" <<'PY'
import sys
a=open(sys.argv[1],'rb').read(); b=open(sys.argv[2],'rb').read()
n=min(len(a),len(b))//8
d=[i for i in range(n) if a[8*i:8*i+8]!=b[8*i:8*i+8]]
print("clang -ffp-contract=fast vs clang default : %d of %d differ (rate %.3g)"
      % (len(d), n, len(d)/n), file=sys.stderr)
print(",".join(str(i) for i in d))
PY
)
[ -n "$IDX" ] || { echo "no divergence in $N signatures; try more"; exit 0; }

echo
echo "## step 4 -- key recovery from each divergent pair"
"$OUT/cl_fast" "$N" -d "$IDX" > "$OUT/dump_fast.txt" 2>/dev/null
"$OUT/cl_def"  "$N" -d "$IDX" > "$OUT/dump_def.txt"  2>/dev/null
cd "$ROOT"
julia --project=falcon falcon/scripts/cref_contract_recover.jl \
      "$OUT/dump_fast.txt" "$OUT/dump_def.txt"
