#!/bin/sh
# falcon_fma.sh -- measure the reference's OWN FALCON_FMA option: does enabling
# it change signatures, at what rate, and does the change recover the key?
#
#     sh falcon/scripts/falcon_fma.sh [nsig]     (default 300000)
#
# Needs an x86-64 host with AVX2+FMA.  config.h documents FALCON_FMA and claims
# it changes signatures "with low probability, less than 2^(-40); produced
# signatures are still safe and interoperable."  This measures the real rate and
# runs the key recovery, which is what makes "safe and interoperable" false in
# the deterministic / cross-implementation setting.
set -e
N=${1:-300000}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-/tmp/falcon_fma}; mkdir -p "$OUT"
SRC="$HERE/cref/codec.c $HERE/cref/common.c $HERE/cref/falcon.c $HERE/cref/fft.c \
$HERE/cref/fpr.c $HERE/cref/keygen.c $HERE/cref/rng.c $HERE/cref/shake.c \
$HERE/cref/sign.c $HERE/cref/vrfy.c"
gcc -O2 -march=native -DFALCON_FPNATIVE=1 -DFALCON_AVX2=1 -DFALCON_FMA=1 -I "$HERE/cref" -o "$OUT/fma"   "$HERE/cref_contract_diff.c" $SRC -lm
gcc -O2 -march=native -DFALCON_FPNATIVE=1 -DFALCON_AVX2=1                 -I "$HERE/cref" -o "$OUT/nofma" "$HERE/cref_contract_diff.c" $SRC -lm
echo "fma insns: fma=$(objdump -d "$OUT/fma"|grep -cE 'vfmadd|vfmsub') default=$(objdump -d "$OUT/nofma"|grep -cE 'vfmadd|vfmsub')"
"$OUT/fma" "$N" > "$OUT/fma.bin" 2>"$OUT/fma.err"; "$OUT/nofma" "$N" > "$OUT/nofma.bin" 2>"$OUT/nofma.err"
IDX=$(python3 - "$OUT/fma.bin" "$OUT/nofma.bin" <<'PY'
import sys,math
a=open(sys.argv[1],'rb').read(); b=open(sys.argv[2],'rb').read(); n=min(len(a),len(b))//8
d=[i for i in range(n) if a[8*i:8*i+8]!=b[8*i:8*i+8]]; r=len(d)/n
print("FALCON_FMA=1 vs default: %d of %d differ (rate %.3g = 2^%.1f); config.h claims < 2^-40"%(len(d),n,r,math.log2(r)), file=sys.stderr)
print(",".join(str(i) for i in d[:6]))
PY
)
MX=$(echo "$IDX"|tr ',' '\n'|sort -n|tail -1)
"$OUT/fma" $((MX+1)) -d "$IDX" > "$OUT/dump_fma.txt" 2>/dev/null
"$OUT/nofma" $((MX+1)) -d "$IDX" > "$OUT/dump_nofma.txt" 2>/dev/null
cd "$ROOT"; julia --project=falcon falcon/scripts/cref_contract_recover.jl "$OUT/dump_fma.txt" "$OUT/dump_nofma.txt"
