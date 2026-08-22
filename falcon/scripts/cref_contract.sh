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
echo "(step 3b re-signs a smaller batch per build; its streams are truncated to match)"


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
echo "## step 3b -- the negative control: do ORDINARY build choices diverge?"
M=${M:-30000}
build() { nm=$1; shift; "$@" -I "$HERE/cref" -o "$OUT/$nm" "$HERE/cref_contract_diff.c" $SRC -lm; }
build m_g_base  gcc   -O2 -DFALCON_FPNATIVE=1
build m_g_mnat  gcc   -O2 -march=native -DFALCON_FPNATIVE=1
build m_g_avx2  gcc   -O2 -march=native -DFALCON_FPNATIVE=1 -DFALCON_AVX2=1
build m_g_emu   gcc   -O2 -DFALCON_FPEMU=1
build m_g_emun  gcc   -O2 -march=native -DFALCON_FPEMU=1
build m_c_base  clang -O2 -DFALCON_FPNATIVE=1
build m_c_avx2  clang -O2 -march=native -DFALCON_FPNATIVE=1 -DFALCON_AVX2=1
build m_c_emu   clang -O2 -DFALCON_FPEMU=1
# and the two contracted builds of the hand-vectorised path, to see whether the
# AVX2 intrinsics are contracted too and whether they land on the same stream
build x_c_avx2f clang -O2 -march=native -ffp-contract=fast -DFALCON_AVX2=1 -DFALCON_FPNATIVE=1
build x_g_avx2f gcc   -O2 -march=native -ffp-contract=fast -DFALCON_AVX2=1 -DFALCON_FPNATIVE=1
for NAME in m_g_base m_g_mnat m_g_avx2 m_g_emu m_g_emun m_c_base m_c_avx2 m_c_emu \
            x_c_avx2f x_g_avx2f; do
  "$OUT/$NAME" "$M" > "$OUT/$NAME.bin" 2>/dev/null
done
python3 - "$OUT" "$M" <<'PY2'
import sys, os, itertools
O, M = sys.argv[1], int(sys.argv[2])
names = ["m_g_base","m_g_mnat","m_g_avx2","m_g_emu","m_g_emun",
         "m_c_base","m_c_avx2","m_c_emu","x_g_avx2f","cl_fast","x_c_avx2f"]
label = {"m_g_base":"gcc -O2 native", "m_g_mnat":"gcc -O2 -march=native",
         "m_g_avx2":"gcc AVX2 path", "m_g_emu":"gcc -O2 emulated FP",
         "m_g_emun":"gcc -march=native emulated", "m_c_base":"clang -O2 native",
         "m_c_avx2":"clang AVX2 path", "m_c_emu":"clang -O2 emulated FP",
         "x_g_avx2f":"gcc AVX2 -ffp-contract=fast",
         "cl_fast":"clang -ffp-contract=fast",
         "x_c_avx2f":"clang AVX2 -ffp-contract=fast"}
d = {}
for n in names:
    p = os.path.join(O, n + ".bin")
    if os.path.exists(p):
        d[n] = open(p, "rb").read()[:8*M]
n = min(len(v) for v in d.values()) // 8
groups = {}
for k in d: groups.setdefault(d[k][:8*n], []).append(k)
print("%d signatures per build" % n)
for i, (_, v) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
    print("  group %d (byte-identical): %s" % (i+1, ", ".join(label[x] for x in v)))
for a, b in itertools.combinations(list(d), 2):
    x, y = d[a][:8*n], d[b][:8*n]
    c = sum(1 for i in range(n) if x[8*i:8*i+8] != y[8*i:8*i+8])
    if c:
        print("  %-28s vs %-28s : %d of %d differ" % (label[a], label[b], c, n))
PY2

echo
echo "## step 4 -- key recovery from each divergent pair"
"$OUT/cl_fast" "$N" -d "$IDX" > "$OUT/dump_fast.txt" 2>/dev/null
"$OUT/cl_def"  "$N" -d "$IDX" > "$OUT/dump_def.txt"  2>/dev/null
cd "$ROOT"
julia --project=falcon falcon/scripts/cref_contract_recover.jl \
      "$OUT/dump_fast.txt" "$OUT/dump_def.txt"

echo
echo "## step 5 -- the same at FALCON-1024"
K=${K:-60000}
clang -O2 -march=native -ffp-contract=fast -DFALCON_FPNATIVE=1 -DFALCON_LOGN=10 \
      -I "$HERE/cref" -o "$OUT/f10_fast" "$HERE/cref_contract_diff.c" $SRC -lm
clang -O2 -march=native                    -DFALCON_FPNATIVE=1 -DFALCON_LOGN=10 \
      -I "$HERE/cref" -o "$OUT/f10_def"  "$HERE/cref_contract_diff.c" $SRC -lm
"$OUT/f10_fast" "$K" > "$OUT/f10_fast.bin" 2> "$OUT/f10_fast.err"
"$OUT/f10_def"  "$K" > "$OUT/f10_def.bin"  2> "$OUT/f10_def.err"
echo "fast $(head -1 "$OUT/f10_fast.err")   def $(head -1 "$OUT/f10_def.err")"
python3 - "$OUT/f10_fast.bin" "$OUT/f10_def.bin" <<'PY3'
import sys
a=open(sys.argv[1],'rb').read(); b=open(sys.argv[2],'rb').read()
n=min(len(a),len(b))//8
d=[i for i in range(n) if a[8*i:8*i+8]!=b[8*i:8*i+8]]
print("FALCON-1024: %d of %d differ (rate %.3g)" % (len(d), n, len(d)/n))
PY3
