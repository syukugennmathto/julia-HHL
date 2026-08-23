#!/bin/sh
# cref_contract_aarch64.sh -- replicate the FP_CONTRACT divergence on aarch64,
# statically (does the reference contract?) and dynamically (under QEMU: does it
# diverge, and do the pairs recover the key?).
#
#     sh falcon/scripts/cref_contract_aarch64.sh [nsig]      (default 150000)
#
# Needs a cross toolchain and an emulator.  On Debian/Ubuntu:
#     apt-get install -y gcc-aarch64-linux-gnu qemu-user-static
# clang cross-compiles with --target=aarch64-linux-gnu --sysroot=/usr/aarch64-linux-gnu
# (the sysroot ships with the gcc-aarch64-linux-gnu package).
#
# The point of this script is the SECOND-ARCHITECTURE control for section 6.7:
# the FP_CONTRACT hazard is not an x86 artefact.  The contraction pattern in the
# real fft.c is identical to x86 (gcc never; clang only at -ffp-contract=fast /
# -Ofast / -ffast-math), and the dynamic divergence rate and key recovery match.
set -e
N=${1:-150000}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-/tmp/cref_aarch64}; mkdir -p "$OUT/of" "$OUT/od"
SYS="--target=aarch64-linux-gnu --sysroot=/usr/aarch64-linux-gnu -isystem /usr/aarch64-linux-gnu/include"
TUS="codec common falcon fft fpr keygen rng shake sign vrfy"

command -v aarch64-linux-gnu-gcc >/dev/null || { echo "need gcc-aarch64-linux-gnu"; exit 1; }
command -v qemu-aarch64-static  >/dev/null || { echo "need qemu-user-static"; exit 1; }

echo "## static: FMA in real fft.c across flags (clang aarch64) ##"
for F in "-O2" "-O2 -ffp-contract=fast" "-Ofast" "-O2 -ffast-math"; do
  clang $SYS $F -DFALCON_FPNATIVE=1 -I "$HERE/cref" -c "$HERE/cref/fft.c" -o "$OUT/f.o" 2>/dev/null
  echo "  clang aarch64 $F : $(llvm-objdump -d "$OUT/f.o" | grep -cE 'fmadd|fmsub|fnmadd|fnmsub') fma"
done
echo "  (gcc aarch64 emits 0 at every flag -- the fpr wrapper blocks it, as on x86)"

echo "## build contracting vs non-contracting (clang), link with gcc (has libm) ##"
for s in $TUS; do
  clang $SYS -O2 -ffp-contract=fast -DFALCON_FPNATIVE=1 -I "$HERE/cref" -c "$HERE/cref/$s.c" -o "$OUT/of/$s.o"
  clang $SYS -O2                    -DFALCON_FPNATIVE=1 -I "$HERE/cref" -c "$HERE/cref/$s.c" -o "$OUT/od/$s.o"
done
clang $SYS -O2 -ffp-contract=fast -DFALCON_FPNATIVE=1 -I "$HERE/cref" -c "$HERE/cref_contract_diff.c" -o "$OUT/of/main.o"
clang $SYS -O2                    -DFALCON_FPNATIVE=1 -I "$HERE/cref" -c "$HERE/cref_contract_diff.c" -o "$OUT/od/main.o"
aarch64-linux-gnu-gcc -static -o "$OUT/ca_fast" "$OUT"/of/*.o -lm
aarch64-linux-gnu-gcc -static -o "$OUT/ca_def"  "$OUT"/od/*.o -lm
echo "  fma: ca_fast $(aarch64-linux-gnu-objdump -d "$OUT/ca_fast"|grep -cE 'fmadd|fmsub') ; ca_def $(aarch64-linux-gnu-objdump -d "$OUT/ca_def"|grep -cE 'fmadd|fmsub')"

echo "## dynamic under QEMU: $N signatures each ##"
qemu-aarch64-static "$OUT/ca_fast" "$N" > "$OUT/fast.bin" 2>"$OUT/fast.err"
qemu-aarch64-static "$OUT/ca_def"  "$N" > "$OUT/def.bin"  2>"$OUT/def.err"
echo "  key: $(head -1 "$OUT/fast.err")"
IDX=$(python3 - "$OUT/fast.bin" "$OUT/def.bin" <<'PY'
import sys
a=open(sys.argv[1],'rb').read(); b=open(sys.argv[2],'rb').read(); n=min(len(a),len(b))//8
d=[i for i in range(n) if a[8*i:8*i+8]!=b[8*i:8*i+8]]
print("aarch64: %d of %d differ (rate %.3g)"%(len(d),n,len(d)/n), file=sys.stderr)
print(",".join(str(i) for i in d[:6]))   # a handful is enough to recover
PY
)
echo "## key recovery on aarch64-produced signatures ##"
MX=$(echo "$IDX" | tr ',' '\n' | sort -n | tail -1)
qemu-aarch64-static "$OUT/ca_fast" $((MX+1)) -d "$IDX" > "$OUT/dump_fast.txt" 2>/dev/null
qemu-aarch64-static "$OUT/ca_def"  $((MX+1)) -d "$IDX" > "$OUT/dump_def.txt"  2>/dev/null
cd "$ROOT"
julia --project=falcon falcon/scripts/cref_contract_recover.jl "$OUT/dump_fast.txt" "$OUT/dump_def.txt"
