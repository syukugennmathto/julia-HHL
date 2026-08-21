#!/usr/bin/env python3
#
# pyref_bench.py -- time the Python reference implementation's keygen / sign /
# verify, in the same shape as scripts/cref_bench.c and scripts/bench.jl, so
# all three can be put side by side.
#
#     PYTHONPATH=falcon/scripts/pyref python3 falcon/scripts/pyref_bench.py \
#         [logn] [keygen] [sign] [verify]
#
# Output: one "op n median_ms mean_ms min_ms p25_ms p75_ms iters" line per
# operation, plus a header, so scripts/bench_compare.jl can read it next to
# the C and Julia numbers.
#
# The Python reference (Prest's `falcon.py`) exists to be read, not to be
# fast: it carries every ring element as a Python list of Python floats /
# ints, and `ntru_gen` in particular is very slow.  Iteration counts are
# therefore small by default, and keygen especially so.  A large ratio to C
# is the whole point of showing it -- it is the "before" against which both C
# and this project's Julia are the "after".
#
# Fairness notes, matching the C and Julia drivers:
#   - the first call is not special in Python (no JIT), so no warm-up is done;
#   - `median` is quoted, not `mean`, because keygen retries and its
#     distribution has a long right tail;
#   - verify is timed from the serialized public key, decoding on every call,
#     exactly as C's falcon_verify and our bench.jl do.

import sys
import time
from statistics import median, mean

import falcon


def pct(xs, p):
    xs = sorted(xs)
    if not xs:
        return 0.0
    k = (len(xs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def report(op, n, t):
    print("%s %d %.6f %.6f %.6f %.6f %.6f %d" % (
        op, n, median(t), mean(t), min(t), pct(t, 0.25), pct(t, 0.75), len(t)))
    sys.stdout.flush()


def main():
    logn = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    nkg = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    nsg = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    nvf = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    n = 1 << logn

    print("# pyref_bench python=%s n=%d" % (sys.version.split()[0], n))
    print("# op n median_ms mean_ms min_ms p25_ms p75_ms iters")

    scheme = falcon.Falcon(n)
    msg = b"falcon-jl benchmark message"

    # keygen
    t = []
    sk = vk = None
    for _ in range(nkg):
        t0 = time.perf_counter()
        sk, vk = scheme.keygen()
        t.append((time.perf_counter() - t0) * 1e3)
    report("keygen", n, t)

    # sign (the Python reference rebuilds nothing per call: the tree lives in sk,
    # so this is the counterpart of C's sign_tree and our sign_tree)
    t = []
    sig = None
    for _ in range(nsg):
        t0 = time.perf_counter()
        sig = scheme.sign(sk, msg)
        t.append((time.perf_counter() - t0) * 1e3)
    report("sign_tree", n, t)

    # verify, from the serialized public key
    assert scheme.verify(vk, msg, sig)
    t = []
    for _ in range(nvf):
        t0 = time.perf_counter()
        scheme.verify(vk, msg, sig)
        t.append((time.perf_counter() - t0) * 1e3)
    report("verify", n, t)


if __name__ == "__main__":
    main()
