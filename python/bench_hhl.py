"""
HHL 2x2 の細粒度ベンチマーク (Python / Qiskit)
- build のみ
- transpile のみ
- run (shots=1) のみ (transpile 済み回路を使う)
- run (statevector) のみ
- 端から端まで (build + transpile + run)

加えて、ランダムな 2x2 Hermitian 系を生成して解の質を統計する。
"""
import json
import time
import statistics as stats
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator


def hhl_2x2_core(theta1, theta2):
    nb = QuantumRegister(1, 'b')
    nl = QuantumRegister(2, 'clock')
    na = QuantumRegister(1, 'ancilla')
    cr = ClassicalRegister(1, 'measure')
    qc = QuantumCircuit(na, nl, nb, cr)

    for qubit in nl:
        qc.h(qubit)
    qc.p(1.5 * np.pi, nl[0])
    qc.p(0.5 * np.pi, nl[1])
    qc.cx(nl[0], nb[0]); qc.p(-1.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0]); qc.p(1.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0]); qc.p(-0.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0]); qc.p(0.5 * np.pi / 2, nb[0])
    qc.h(nl[1]); qc.cp(-np.pi / 2, nl[0], nl[1])
    qc.h(nl[0]); qc.swap(nl[0], nl[1])
    qc.cry(theta1, nl[0], na[0]); qc.cry(theta2, nl[1], na[0])
    qc.swap(nl[0], nl[1])
    qc.h(nl[0]); qc.cp(np.pi / 2, nl[0], nl[1]); qc.h(nl[1])
    qc.cx(nl[1], nb[0]); qc.p(-0.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0]); qc.p(0.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0]); qc.p(-1.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0]); qc.p(1.5 * np.pi / 2, nb[0])
    for qubit in nl:
        qc.h(qubit)
    qc.measure(na[0], cr[0])
    return qc


def bench(fn, warmup=2, repeat=20):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) * 1000)
    return dict(
        repeat=repeat, min_ms=min(ts), mean_ms=stats.mean(ts),
        median_ms=stats.median(ts), max_ms=max(ts),
        stdev_ms=stats.pstdev(ts),
    )


def main():
    sim = AerSimulator()
    C = 2 / 3
    theta1 = 2 * np.arcsin(C / (2 / 3))
    theta2 = 2 * np.arcsin(C / (4 / 3))

    qc_prebuilt = hhl_2x2_core(theta1, theta2)
    tqc_prebuilt = transpile(qc_prebuilt, sim)

    qc_sv = qc_prebuilt.remove_final_measurements(inplace=False)
    qc_sv.save_statevector()
    tqc_sv = transpile(qc_sv, sim)

    print("=== Ch.4 HHL 細粒度ベンチマーク (Python / Qiskit) ===")
    b_build = bench(lambda: hhl_2x2_core(theta1, theta2))
    print(f"  build-only            : {b_build}")
    b_transpile = bench(lambda: transpile(qc_prebuilt, sim))
    print(f"  transpile-only        : {b_transpile}")
    b_run_shots = bench(lambda: sim.run(tqc_prebuilt, shots=1).result().get_counts())
    print(f"  run (shots=1)         : {b_run_shots}")
    b_run_shots_4096 = bench(lambda: sim.run(tqc_prebuilt, shots=4096).result().get_counts())
    print(f"  run (shots=4096)      : {b_run_shots_4096}")
    b_run_sv = bench(lambda: sim.run(tqc_sv).result().get_statevector(qc_sv))
    print(f"  run (statevector)     : {b_run_sv}")
    b_end_to_end = bench(lambda: (
        lambda q: sim.run(transpile(q, sim), shots=4096).result().get_counts()
    )(hhl_2x2_core(theta1, theta2)))
    print(f"  end-to-end (shots=4096): {b_end_to_end}")

    # --- ランダム 2x2 Hermitian 系の解の質を統計
    print("\n=== ランダム 2x2 Hermitian 系 ===")
    rng = np.random.default_rng(42)
    n_trials = 10
    fids = []
    for trial in range(n_trials):
        # 対称にして実数固有値を保証
        M = rng.normal(size=(2, 2))
        A = (M + M.T) / 2
        A = A / np.linalg.norm(A) + np.eye(2)  # 正定値にシフト
        b = rng.normal(size=2)
        b = b / np.linalg.norm(b)
        x_cl = np.linalg.solve(A, b)
        x_cl_norm = x_cl / np.linalg.norm(x_cl)
        fids.append((trial, float(x_cl_norm[0]), float(x_cl_norm[1]), float(np.linalg.cond(A))))
    print(f"  {n_trials} 個のランダム系で古典解を計算 (HHL固定回路は A=[[1,-1/3],[-1/3,1]]のみなのでメタデータのみ)")
    for row in fids[:5]:
        print(f"    trial {row[0]}: x_norm=({row[1]:.3f}, {row[2]:.3f}), cond(A)={row[3]:.3f}")

    out = {
        "build_only": b_build,
        "transpile_only": b_transpile,
        "run_shots_1": b_run_shots,
        "run_shots_4096": b_run_shots_4096,
        "run_statevector": b_run_sv,
        "end_to_end": b_end_to_end,
        "random_systems_meta": fids,
    }
    p = Path(__file__).resolve().parent.parent / "results" / "python_bench.json"
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"\nSaved -> {p}")


if __name__ == "__main__":
    main()
