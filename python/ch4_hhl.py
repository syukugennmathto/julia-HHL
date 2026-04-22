"""第4章: HHL本体 (2x2) - Qiskit 実装 + ノイズシミュレーション + ベンチマーク"""
import json
import time
import statistics as stats
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error


def hhl_2x2():
    """
    A = [[1, -1/3], [-1/3, 1]], b = [1, 0]
    固有値 lambda1 = 2/3, lambda2 = 4/3
    古典解 x = [9/8, 3/8] (正規化前)
    """
    nb = QuantumRegister(1, 'b')
    nl = QuantumRegister(2, 'clock')
    na = QuantumRegister(1, 'ancilla')
    cr = ClassicalRegister(1, 'measure')
    qc = QuantumCircuit(na, nl, nb, cr)

    # QPE
    for qubit in nl:
        qc.h(qubit)

    qc.p(1.5 * np.pi, nl[0])
    qc.p(0.5 * np.pi, nl[1])
    qc.cx(nl[0], nb[0])
    qc.p(-1.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0])
    qc.p(1.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0])
    qc.p(-0.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0])
    qc.p(0.5 * np.pi / 2, nb[0])

    # QFT^{-1}
    qc.h(nl[1])
    qc.cp(-np.pi / 2, nl[0], nl[1])
    qc.h(nl[0])
    qc.swap(nl[0], nl[1])

    # 制御回転
    C = 2 / 3
    theta_1 = 2 * np.arcsin(C / (2 / 3))
    theta_2 = 2 * np.arcsin(C / (4 / 3))
    qc.cry(theta_1, nl[0], na[0])
    qc.cry(theta_2, nl[1], na[0])

    # 逆QPE
    qc.swap(nl[0], nl[1])
    qc.h(nl[0])
    qc.cp(np.pi / 2, nl[0], nl[1])
    qc.h(nl[1])

    qc.cx(nl[1], nb[0])
    qc.p(-0.5 * np.pi / 2, nb[0])
    qc.cx(nl[1], nb[0])
    qc.p(0.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0])
    qc.p(-1.5 * np.pi / 2, nb[0])
    qc.cx(nl[0], nb[0])
    qc.p(1.5 * np.pi / 2, nb[0])

    for qubit in nl:
        qc.h(qubit)

    qc.measure(na[0], cr[0])
    return qc, nb, nl, na, cr


def hhl_2x2_with_save(save_statevector=True):
    """State vector を取り出すため、測定を省いて save_statevector を付ける。"""
    qc, nb, nl, na, cr = hhl_2x2()
    qc2 = qc.remove_final_measurements(inplace=False)
    if save_statevector:
        qc2.save_statevector()
    return qc2


def extract_solution_from_statevector(sv_array, n_b=1, n_l=2, n_a=1):
    """
    Qiskit のレジスタ順は [ancilla, clock, b]。状態ベクトルのインデックス i を
    ビット列にした時 (Qiskit は i の LSB がレジスタ宣言順の最上位?)
    実際の並びは: bit0..bit_{na-1} = ancilla, bit_{na..}..= clock, bit_{...} = b
    -> idx = b_bit << (na+nl) | clock_bits << na | ancilla_bit
    アンシラ=1, クロック=00 の部分から b の振幅を取り出す。
    """
    n = n_a + n_l + n_b
    sol = np.zeros(2 ** n_b, dtype=complex)
    for b_val in range(2 ** n_b):
        # clock=0, ancilla=1
        idx = (b_val << (n_a + n_l)) | (0 << n_a) | 1
        sol[b_val] = sv_array[idx]
    norm = np.linalg.norm(sol)
    if norm < 1e-12:
        return sol
    return sol / norm


def benchmark(func, *, warmup=1, repeat=5):
    for _ in range(warmup):
        func()
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        func()
        times.append((time.perf_counter() - t0) * 1000)
    return {
        "repeat": repeat,
        "min_ms": min(times),
        "mean_ms": stats.mean(times),
        "median_ms": stats.median(times),
        "max_ms": max(times),
        "stdev_ms": stats.pstdev(times),
    }


def main():
    results = {}

    # --- 1. 回路の描画
    qc_meas, *_ = hhl_2x2()
    print("=== HHL 2x2 回路 (ゲート数) ===")
    print(f"  depth = {qc_meas.depth()}, size = {qc_meas.size()}")

    # --- 2. 測定ベース実行
    sim = AerSimulator()
    tqc = transpile(qc_meas, sim)
    job = sim.run(tqc, shots=4096)
    counts = job.result().get_counts()
    print(f"\n[測定結果 shots=4096] {counts}")
    p_success = counts.get('1', 0) / 4096
    print(f"  後選択成功確率 P(ancilla=1) ≈ {p_success:.4f}")
    results["ancilla_success_prob"] = p_success

    # --- 3. 状態ベクトルから解を取り出す
    qc_sv = hhl_2x2_with_save(save_statevector=True)
    sv_job = sim.run(transpile(qc_sv, sim))
    sv = sv_job.result().get_statevector(qc_sv)
    sv_array = np.asarray(sv)
    sol = extract_solution_from_statevector(sv_array)
    print(f"\n[量子解(state vector から)] |x> ∝ {np.abs(sol)}")

    # 古典解
    A = np.array([[1.0, -1.0 / 3], [-1.0 / 3, 1.0]])
    b = np.array([1.0, 0.0])
    x_cl = np.linalg.solve(A, b)
    x_cl_norm = x_cl / np.linalg.norm(x_cl)
    print(f"[古典解] x = {x_cl}, 正規化後 = {x_cl_norm}")

    fidelity = abs(np.vdot(sol, x_cl_norm.astype(complex))) ** 2
    print(f"[Fidelity |<x_q|x_cl>|^2] = {fidelity:.6f}")
    results["fidelity"] = fidelity
    results["classical_solution"] = x_cl.tolist()
    results["quantum_solution_abs"] = np.abs(sol).tolist()

    # --- 4. ノイズシミュレーション
    print("\n=== ノイズあり ===")
    noise_model = NoiseModel()
    err1 = depolarizing_error(0.001, 1)
    err2 = depolarizing_error(0.01, 2)
    noise_model.add_all_qubit_quantum_error(err1, ['u1', 'u2', 'u3', 'p', 'h', 'x', 'ry', 'rz', 'sx'])
    noise_model.add_all_qubit_quantum_error(err2, ['cx', 'cp', 'cry', 'swap'])

    noisy_sim = AerSimulator(noise_model=noise_model)
    tqc_noisy = transpile(qc_meas, noisy_sim)
    job_n = noisy_sim.run(tqc_noisy, shots=4096)
    noisy_counts = job_n.result().get_counts()
    p_succ_n = noisy_counts.get('1', 0) / 4096
    print(f"  ノイズあり counts={noisy_counts}, P(anc=1)={p_succ_n:.4f}")
    results["noisy_ancilla_success_prob"] = p_succ_n

    # --- 5. ベンチマーク (shots=4096, 測定あり版の build+transpile+run)
    print("\n=== ベンチマーク (Python / Qiskit) ===")

    def run_measure_4096():
        qc, *_ = hhl_2x2()
        tqc = transpile(qc, sim)
        sim.run(tqc, shots=4096).result().get_counts()

    def run_statevector():
        qc = hhl_2x2_with_save(save_statevector=True)
        r = sim.run(transpile(qc, sim)).result()
        r.get_statevector(qc)

    bench_meas = benchmark(run_measure_4096, warmup=1, repeat=5)
    bench_sv = benchmark(run_statevector, warmup=1, repeat=5)
    print(f"  build+transpile+run(shots=4096): {bench_meas}")
    print(f"  build+transpile+run(statevector): {bench_sv}")
    results["benchmark_measure_4096"] = bench_meas
    results["benchmark_statevector"] = bench_sv

    # 保存
    out = Path(__file__).resolve().parent.parent / "results" / "python_hhl.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
