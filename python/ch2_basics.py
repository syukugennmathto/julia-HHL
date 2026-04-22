"""第2章: 量子コンピュータの基礎 - 全スニペットを実行して結果を出力"""
import time
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

SHOTS = 1024
sim = AerSimulator()


def run_and_report(title, qc):
    t0 = time.perf_counter()
    job = sim.run(qc, shots=SHOTS)
    result = job.result()
    counts = result.get_counts()
    dt = (time.perf_counter() - t0) * 1000
    print(f"[{title}] ({dt:.2f} ms) counts={counts}")
    return counts, dt


# 1量子ビットの測定回路
qc = QuantumCircuit(1, 1)
qc.measure(0, 0)
run_and_report("1Q measure |0>", qc)

# 4量子ビットの測定回路
qc = QuantumCircuit(4, 4)
qc.measure([0, 1, 2, 3], [0, 1, 2, 3])
run_and_report("4Q measure |0000>", qc)

# Xゲート
qc = QuantumCircuit(1, 1)
qc.x(0)
qc.measure(0, 0)
run_and_report("X gate", qc)

# H then Z then measure
qc = QuantumCircuit(1, 1)
qc.h(0)
qc.z(0)
qc.measure(0, 0)
run_and_report("H-Z", qc)

# 1Q Hadamard
qc = QuantumCircuit(1, 1)
qc.h(0)
qc.measure(0, 0)
run_and_report("H (1Q)", qc)

# 2Q Hadamard
qc = QuantumCircuit(2, 2)
qc.h(0)
qc.h(1)
qc.measure([0, 1], [0, 1])
run_and_report("H (2Q)", qc)

print("\nCh.2 done.")
