"""第3章: HHLに必要な量子回路 - エンコーディングとQPE骨格"""
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

# --- 基底エンコーディング: 5 = 0b0101 -> |0101>
qc = QuantumCircuit(4)
qc.x(0)
qc.x(2)
print("=== 基底エンコーディング (5 = |0101>) ===")
print(qc.draw())

# --- テンソル積エンコーディング: |+> x |+>
qc = QuantumCircuit(2)
qc.h(0)
qc.h(1)
print("\n=== テンソル積エンコーディング (|+>|+>) ===")
print(qc.draw())

# --- 量子カーネル回路
n_features = 2
x = ParameterVector('x', n_features)
qc = QuantumCircuit(n_features)
for i in range(n_features):
    qc.h(i)
    qc.rz(2.0 * x[i], i)
print("\n=== 量子カーネル回路 ===")
print(qc.draw())


# --- 位相推定アルゴリズム（骨格）
def quantum_phase_estimation(num_ancilla):
    n_qubits = 1
    qc = QuantumCircuit(num_ancilla + n_qubits, num_ancilla)
    for i in range(num_ancilla):
        qc.h(i)
    qc.measure(range(num_ancilla), range(num_ancilla))
    return qc


print("\n=== QPE骨格 (num_ancilla=3) ===")
print(quantum_phase_estimation(3).draw())

print("\nCh.3 done.")
