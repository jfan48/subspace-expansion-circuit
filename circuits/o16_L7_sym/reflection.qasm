OPENQASM 2.0;
include "qelib1.inc";
// 16O at L = 7, species-symmetric register (7 system qubits, 109 valid codes, orbit basis), 9-qubit frame of the GQSP signals
// reflection R = 2|0><0|_a - 1 = Z on the walk ancilla (qubit 7); not transpiled separately by the build: it enters as untranspiled z gates in Sinv = Z_a cW Z_a, which is transpiled as part of the signal sig_i
// gate u(theta,phi,lambda) = Qiskit u (= u3); load with qiskit.qasm2.loads(text, custom_instructions=qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
qreg q[9];
z q[7];
