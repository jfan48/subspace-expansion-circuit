OPENQASM 2.0;
include "qelib1.inc";
// 4He at L = 3, symmetry-adapted P = 0 register (9 system qubits, 157 valid codes), 14-qubit frame of the GQSP signals; acts on qubits 9, 10 only
// reflection R = 2|00><00| - 1 on the walk ancillas (a 9, b 10), up to the dropped global phase (= pi)
// inverse signal: Sinv = R cW R = |0><0|_c x 1 + |1><1|_c x W^dag (cW = controlled walk, c = qubit 13)
// gate u(theta,phi,lambda) = Qiskit u (= u3); load with qiskit.qasm2.loads(text, custom_instructions=qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
qreg q[14];
u(pi,0,pi) q[9];
u(pi/2,-pi,-pi) q[10];
cx q[9],q[10];
u(pi,0,pi) q[9];
u(pi/2,0,0) q[10];
