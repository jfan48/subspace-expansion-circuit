# Site encoding (Jordan-Wigner): deuteron and triton on the L = 2 lattice

Hand-built circuits of the site encoding: the gate units, the controlled product-formula steps and the deuteron's
controlled walk of Table XI, and the Hadamard tests of the fourth-order step behind Figs. 13 and 14. No transpiler
was used: every gate was written by `trotter_circuits.py` (steps, units) and
`qubitized_walk.py` (walk), scripts not part of this release. `manifest.json` lists every file with its qubit and CNOT counts, the
paper's count where Table XI prints one, and a machine-readable `readout`; reference values are in
`data/site_encoding/`.

## Gate set and loading

The files use `h`, `s`, `sdg`, `rz`, `ry`, `u1`, `x` and `cx` from `qelib1.inc` (no `u`), with
rz(t) = exp(-i t Z/2) and u1(f) = diag(1, e^{if}). They load with `qiskit.qasm2.loads(text)`. The triton
Hadamard tests define the controlled step once as the gate `cs4_step` and call it n times. No global phase is lost: the identity terms of the Hamiltonian are applied as one
`u1` gate on the ancilla at the end of each controlled step.

## Qubits, bit order and reference state

- Basis index = sum_q bit_q 2^q (Qiskit little-endian).
- n_orb spin-isospin orbitals (deuteron 2, triton 3), 8 sites each. Orbital o occupies qubits 8o .. 8o+7, and site
  (i, j, k), i, j, k in {0, 1}, is qubit 8o + i + 2j + 4k (occupation 1 = |1>).
- Jordan-Wigner ordering follows the qubit index (orbital-major). A hop between sites p < q of one orbital carries
  the string Z on the qubits strictly between them.
- The controlled steps and the Hadamard tests use n_qubits = 8 n_orb + 1. The ancilla (control, Hadamard-test
  ancilla) is qubit 8 n_orb: 16 for the deuteron, 24 for the triton.
- Reference state: every nucleon on site (0,0,0), that is qubits 8o occupied, prepared by `x` on qubits 0, 8 (and 16
  for the triton). Its basis index is 257 (deuteron) or 65 793 (triton). This is a localized basis state, not a
  Hartree-Fock state.
- Hamiltonian (convention A of the paper, doubled L = 2 hops), per orbital:
  H_kin = sum_q 6 kappa n_q - kappa sum_<ij> (X_i Zbar X_j + Y_i Zbar Y_j). Interaction:
  H_I = V2 sum_sites sum_{o<o'} n_o n_o', plus V3 n_1 n_2 n_3 for the triton. Constants: kappa = 5.18, V2 = -14.25,
  V3 = -17.04 MeV.

## Files

| files | content | qubits | CNOT (Table XI) |
|---|---|---|---|
| `unit_hop_s{0,1,3}` | controlled hop C_a-exp(+i theta (X_i Zbar X_j + Y_i Zbar Y_j)), i = 0, j = s + 1, string qubits 1..s, control a = s + 2, theta = kappa x 0.01 | s + 3 | 5 / 7 / 11 |
| `unit_site_diagonal_{2,3}orb` | controlled phase polynomial of one site, exp(-i dt sum_S c_S Z_S) (identity term excluded), site qubits 0..n-1, control n, dt = 0.01 | 3 / 4 | 6 / 14 |
| `deuteron_C{S1,S2,S4}_dt0.01_n1` | controlled first-, second-, fourth-order step, dt = 0.01 MeV^-1 | 17 | 220 / 392 / 1 272 |
| `triton_C{S1,S2,S4}_dt0.01_n1` | the same for the triton | 25 | 370 / 628 / 2 108 |
| `{deuteron,triton}_CS2_dt0.01_n3_merged` | three consecutive S2 steps, adjacent hopping half-layers merged | 17 / 25 | 832 / 1 368 |
| `deuteron_controlled_walk` | controlled qubitized walk of the 72-Pauli-string LCU: control 0, selection register 1-7, scratch 8-14, system 15-30 (system qubit q on 15 + q) | 31 | 549 |
| `hadamard_deuteron_S4/deuteron_S4_dt0.01_n01..n50` | Hadamard test of n controlled S4 steps (one merged n-step circuit per file), Fig. 13 | 17 | 1 272 + 1 100 (n - 1) |
| `hadamard_triton_S4/triton_S4_dt0.01_n01..n70` | Hadamard test of n controlled S4 steps: `cs4_step` called n times, as simulated; Fig. 14 | 25 | 2 108 n |

The step files hold the controlled step only. The Hadamard-test files also contain the preparation: `x` on the
reference qubits, `h` on the ancilla.

## Readout

C(n) = 2 conj(A0) A1, with A0 = <0|<ref|psi> and A1 = <1|<ref|psi> the amplitudes of the ancilla value 0 and 1 with
the system in the reference state, read before any final gate. C(n) is <ref|U^n|ref>. On hardware, a final `h` on
the ancilla gives Re C = P(0, ref) - P(1, ref), and `sdg` then `h` gives Im C.

Reference values:
- `data/site_encoding/deuteron_S4_dt0.01_overlaps.json` and `triton_S4_dt0.01_overlaps.json`: the overlaps of the
  paper's runs (results_simulation.json).
- `data/site_encoding/table_xi_steps.json`: the ideal overlaps <ref|U^n|ref> of the S1, S2 and S4 steps, computed in
  the fixed-particle-number sector independently of the circuits.

For the step files the runner prepends the preparation (`readout.prep`). The gate units are compared as operators
with the exact controlled exponentials built from the parameters in their `readout`. The 31-qubit walk is counted
only.
