# 4He at L = 3: symmetry-adapted register of total momentum zero (9 system qubits)

"The build" means the scripts that constructed and verified these circuits for the paper; they are not part of this release.

Gate-level circuits of the 9-qubit register of the 157 states of total momentum zero in the two-excitation space
(the "symmetric 4He" column of Table VII). Every file is the transpiled circuit the build simulated
({u, cx}, Qiskit 2.3.1, `transpile(..., basis_gates=["u", "cx"], optimization_level=1, qubits_initially_zero=False)`),
written with `qiskit.qasm2.dumps`. Files above 1 MB are gzipped (`*.qasm.gz`). `manifest.json` lists every file with
its qubit count, CNOT and `u` counts, depth, the build's and the paper's count, dt / T / n, the global phase of each
piece, the machine-readable `readout` (type, basis indices, phase, pointer to the reference value) that `simulate/` runs, and the result of the export check (`verification`). Reference values: `data/he4_L3_sym/`.

## Loading

`u(theta, phi, lambda)` is Qiskit's `u` (= `u3`). Every file is self-contained. Load with

```python
from qiskit import qasm2
qc = qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)   # text: the file, gunzipped if .gz
```

or with `simulate/circuits.load(path)`, which also reads `.gz` files.

Two layouts:
- flat: one transpiled circuit verbatim (steps, S2 Hadamard tests, walk, U_H, controlled walk, reflection);
- sequence (`cheb_walk_chain_*`, `gqsp_hadamard_*`): each transpiled piece is defined once as an OpenQASM 2 `gate`
  acting on all qubits of the file, and the circuit is the list of calls. `qc.decompose()` unrolls it.

OpenQASM 2 drops global phases. They are listed per piece in `manifest.json` (`pieces[].global_phase`); a piece of
the build equals `exp(i global_phase)` times the gates in the file.

## Qubits and bit order

Basis index = sum_q bit_q 2^q (Qiskit little-endian; qubit 0 is the least significant bit).

| qubits | role |
|---|---|
| 0-2 | pair label: 0..5 = species pairs (01), (02), (03), (12), (13), (23); 7 = HF |
| 3-8 | momentum p of the pair, three digits (k_c mod 3) on (3+2c, 4+2c), value a + 2b (code 11 unused) |
| system code | label + 8 x (6-bit momentum code); the 157 valid codes are in `data/he4_L3_sym/register.json` |
| steps: 9 / 10, 11 | control (Hadamard tests: test qubit) / work qubits (enter and leave at 0) |
| walk: 9, 10 / 11, 12 | block-encoding qubit a, branch qubit b / work qubits |
| controlled walk: 13 | control |
| GQSP: 13 / 14 | processing qubit c / Hadamard-test ancilla h |

Reference state HF = system code 7 (label 7, p = 0), a basis state: `x` on qubits 0, 1, 2. The S2 Hadamard tests contain
these `x` gates inside the transpiled circuit; the chain and GQSP files apply them as plain `x` statements first
(the build set this basis state directly).

## Files and readout

| files | circuit | readout |
|---|---|---|
| `controlled_S2_dt*_n{1,2,3}`, `controlled_S4_dt0.0075_n{1,2}` | n merged controlled S2 / S4 steps (pieces A, B), 12 qubits | control 1: S^n on the valid codes; control 0: identity (times exp(-i phase) unless the phase is restored) |
| `hadamard_S2_dt{0.0075,0.00375}_n{1,2,3}` | Hadamard test of n S2 steps, 12 qubits | A0 = amp[7], A1 = amp[7 + 2^9] before any final gate; C(n) = 2 conj(A0) A1 = <HF\|S2^n\|HF> |
| `walk`, `U_H`, `controlled_walk`, `reflection` | walk W = R U_H, block encoding U_H (block (H - c)/alpha), W controlled on qubit 13, R = 2\|00><00\| - 1 on (9, 10) | see below |
| `cheb_walk_chain_k20` | 20 calls of the walk from \|0, HF> | psi_k = state after the k-th call; mu_2k = 2 P(a = b = 0) - 1; mu_2k+1 = Re[exp(i phi_UH) <psi_k\|U_H'\|psi_k>] |
| `gqsp_hadamard_n01..n14`, `_T0.05`, `_T0.3` | GQSP evolution e^{-i(H - c)T}, Hadamard-test version, 15 qubits | A0 = amp[7], A1 = amp[7 + 2^14]; C = 2 conj(A0) A1 = <HF\|F(x)/s\|HF> ~ <HF\|e^{-i(H - c)T}\|HF> |

- Time steps. Steps and S2 Hadamard tests: dt = 0.0075 MeV^-1 (a1/a2 of the build: the step bound
  0.99 pi/(E_max - E0) = 0.0074984 rounded) and dt/2 = 0.00375. GQSP: dt = 0.99 pi/(E_max - E0) = 0.00749837 MeV^-1,
  T = n dt, n = 1..14, plus T = 0.05 and 0.3 (Table VII). E0 = -145.06733, E_max = 269.71307,
  alpha = 229.77885, c = 55.46992 MeV.
- C(n) = 2 conj(A0) A1 from the amplitudes of |0>|ref> and |1>|ref> before any final gate (on hardware: final H on
  the test qubit gives Re C = P(0, HF) - P(1, HF); S^dag then H gives Im C). The global phases cancel in C. The GQSP
  value approximates <HF|e^{-i(H - c)T}|HF>, so <HF|e^{-iHT}|HF> ~ exp(-i c T) C.
- Chebyshev moments: mu_n = <HF|T_n((H - c)/alpha)|HF>, n = 0..41. Even orders are probabilities (no phase).
  Odd orders need U_H with its phase: U_H = exp(i phi_UH) U_H' with U_H' the gates of `U_H.qasm` and
  phi_UH = 2.1714688825171677 (`manifest.json`); the walk's phase does not enter. A moment read as an amplitude,
  <0, HF|W^k|0, HF> = T_k moment, needs exp(i k phi_W), phi_W = 5.313061536106961.
- Controlled steps and the controlled walk: the dropped phase multiplies both control branches; it matters only
  when amplitudes are compared with the operator, not inside a Hadamard test.
- The reflection file is 1 - 2|00><00|; R = exp(i pi) times it. In Sinv = R cW R the two phases cancel.

## Provenance (build scripts, not part of this release)

`a1_p0_steps.py` (steps), `a2_p0_krylov.py` (S2 Hadamard tests), `a3_p0_walk.py` (walk, U_H, controlled
walk), `a4_p0_cheb.py` (walk chain), `a3b_gqsp.py` (GQSP at T = 0.05, 0.3; T = 0.3 is written here piecewise, as
`run_qkrylov.py` builds it, with the same CNOT count 1 066 902 as the build's whole-circuit transpile),
`run_qkrylov.py` (GQSP n = 1..14 and T = 0.05). Every file was reloaded, unrolled and simulated (Aer
statevector); the CNOT counts equal the build's and the values agree with the build's (`verification` in
`manifest.json`): GQSP overlaps to 6.7e-12, S2 Hadamard-test overlaps to 1.9e-13, Chebyshev moments (orders 0..41) to
2.0e-12; the 157 x 157 blocks of the controlled steps equal the exact step operators to 4.3e-11 (build: 4.1e-11).
The small excess over the build's own numbers comes from `qasm2.dumps` writing angles within 1e-12 of simple
multiples of pi symbolically (largest parameter change 2e-13).
