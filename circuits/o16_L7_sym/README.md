# 16O at L = 7: species-symmetric register (7 system qubits)

"The build" means the scripts that constructed and verified these circuits for the paper; they are not part of this release.

Gate-level circuits of the 7-qubit register of the 109 states of the two-excitation space that are symmetric under
permutations of the four spin-isospin species (the "symmetric 16O" column of Table VII). Every file is the
transpiled circuit the build simulated ({u, cx}, Qiskit 2.3.1, `transpile(..., basis_gates=["u", "cx"],
optimization_level=1, qubits_initially_zero=False)`), written with `qiskit.qasm2.dumps`. Files above 1 MB are
gzipped (`*.qasm.gz`). `manifest.json` lists every file with its qubit count, CNOT and `u` counts, depth, the build's
and the paper's count, dt / T / n, the global phase of each piece, the machine-readable `readout` (type, basis indices, phase, pointer to the reference value) that `simulate/` runs, and the result of the export check (`verification`). Reference values: `data/o16_L7_sym/`.

## Loading

`u(theta, phi, lambda)` is Qiskit's `u` (= `u3`). Every file is self-contained. Load with

```python
from qiskit import qasm2
qc = qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)   # text: the file, gunzipped if .gz
```

or with `simulate/circuits.load(path)`, which also reads `.gz` files.

Flat files hold one transpiled circuit verbatim; in `cheb_walk_chain_k40` and `gqsp_hadamard_*` each transpiled piece
is defined once as an OpenQASM 2 `gate` on all qubits of the file and the circuit is the list of calls
(`qc.decompose()` unrolls it). OpenQASM 2 drops global phases; they are listed per piece in `manifest.json`
(`pieces[].global_phase`): a piece of the build equals `exp(i global_phase)` times the gates in the file.

## Qubits and bit order

Basis index = sum_q bit_q 2^q (Qiskit little-endian; qubit 0 is the least significant bit).

| qubits | role |
|---|---|
| 0-6 | system: code = index of the species-permutation orbit (orbit basis: normalized sum of the V_2 states related by species permutations); codes 0..108 valid, 109..127 unused |
| evolution: 7 | control (Hadamard tests: test qubit) |
| walk: 7 | walk ancilla a (block-encoding qubit) |
| controlled walk: 8 | control |
| GQSP: 7 / 8 / 9 | walk ancilla a / processing qubit c / Hadamard-test ancilla h |

Reference state HF = code 0: the all-zero basis state, no preparation.

## Files and readout

| files | circuit | readout |
|---|---|---|
| `controlled_evolution_tau0.030253` | controlled exact-eigenbasis evolution exp(-i dt H_s): W^T, controlled diagonal, W (W = eigenvectors on codes 0..108); 14 892 CNOT for any tau | control 1: exp(-i dt H_s); control 0: identity (up to the dropped phase) |
| `hadamard_evolution_chain_n{1,2,3}`, `_single_n{2,3}` | Hadamard test of n consecutive evolutions of dt (chain) or one evolution of n dt (single; single n = 1 is chain n = 1) | A0 = amp[0], A1 = amp[2^7] before any final gate; C(n) = 2 conj(A0) A1 = <HF\|e^{-i n dt H_s}\|HF> |
| `walk`, `U_H`, `controlled_walk`, `reflection` | W = Z_a U_H; U_H = W (Ry(theta) Z_a) W^T with block (H_s - c)/alpha; W controlled on qubit 8; R = Z_a | see below |
| `cheb_walk_chain_k40` | 40 calls of the walk from \|0_a, HF> | psi_k = state after the k-th call; mu_2k = 2 P(a = 0) - 1; mu_2k+1 = Re[exp(i phi_UH) <psi_k\|U_H'\|psi_k>] |
| `gqsp_hadamard_n01..n37`, `_T0.05`, `_T0.3` | GQSP evolution e^{-i(H_s - c)T}, Hadamard-test version, 10 qubits | A0 = amp[0], A1 = amp[2^9]; C = 2 conj(A0) A1 = <HF\|F(x)/s\|HF> ~ <HF\|e^{-i(H_s - c)T}\|HF> |

- Time step: dt = 0.99 pi/(E_max - E0) = 0.0302532944 MeV^-1 with the spectral range of this register,
  E0 = -88.93615, E_max = 13.86841 MeV (the build's a1, a2 and `run_qkrylov.py`); GQSP at T = n dt, n = 1..37, plus
  T = 0.05 and 0.3 (Table VII). alpha = 51.40228, c = -37.53387 MeV.
- C(n) = 2 conj(A0) A1 from the amplitudes of |0>|ref> and |1>|ref> before any final gate (on hardware: final H on
  the test qubit gives Re C = P(0, HF) - P(1, HF); S^dag then H gives Im C). Global phases cancel in C. The GQSP value
  approximates <HF|e^{-i(H_s - c)T}|HF>, so <HF|e^{-iH_s T}|HF> ~ exp(-i c T) C.
- Chebyshev moments mu_n = <HF|T_n((H_s - c)/alpha)|HF>, n = 0..81. Even orders are probabilities (no phase). Odd
  orders need U_H with its phase: U_H = exp(i phi_UH) U_H', phi_UH = 5.066232 (exact value in `manifest.json`); the
  walk's phase does not enter. A moment read as an amplitude <0, HF|W^k|0, HF> needs exp(i k phi_W), phi_W = phi_UH.
- Controlled evolution and controlled walk: the dropped phase multiplies both control branches; it matters only
  when amplitudes are compared with the operator, not inside a Hadamard test.
- The reflection R = Z_a was not transpiled separately by the build: it enters as plain `z` gates in
  Sinv = Z_a cW Z_a, which is transpiled as part of the signal `sig_i`; `reflection.qasm` is that `z` gate.

## Provenance (build scripts, not part of this release)

`a1_sym_circuits.py` (controlled evolution, walk, U_H, controlled walk), `a2_sym_tests.py` (Hadamard tests,
walk chain; the 82 moments of `results_a2.json` feed Fig. 15 (b)), `a3_sym_gqsp.py` (GQSP at T = 0.05, 0.3; T = 0.3
is written here piecewise, as `run_qkrylov.py` builds it, with the same CNOT count 715 080 as the build's
whole-circuit transpile), `run_qkrylov.py` (GQSP n = 1..37 and T = 0.05). Every file was reloaded, unrolled and
simulated (Aer statevector); the CNOT counts equal the build's and the values agree with the build's
(`verification` in `manifest.json`): GQSP overlaps to 2.3e-11, Hadamard-test overlaps of the evolution to 4.8e-13,
Chebyshev moments (orders 0..81) to 1.4e-11 (to 2.2e-12 from the exact moments); the evolution, walk, U_H and
controlled walk equal their exact operators on all basis inputs to 1.2e-11 (build: 1.1e-11). The excess over the
build's own numbers comes from `qasm2.dumps` writing angles within 1e-12 of simple multiples of pi symbolically
(largest parameter change 2.5e-13; the same pieces simulated before writing reproduce the build to 8e-14).
