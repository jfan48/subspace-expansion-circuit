# 4He at L = 3: species register, controlled projected pair piece (35 qubits)

"The build" means the scripts that constructed and verified these circuits for the paper; they are not part of this release.

The full-size controlled projected pair piece of the species register (Table VII, "pair piece" row, 4He species
column: 5 436 CNOT): controlled exp(-i tau P A^(01) P) on registers 0, 1, with registers 2, 3 as spectators,
tau = 0.0075 MeV^-1, A = c2_eff V2 + (u x 1 + 1 x u)/3; P projects on the two-excitation space, and the parts
n = 0, 1, 2 act when 0, 1, 2 of registers 2, 3 are excited. Transpiled with Qiskit 2.3.1 ({u, cx}, optimization
level 1, `qubits_initially_zero=False`), written with `qiskit.qasm2.dumps`. Load with
`qiskit.qasm2.loads(text, custom_instructions=qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS)` (`u` = Qiskit's `u` = `u3`).

| file | content | CNOT |
|---|---|---|
| `pair_piece_01_sections.qasm` | gates `pre` (compute z_u, z_v), `mid` (parts n = 0, 2, 1, own flags z_s, z_t, condition flags), `post` (uncompute z_u, z_v), each transpiled separately in the 35-qubit register; the circuit is pre, mid, post | 84 + 5 268 + 84 = 5 436 |
| `pair_piece_01_whole_transpile.qasm` | the same circuit pre + mid + post transpiled in one call; the count printed in Table VII | 5 436 |

## Verification: in sections, not as a whole statevector

A 35-qubit statevector is beyond reach, so the build (`b6_pair_piece.py`) verified the full-size circuit
by simulating its sections with registers 2, 3 in basis states, for which the 35-qubit state factorizes between the
sections: (1) `pre` and `post`, restricted to the qubits they touch (registers 2, 3, z_u, z_v and the idle qubits
0-3 that the multi-controlled X synthesis borrows as dirty ancillas), on all 64 x 64 codes of registers 2, 3 with a
random state on the borrowed qubits; (2) `mid`, which touches no qubit of registers 2, 3, remapped to its 23 qubits
and simulated per class (none / one / both of registers 2, 3 excited) against exp(-i tau Pi_n A Pi_n). These checks
were rerun on the reloaded gate definitions of `pair_piece_01_sections.qasm` (`verification` in `manifest.json`).
The whole-transpile file was counted only; neither the build nor this export simulated it.

## Qubits and bit order

Basis index = sum_q bit_q 2^q (little-endian). Register r (species r) on qubits 6r..6r+5: digit c of the momentum k_r
(k_c mod 3) on (6r+2c, 6r+2c+1) with value a + 2b (code 11 unused); k = 0 is 000000, so HF is the all-zero state.
Control 24; canonicalization work qubits 25-29; flags z_s 30, z_t 31, z_u 32, z_v 33 (z_r = [register r not at
k = 0]); condition flag 34. Work and flag qubits enter and leave at 0. Global phases (dropped by OpenQASM 2) are in
`manifest.json`: `mid` and the whole transpile carry 2.419057 rad; they matter only for a direct comparison of
amplitudes with the operator. Reference data: `data/he4_L3_species/pair_piece.json` (layout, counts, the build's
section-check deviations).

`simulate/check_all.py` reruns these section checks standalone (`readout` of the sections file and `data/he4_L3_species/pair_piece_refs.json`: random inputs of each class and their exact outputs); `--quick` uses 64 random (u, v) codes and one input per class, `--full` all 4 096 codes and every input.
