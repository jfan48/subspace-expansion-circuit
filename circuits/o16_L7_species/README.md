# 16O at L = 7: species register, controlled projected pair piece (30 qubits)

"The build" means the scripts that constructed and verified these circuits for the paper; they are not part of this release.

The full-size controlled projected pair piece of the species register (Table VII, "pair piece" row, 16O species
column: 61 604 CNOT): controlled exp(-i tau P A^(01) P) on registers 0, 1, with registers 2, 3 as spectators,
tau = 0.99 pi/(21.17817 + 88.93615) = 0.0282449 MeV^-1 (step bound of the two-excitation space); A = A2' +
(O' x 1 + 1 x O')/3, the folded pair operator; P projects on the two-excitation space, and the parts n = 0 (205
pair states, through the pair compression and the 8-qubit eigenbasis), n = 1 (25 states) and n = 2 (phase) act when
the other two registers have total rank 0, 1, 2. Transpiled with Qiskit 2.3.1 ({u, cx}, optimization level 1,
`qubits_initially_zero=False`), written with `qiskit.qasm2.dumps`, gzipped. Load with
`qiskit.qasm2.loads(text, custom_instructions=qiskit.qasm2.LEGACY_CUSTOM_INSTRUCTIONS)` (`u` = Qiskit's `u` = `u3`).

| file | content | CNOT |
|---|---|---|
| `pair_piece_01_sections.qasm.gz` | gates `pre` (compute z_u, z_v), `mid` (parts n = 0, 2, 1, own flags, condition flags), `post` (uncompute z_u, z_v), each transpiled separately in the 30-qubit register; the circuit is pre, mid, post | 84 + 61 436 + 84 = 61 604 |
| `pair_piece_01_whole_transpile.qasm.gz` | the same circuit pre + mid + post transpiled in one call; the count printed in Table VII | 61 604 |

## Verification: in sections, not as a whole statevector

A 30-qubit statevector (16 GB) was not simulated. The build (`b2_pair_piece.py`, `results_b2_01.json`)
verified the full-size circuit by simulating its sections with registers 2, 3 in basis states, for which the 30-qubit
state factorizes between the sections: (1) `pre` and `post`, restricted to the qubits they touch (registers 2, 3,
z_u, z_v and the idle qubits 0-3 borrowed as dirty ancillas by the multi-controlled X synthesis), on all 64 x 64 codes
of registers 2, 3 with a random state on the borrowed qubits; (2) `mid`, which reads registers 2, 3 only through bit 5
(qubits 17, 23), remapped to its 20 qubits and simulated per class of (register 2, register 3) basis states, two
representatives per class, against exp(-i tau Pi_n A Pi_n). These checks were rerun on the reloaded gate definitions
of `pair_piece_01_sections.qasm.gz` (`verification` in `manifest.json`). The whole-transpile file was counted only;
neither the build nor this export simulated it.

## Qubits and bit order

Basis index = sum_q bit_q 2^q (little-endian). Register r (species r) on qubits 6r..6r+5, bit j of the 6-bit string
code on qubit 6r+j: HF = 000000; rank 1 (hole h in 0..3, particle p in 0..2): bits 1..0 = p + 1, bits 3..2 = h,
bits 5..4 = 0; rank 2 (hole pair hh in 0..5, particle pair pp in 0..2): bit 5 = 1, bits 4..2 = hh, bits 1..0 = pp.
HF of the nucleus is the all-zero state. Control 24; flags z_s 25, z_t 26, z_u 27, z_v 28 (z_r = [register r not
HF]); condition flag 29; flags enter and leave at 0. Global phases (dropped by OpenQASM 2) are in `manifest.json`:
`mid` and the whole transpile carry 0.205942 rad; they matter only for a direct comparison of amplitudes with the
operator. Reference data: `data/o16_L7_species/pair_piece.json` (layout, counts, the build's section-check
deviations).

`simulate/check_all.py` reruns these section checks standalone (`readout` of the sections file and `data/o16_L7_species/pair_piece_refs.json`: random inputs of each class and their exact outputs); `--quick` uses 64 random (u, v) codes and one input per class, `--full` all 4 096 codes and every input.
