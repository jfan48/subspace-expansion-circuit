# subspace-expansion-circuit

Gate-level circuits of the paper

> J. Fan and T. Papenbrock, *Ground-state energies from real-time and Chebyshev subspaces for nuclear lattice
> Hamiltonians* (2026), arXiv link to be added.

The repository holds every circuit that the paper simulated gate by gate, as OpenQASM 2 files, together with the
values each circuit reproduces and a script that reruns them. Table and figure numbers below refer to the paper.

## Contents

| folder | system and encoding | qubits | circuits | paper |
|---|---|---|---|---|
| `circuits/site_encoding` | deuteron and triton on the L = 2 lattice, second quantization, one qubit per site and spin-isospin state (Jordan-Wigner) | 17, 25 (walk 31) | gate units; controlled first-, second- and fourth-order product-formula steps; controlled walk of the deuteron; Hadamard tests of the fourth-order step at dt = 0.01 MeV^-1 for n <= 50 (deuteron) and n <= 70 (triton) steps | Table XI; Figs. 13, 14, Table VI |
| `circuits/he4_L3_sym` | 4He at L = 3, first-quantized register of the 157 states of total momentum zero (9 system qubits) | 12-15 | controlled second- and fourth-order steps and their Hadamard tests; walk, block encoding, controlled walk; chain of walks for the Chebyshev moments of orders 0-41; Hadamard tests of qubitized evolution (GQSP) for n = 1-14 time steps and for T = 0.05, 0.3 MeV^-1 | Table VII; Fig. 15 (a), Table VIII |
| `circuits/o16_L7_sym` | 16O at L = 7, first-quantized register of the 109 states symmetric under permutations of the spin-isospin states (7 system qubits) | 8-10 | controlled evolution in the exact eigenbasis and its Hadamard tests; walk, block encoding, controlled walk; chain of walks for the Chebyshev moments of orders 0-81; GQSP Hadamard tests for n = 1-37 and T = 0.05, 0.3 MeV^-1 | Table VII; Fig. 15 (b), Table VIII |
| `circuits/he4_L3_species`, `circuits/o16_L7_species` | species registers of the same nuclei (24 system qubits) | 35, 30 | the controlled projected pair piece, in three sections and as one transpiled circuit | Table VII |

Each folder has a `README.md` (qubit roles, reference state, time steps, readout) and a `manifest.json` with one
entry per file: description, qubit and CNOT counts, the count the paper prints, the global phase of each piece, and
a machine-readable `readout` that points to the reference value in `data/`.

Not included: the circuits of the 8Be register, which the paper counts in parts but does not simulate, and the
scripts that build the Hamiltonians and the circuits and analyze the estimates.

## Format

- OpenQASM 2.0, one self-contained file per circuit; files above 1 MB are gzipped (`*.qasm.gz`).
- Register circuits: gate set {u, cx}, transpiled with Qiskit 2.3.1 (`optimization_level=1`). `u` is Qiskit's `u`
  (= `u3`), so load them with
  ```python
  from qiskit import qasm2
  qc = qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
  ```
- Site-encoding circuits: written gate by gate without a transpiler, gate set h, s, sdg, rz, ry, u1, x, cx of
  `qelib1.inc`; `qasm2.loads(text)` loads them.
- A circuit that repeats a piece (GQSP sequences, chains of walks, the triton step) defines the piece once as an
  OpenQASM 2 `gate` and calls it; `qc.decompose()` unrolls it. The CNOT counts are those of the unrolled circuits.
- Bit order: Qiskit little-endian, basis index = sum over qubits q of bit_q 2^q.
- OpenQASM 2 has no global phase. The manifests list the global phase of every transpiled piece. It cancels in every
  Hadamard test and enters only the odd Chebyshev moments and direct comparisons with an operator, as the folder
  READMEs explain.

## Readout

A Hadamard-test file contains the preparation of the reference state, the Hadamard gate on the test qubit, and the
controlled evolution, without the final basis change. The overlap is
C = 2 conj(A0) A1, where A0 and A1 are the amplitudes of |0>|ref> and |1>|ref> (test qubit first). On a device, a
final H on the test qubit gives Re C = P(0, ref) - P(1, ref), and S^dagger followed by H gives Im C, with P the
joint probabilities of the test-qubit outcome and the reference state. The Chebyshev moments are read from the
states of a chain of walks as described in the folder READMEs.

## Rerunning

```
pip install -r requirements.txt
python3 simulate/check_all.py --full --cx-only   # CNOT counts of all 221 files against the manifests and the paper (~15 s)
python3 simulate/check_all.py --quick            # counts plus statevector simulation of a subset (~4 min on a laptop)
python3 simulate/check_all.py --full             # simulation of every file that fits in a statevector
python3 simulate/reproduce.py circuits/o16_L7_sym/gqsp_hadamard_n05.qasm.gz
```

`simulate/circuits.py` loads, counts and simulates the files (statevector, Qiskit Aer); long circuits are simulated
one piece at a time, so the GQSP circuits with up to two million CNOTs run without being unrolled. Three files are
too large for a statevector and are only counted: the 31-qubit controlled walk of the deuteron, and the 35- and
30-qubit pair pieces as single circuits (their sections are simulated).

## Reference values

`data/<folder>/` holds, as JSON, the values from which the paper's tables and figures were made: overlaps, Chebyshev
moments, operator blocks, and the energies and normalizations of the registers. Each manifest entry names the file
and key of its reference value.

## License

The code in `simulate/` is released under the MIT License (`LICENSE`). The circuit files in `circuits/` and the
reference values in `data/` are released under the Creative Commons Attribution 4.0 International License
(`LICENSE-CC-BY-4.0.txt`); please cite the paper above when you use them.
