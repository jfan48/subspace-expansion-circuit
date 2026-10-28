"""Load, count and simulate the OpenQASM 2 circuits of this repository (needs only qiskit, qiskit-aer, numpy).

    load(path)                 QuantumCircuit; reads .qasm or .qasm.gz, accepts Qiskit's `u` gate
    count_cx(qc)               CNOT count, recursive over custom gates (memoized by gate name), no unrolling
    run(qc, psi0, ...)         statevector (Aer) of qc applied to psi0; piecewise: each custom gate call is one
                               piece (its definition), runs of standard gates are grouped into one piece, so the
                               two-million-CNOT GQSP circuits run without being unrolled
    hadamard_overlap(psi, ref, test)
                               C = 2 conj(A0) A1, A0 = <0|_test <ref| psi>, A1 = <1|_test <ref| psi>
    chebyshev_moments(chain, uh, uh_phase, ancillas, n_max)
                               mu_2k = 2 P(ancillas 0 | psi_k) - 1, mu_2k+1 = Re[exp(i uh_phase) <psi_k|U_H'|psi_k>],
                               psi_k = state after the k-th call of the walk gate in the chain file
    evaluate(entry, level)     runs the machine-readable "readout" of a manifest entry; returns a dict with the
                               computed values, the reference values from data/ and the largest deviation

Basis index = sum_q bit_q 2^q (Qiskit little-endian).  OpenQASM 2 drops global phases; where a value depends on
one, the readout carries it ("phase") and amplitudes are multiplied by exp(i phase).
"""
import gzip
import json
import os

import numpy as np
from qiskit import QuantumCircuit, qasm2
import qiskit_aer  # noqa: F401  (registers set_statevector / save_statevector on QuantumCircuit)
from qiskit_aer import AerSimulator

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASIC = {"u", "u1", "u2", "u3", "p", "rz", "ry", "rx", "h", "s", "sdg", "t", "tdg", "x", "y", "z", "id", "sx",
         "cx", "cz", "swap", "ccx", "barrier"}
_SIM = {}


def sim():
    if "s" not in _SIM:
        _SIM["s"] = AerSimulator(method="statevector", precision="double")
    return _SIM["s"]


def repo_path(path):
    return path if os.path.isabs(path) else os.path.join(REPO, path)


def read_text(path):
    path = repo_path(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    with open(path) as f:
        return f.read()


def load(path):
    path = repo_path(path)
    folder = os.path.dirname(path)
    return qasm2.loads(read_text(path), include_path=(folder,),
                       custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)


# ---------------------------------------------------------------------------------------------- counting
def count_cx(qc, _memo=None):
    memo = {} if _memo is None else _memo
    n = 0
    for inst in qc.data:
        op = inst.operation
        if op.name == "cx":
            n += 1
        elif op.num_qubits < 2 or op.name == "barrier":
            continue
        else:
            if op.name not in memo:
                memo[op.name] = count_cx(op.definition, memo)
            n += memo[op.name]
    return n


def count_ops(qc, _memo=None):
    """Number of gates of the unrolled circuit (no unrolling done)."""
    memo = {} if _memo is None else _memo
    n = 0
    for inst in qc.data:
        op = inst.operation
        if op.name in BASIC:
            n += 1
        else:
            if op.name not in memo:
                memo[op.name] = count_ops(op.definition, memo)
            n += memo[op.name]
    return n


# ---------------------------------------------------------------------------------------------- simulation
def basis_state(nq, idx):
    v = np.zeros(2 ** nq, dtype=complex)
    v[idx] = 1.0
    return v


def _run_flat(circ, psi):
    """One Aer run of a circuit of basic gates on psi (None: start from |0...0>)."""
    c = QuantumCircuit(circ.num_qubits)
    if psi is not None:
        c.set_statevector(psi)
    c.compose(circ, inplace=True)
    c.save_statevector()
    return np.asarray(sim().run(c, shots=1).result().get_statevector())


def pieces(qc):
    """Split qc into pieces: ('std', circuit of basic gates) or (gate name, definition circuit, qubits)."""
    out, block = [], None
    for inst in qc.data:
        op = inst.operation
        qs = [qc.find_bit(q).index for q in inst.qubits]
        if op.name in BASIC:
            if block is None:
                block = QuantumCircuit(qc.num_qubits)
            block.append(op, qs)
        else:
            if block is not None:
                out.append(("std", block, None))
                block = None
            out.append((op.name, op.definition, qs))
    if block is not None:
        out.append(("std", block, None))
    return out


def run(qc, psi0=None, callback=None, unroll_below=200_000):
    """Statevector of qc applied to psi0 (None = |0...0>).  callback(name, psi) after every custom-gate call.
    Circuits whose unrolled size is below `unroll_below` gates and that need no callback are unrolled and run in
    one Aer call (cheaper for many qubits); larger ones run piece by piece."""
    if callback is None and count_ops(qc) < unroll_below:
        names = [nm for nm in qc.count_ops() if nm not in BASIC]
        flat = qc
        while names:
            flat = flat.decompose(gates_to_decompose=names)
            names = [nm for nm in flat.count_ops() if nm not in BASIC]
        return _run_flat(flat, psi0)
    psi = psi0
    cache = {}
    for name, circ, qs in pieces(qc):
        if name == "std":
            psi = _run_flat(circ, psi)
        else:
            key = (name, tuple(qs))
            if key not in cache:
                sub = QuantumCircuit(qc.num_qubits)
                sub.compose(circ, qubits=qs, inplace=True)
                cache[key] = sub
            sub = cache[key]
            if set(sub.count_ops()) <= BASIC:
                psi = _run_flat(sub, psi)
            else:
                psi = run(sub, psi, unroll_below=0)
            if callback is not None:
                callback(name, psi)
    if psi is None:
        psi = basis_state(qc.num_qubits, 0)
    return psi


# ---------------------------------------------------------------------------------------------- readouts
def hadamard_overlap(psi, ref, test):
    A0, A1 = psi[ref], psi[ref + (1 << test)]
    return complex(2 * np.conj(A0) * A1)


def chebyshev_moments(chain, uh, uh_phase, ancillas, n_max, walk_gate="walk"):
    """Moments mu_0..mu_n_max from a chain file (prelude of basic gates, then calls of `walk_gate`) and the U_H
    circuit (gates of U_H.qasm, global phase uh_phase)."""
    nq = chain.num_qubits
    mask = sum(1 << q for q in ancillas)
    idx0 = np.nonzero((np.arange(2 ** nq) & mask) == 0)[0]
    mu = np.full(n_max + 1, np.nan)

    def moments(k, psi):
        if 2 * k <= n_max:
            mu[2 * k] = 2 * np.sum(np.abs(psi[idx0]) ** 2) - 1
        if 2 * k + 1 <= n_max:
            mu[2 * k + 1] = (np.exp(1j * uh_phase) * np.vdot(psi, run(uh, psi))).real
    pcs = pieces(chain)
    psi, k, cache = None, 0, {}
    for name, circ, qs in pcs:
        if name == "std":
            psi = _run_flat(circ, psi)
            continue
        if k == 0:
            psi = psi if psi is not None else basis_state(nq, 0)
            moments(0, psi)
        assert name == walk_gate, name
        if tuple(qs) not in cache:
            sub = QuantumCircuit(nq); sub.compose(circ, qubits=qs, inplace=True)
            cache[tuple(qs)] = sub
        psi = _run_flat(cache[tuple(qs)], psi)
        k += 1
        moments(k, psi)
    return mu, k


# ---------------------------------------------------------------------------------------------- references
def _get(obj, path):
    for p in path:
        obj = obj[p]
    return obj


def resolve(refspec):
    d = json.load(open(repo_path(refspec["file"])))
    if "path" in refspec:
        return _get(d, refspec["path"])
    if "re_path" in refspec:
        return [_get(d, refspec["re_path"]), _get(d, refspec["im_path"])]
    return d


def as_complex(v):
    if isinstance(v, (list, tuple)) and len(v) == 2 and not isinstance(v[0], (list, tuple)):
        return complex(v[0], v[1])
    return complex(v)


def as_complex_array(v):
    return np.array([as_complex(z) for z in v])


# ---------------------------------------------------------------------------------------------- exact units
_PAULI = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1])}


def pauli_dense(nq, ops):
    M = np.ones((1, 1), dtype=complex)
    for q in reversed(range(nq)):
        M = np.kron(M, _PAULI[ops.get(q, "I")])
    return M


def controlled_dense(U, a, nq):
    """|0><0|_a x 1 + |1><1|_a x U, a the top qubit (a = nq - 1), U on qubits 0..nq-2."""
    assert a == nq - 1
    d = 2 ** (nq - 1)
    C = np.eye(2 * d, dtype=complex)
    C[d:, d:] = U
    return C


def unit_exact(ro):
    nq, a = ro["nq"], ro["a"]
    ns = nq - 1
    if ro["kind"] == "controlled_hop":
        i, j, th = ro["i"], ro["j"], ro["theta"]
        zs = {r: "Z" for r in range(i + 1, j)}
        P = pauli_dense(ns, {i: "X", j: "X", **zs}) + pauli_dense(ns, {i: "Y", j: "Y", **zs})
        w, V = np.linalg.eigh(P)
        U = (V * np.exp(1j * th * w)) @ V.conj().T
    else:
        dt = ro["dt"]
        h = np.zeros(2 ** ns)
        x = np.arange(2 ** ns)
        for orbs, c in ro["coefficients"]:
            par = np.zeros(2 ** ns, dtype=int)
            for o in orbs:
                par ^= (x >> ro["site_qubits"][o]) & 1
            h += c * (1 - 2 * par)
        U = np.diag(np.exp(-1j * dt * h))
    return controlled_dense(U, a, nq)


# ---------------------------------------------------------------------------------------------- pair pieces
def _remap(circ, qubits):
    loc = {q: i for i, q in enumerate(qubits)}
    out = QuantumCircuit(len(qubits))
    for inst in circ.data:
        out.append(inst.operation, [loc[circ.find_bit(q).index] for q in inst.qubits])
    return out


def _touched(circ):
    return sorted({circ.find_bit(q).index for inst in circ.data for q in inst.qubits})


def pair_piece_sections(qc, ro, level="quick", seed=7):
    """Section checks of a controlled projected pair piece (gates pre, mid, post): the flag sections on (u, v)
    codes with a random state on the borrowed qubits, and the mid section on the reference vectors."""
    defs = {inst.operation.name: inst.operation.definition for inst in qc.data}
    assert set(defs) == {"pre", "mid", "post"}, set(defs)
    rng = np.random.default_rng(seed)
    regu, regv, (zu, zv) = ro["register_u"], ro["register_v"], ro["flag_qubits"]
    res = {}
    for name, flags_in in (("pre", False), ("post", True)):
        sec = defs[name]
        borrowed = [q for q in _touched(sec) if q not in regu + regv + [zu, zv]]
        order = regu + regv + [zu, zv] + borrowed
        loc = _remap(sec, order)
        nloc, nb = len(order), len(borrowed)
        b = rng.normal(size=2 ** nb) + 1j * rng.normal(size=2 ** nb)
        b /= np.linalg.norm(b)
        pairs = [(u, v) for u in range(64) for v in range(64)]
        if level == "quick":
            pairs = [pairs[i] for i in rng.choice(len(pairs), 64, replace=False)]
        dev = 0.0
        ph = np.exp(1j * ro["phases"][name])
        for i0 in range(0, len(pairs), 64):
            circs, outs_idx = [], []
            for (u, v) in pairs[i0:i0 + 64]:
                fu, fv = int(u != 0), int(v != 0)
                bin_ = u + (v << 6) + ((fu if flags_in else 0) << 12) + ((fv if flags_in else 0) << 13)
                bout = u + (v << 6) + ((0 if flags_in else fu) << 12) + ((0 if flags_in else fv) << 13)
                psi = np.zeros(2 ** nloc, dtype=complex)
                psi[bin_ + (np.arange(2 ** nb) << 14)] = b
                c = QuantumCircuit(nloc); c.set_statevector(psi); c.compose(loc, inplace=True)
                c.save_amplitudes([int(x) for x in bout + (np.arange(2 ** nb) << 14)], label="a")
                circs.append(c)
            r = sim().run(circs, shots=1).result()
            for i in range(len(circs)):
                a = np.asarray(r.data(i)["a"]) * ph
                dev = max(dev, float(np.abs(a - b).max()))
        res[f"{name}: {len(pairs)} (u, v) codes"] = dev
    refs = resolve(ro["reference"])
    local, cl = refs["local_qubits"], refs["control_local"]
    mid = _remap(defs["mid"], local)
    ph = np.exp(1j * ro["phases"]["mid"])
    reps = refs["representatives"]
    if level == "quick":
        seen, sel = set(), []
        for r in reps:
            if r["cls"] not in seen:
                seen.add(r["cls"]); sel.append(r)
        reps = sel
    for r in reps:
        codes = np.array(r["codes"])
        ins = r["inputs"][:1] if level == "quick" else r["inputs"]
        devs = []
        for vin, vout in zip(ins, r["expected"]):
            vin, vout = as_complex_array(vin), as_complex_array(vout)
            psi = np.zeros(2 ** len(local), dtype=complex)
            psi[codes] = vin / np.sqrt(2)
            psi[codes + (1 << cl)] = vin / np.sqrt(2)
            out = _run_flat(mid, psi) * ph * np.sqrt(2)
            devs.append(float(max(np.abs(out[codes] - vin).max(), np.abs(out[codes + (1 << cl)] - vout).max())))
        res[f"mid: {r['cls']}, {r['representative']}"] = max(devs)
    return res


# ---------------------------------------------------------------------------------------------- evaluate
def evaluate(entry, level="full"):
    """Run the readout of a manifest entry.  Returns dict(values=..., reference=..., max_dev=..., detail=...)."""
    ro = entry["readout"]
    t = ro["type"]
    if t == "none":
        return dict(skipped=ro["reason"])
    qc = load(entry["file"])
    nq = qc.num_qubits
    if t == "hadamard":
        psi0 = None
        if "prep" in ro:
            pre = QuantumCircuit(nq)
            for q in ro["prep"].get("x", []):
                pre.x(q)
            for q in ro["prep"].get("h", []):
                pre.h(q)
            pre.compose(qc, inplace=True)
            qc = pre
        psi = run(qc, psi0)
        C = hadamard_overlap(psi, ro["ref"], ro["test"])
        Cref = as_complex(resolve(ro["reference"]))
        return dict(values=[C.real, C.imag], reference=[Cref.real, Cref.imag], max_dev=abs(C - Cref))
    if t == "chebyshev_chain":
        uh = load(ro["U_H_file"])
        mu, k = chebyshev_moments(qc, uh, ro["U_H_phase"], ro["ancilla_qubits"], ro["n_max"], ro["walk_gate"])
        ref = np.array(resolve(ro["reference"]), dtype=float)
        return dict(values=mu.tolist(), reference=ref.tolist(), max_dev=float(np.abs(mu - ref).max()), walks=k)
    if t == "operator_column":
        ph = np.exp(1j * ro["phase"])
        out = run(qc, basis_state(nq, ro["input"])) * ph
        col = out[ro["outputs"]]
        ref = as_complex_array(resolve(ro["reference"]))
        idev = 0.0
        if "identity_input" in ro:
            o0 = run(qc, basis_state(nq, ro["identity_input"])) * ph
            idev = abs(o0[ro["identity_input"]] - 1)
        return dict(values=f"column of {len(col)} amplitudes", max_dev=float(max(np.abs(col - ref).max(), idev)),
                    detail=dict(column_dev=float(np.abs(col - ref).max()), identity_dev=float(idev)))
    if t == "walk_power":
        ph = np.exp(1j * ro["phase"])
        psi = basis_state(nq, ro["ref"])
        amps = [1.0 + 0j]
        for _ in range(ro["k_max"]):
            psi = run(qc, psi) * ph
            amps.append(psi[ro["ref"]])
        ref = np.array(resolve(ro["reference"])[:ro["k_max"] + 1], dtype=complex)
        amps = np.array(amps)
        return dict(values=[[z.real, z.imag] for z in amps], reference=ref.real.tolist(),
                    max_dev=float(np.abs(amps - ref).max()))
    if t == "amplitudes":
        ph = np.exp(1j * ro["phase"])
        vals, refs, dev = [], [], 0.0
        for it in ro["items"]:
            a = run(qc, basis_state(nq, it["input"]))[it["output"]] * ph
            r = as_complex(it["expect"]) if "expect" in it else as_complex(resolve(it["reference"]))
            vals.append([a.real, a.imag]); refs.append([r.real, r.imag])
            dev = max(dev, abs(a - r))
        return dict(values=vals, reference=refs, max_dev=float(dev))
    if t == "unit_operator":
        from qiskit.quantum_info import Operator
        U = Operator(qc).data
        E = unit_exact(ro)
        return dict(values="operator", max_dev=float(np.abs(U - E).max()))
    if t == "pair_piece_sections":
        res = pair_piece_sections(qc, ro, level)
        return dict(values="section checks", detail=res, max_dev=float(max(res.values())))
    raise ValueError(t)


def manifests():
    """All manifest entries of the repository: list of (folder, entry)."""
    out = []
    base = os.path.join(REPO, "circuits")
    for folder in sorted(os.listdir(base)):
        p = os.path.join(base, folder, "manifest.json")
        if os.path.exists(p):
            for e in json.load(open(p)):
                out.append((folder, e))
    return out
