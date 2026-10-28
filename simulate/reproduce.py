"""Run one circuit file and print its value(s) next to the reference value in data/.

    python3 simulate/reproduce.py circuits/he4_L3_sym/gqsp_hadamard_n03.qasm
    python3 simulate/reproduce.py circuits/o16_L7_sym/cheb_walk_chain_k40.qasm.gz

The file is looked up in the manifests (circuits/*/manifest.json); its "readout" says what to compute and its
"readout.reference" where the reference value is stored.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuits as CI  # noqa: E402


def find(path):
    p = os.path.relpath(os.path.abspath(path), CI.REPO) if os.path.exists(path) else path
    for folder, e in CI.manifests():
        if e["file"] == p or e["file"] == p + ".gz" or e["file"].replace(".gz", "") == p:
            return e
    raise SystemExit(f"{path}: not in any manifest")


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    e = find(sys.argv[1])
    level = "full" if "--full" in sys.argv else "quick"
    t0 = time.time()
    qc = CI.load(e["file"])
    cx = CI.count_cx(qc)
    print(f"{e['file']}: {qc.num_qubits} qubits, {cx} CNOT (manifest {e['cx']}"
          + (f", paper {e['paper_cx']}" if "paper_cx" in e else "") + ")")
    print(f"  {e.get('description', '')}")
    r = CI.evaluate(e, level)
    if "skipped" in r:
        print(f"  not simulated: {r['skipped']}")
        return
    ro = e["readout"]
    if ro["type"] == "hadamard":
        print(f"  C = 2 conj(A0) A1 = {complex(*r['values']):.12f}")
        print(f"  reference         = {complex(*r['reference']):.12f}   ({ro['reference']['file']})")
    elif ro["type"] == "chebyshev_chain":
        mu, ref = np.array(r["values"]), np.array(r["reference"])
        print(f"  {r['walks']} walks; moments mu_n, n = 0..{len(mu)-1} (reference {ro['reference']['file']}):")
        for n in range(len(mu)):
            print(f"    n = {n:2d}: {mu[n]:+.12f}   reference {ref[n]:+.12f}")
    elif ro["type"] in ("amplitudes", "walk_power"):
        refs = r["reference"]
        for i, v in enumerate(r["values"]):
            rv = refs[i]
            print(f"  {complex(*v) if isinstance(v, list) else v:.12f}   reference "
                  f"{complex(*rv) if isinstance(rv, list) else rv:.12f}")
    elif "detail" in r:
        for k, v in r["detail"].items():
            print(f"  {k}: deviation {v:.2e}")
    print(f"  largest deviation from the reference: {r['max_dev']:.2e}   ({time.time()-t0:.1f} s)")


if __name__ == "__main__":
    main()
