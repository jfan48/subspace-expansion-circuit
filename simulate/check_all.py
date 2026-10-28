"""Check every circuit of the repository against its manifest and reference data.

    python3 simulate/check_all.py --quick      CNOT counts of all files + a representative subset simulated (< 10 min)
    python3 simulate/check_all.py --full       CNOT counts of all files + every file that fits a statevector
    python3 simulate/check_all.py --full --cx-only     CNOT counts only

CNOT: the count of the loaded file (recursive over its gate definitions) must equal the manifest's "cx", the
paper's count where the manifest gives one ("paper_cx"; "paper_cx_added" = increase over the file named in
"paper_cx_added_relative_to"; "paper_cx_approx" = the paper's three-digit value).
Values: the manifest entry's "readout" is simulated (Aer statevector, double precision) and compared with the
reference in data/; a deviation above --tol (default 1e-9) is a failure.  Exit status 1 on any failure.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuits as CI  # noqa: E402

QUICK_GQSP = {"gqsp_hadamard_n01", "gqsp_hadamard_n02", "gqsp_hadamard_n03", "gqsp_hadamard_T0.05"}


def stem(e):
    return os.path.basename(e["file"]).replace(".gz", "").replace(".qasm", "")


def plan(folder, e, level):
    """None = simulate; otherwise the reason for not simulating."""
    ro = e["readout"]
    if ro["type"] == "none":
        return "not simulable: " + ro["reason"]
    if e.get("system") == "triton" and e["family"].startswith("Hadamard") and e["n"] > (2 if level == "quick" else 3):
        return ("triton n > 3: 25 qubits, n controlled S4 steps per file (about 30 min for the 70 files); the files "
                "are the same step definition called n times, n <= 3 checked" if level == "full" else "quick subset")
    if level == "full":
        return None
    s = stem(e)
    if e["family"] == "Table XI step" and e["qubits"] > 20:
        return "quick subset (25-qubit Table XI step files run in --full)"
    if e.get("system") == "deuteron" and e["family"].startswith("Hadamard"):
        return None if e["n"] <= 3 else "quick subset"
    if s.startswith("gqsp_hadamard"):
        return None if s in QUICK_GQSP else "quick subset"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--cx-only", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--only", default=None, help="substring of the file path")
    args = ap.parse_args()
    level = "full" if args.full else "quick"
    T0 = time.time()
    entries = CI.manifests()
    if args.only:
        entries = [(f, e) for f, e in entries if args.only in e["file"]]
    fails = []
    # ------------------------------------------------------------------ CNOT counts
    cx = {}
    print(f"CNOT counts of {len(entries)} files", flush=True)
    for folder, e in entries:
        qc = CI.load(e["file"])
        n = CI.count_cx(qc)
        cx[e["file"]] = (n, qc.num_qubits)
    fam_rows = {}
    for folder, e in entries:
        n, nq = cx[e["file"]]
        ok = n == e["cx"] and nq == e["qubits"]
        paper = ""
        if "paper_cx" in e:
            ok &= n == e["paper_cx"]
            paper = str(e["paper_cx"])
        if "paper_cx_added" in e:
            rel = cx.get(e["paper_cx_added_relative_to"], (None,))[0]
            if rel is not None:
                ok &= n - rel == e["paper_cx_added"]
            paper = f"+{e['paper_cx_added']}"
        if "paper_cx_approx" in e:
            ok &= float(f"{n:.3g}") == e["paper_cx_approx"]
            paper = f"{e['paper_cx_approx']:.3g}"
        if not ok:
            fails.append((e["file"], f"CNOT {n} / manifest {e['cx']} / paper {paper or '-'}; qubits {nq} / {e['qubits']}"))
        key = (folder, e["family"])
        r = fam_rows.setdefault(key, dict(n=0, ok=0, cx=[], nq=set(), paper=[]))
        r["n"] += 1; r["ok"] += ok; r["cx"].append(n); r["nq"].add(nq)
        if paper:
            r["paper"].append(f"{os.path.basename(e['file']).split('.qasm')[0]}={paper}")
    print(f"{'folder':16s} {'family':38s} {'files':>5s} {'qubits':>9s} {'CNOT range':>21s}  ok   paper counts checked")
    for (folder, fam), r in fam_rows.items():
        q = ",".join(str(x) for x in sorted(r["nq"]))
        print(f"{folder:16s} {fam[:38]:38s} {r['n']:5d} {q:>9s} {min(r['cx']):>10d}..{max(r['cx']):<10d} "
              f"{r['ok']:3d}/{r['n']:<3d} {'; '.join(r['paper'])[:110]}")
    print(f"CNOT pass: {time.time()-T0:.0f} s", flush=True)
    if args.cx_only:
        _finish(fails, T0)
        return
    # ------------------------------------------------------------------ values
    print(f"\nvalues ({level}); tolerance {args.tol:g}", flush=True)
    print(f"{'file':58s} {'qubits':>6s} {'CNOT':>9s} {'readout':18s} {'max dev':>9s} {'time':>7s}  status")
    skipped = {}
    for folder, e in entries:
        why = plan(folder, e, level)
        if why is not None:
            skipped.setdefault(why, []).append(e["file"])
            continue
        t1 = time.time()
        try:
            r = CI.evaluate(e, level)
            dev = r["max_dev"]
            ok = dev <= args.tol
            st = "ok" if ok else "FAIL"
            if not ok:
                fails.append((e["file"], f"deviation {dev:.2e} > {args.tol:g}"))
        except Exception as exc:                     # noqa: BLE001
            dev, st = float("nan"), "ERROR"
            fails.append((e["file"], f"{type(exc).__name__}: {exc}"))
        print(f"{e['file'][9:][:58]:58s} {cx[e['file']][1]:6d} {cx[e['file']][0]:9d} {e['readout']['type'][:18]:18s} "
              f"{dev:9.2e} {time.time()-t1:6.1f}s  {st}", flush=True)
    for why, fl in skipped.items():
        print(f"not simulated ({len(fl)} files): {why}")
    _finish(fails, T0)


def _finish(fails, T0):
    print(f"\ntotal time {time.time()-T0:.0f} s; failures: {len(fails)}")
    for f, msg in fails:
        print(f"  FAIL {f}: {msg}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
