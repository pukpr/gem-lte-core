#!/usr/bin/env python3
"""Earlier Oct2026_quad5 fit (seed_scores.out) vs the same fit after lpap_refine.py (run.out), for every finished
centre (DONE): training and held-out 2000-2005 CC per leg, whether the refinement was applied, dLOD.
Writes compare_own.csv and prints a digest.   python3 compare_own.py"""
import csv, re, statistics as st
from pathlib import Path
HERE = Path(__file__).resolve().parent
pat = re.compile(r"^(\S+): \S+ \| train ([+-][\d.]+).*?held out \S+ ([+-][\d.]+)", re.M)
rows = []
for d in sorted(p for p in HERE.iterdir() if p.is_dir() and (p / "DONE").exists()):
    a = pat.findall((d / "seed_scores.out").read_text()); b = pat.findall((d / "run.out").read_text())
    ref = (d / "refine.out").read_text() if (d / "refine.out").exists() else ""
    applied = "applied:" in ref and "not applied" not in ref
    dl = re.findall(r"dLOD ([\d.]+) -> ([\d.]+)", ref)
    for (n, t0, h0), (_, t1, h1) in zip(a, b):
        rows.append(dict(centre=d.name, leg=n, applied=int(applied), train0=float(t0), train1=float(t1), held0=float(h0), held1=float(h1),
                         dlod0=float(dl[-1][0]) if dl else "", dlod1=float(dl[-1][1]) if dl else ""))
if not rows: raise SystemExit("no finished centres yet")
with open(HERE / "compare_own.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
cs = sorted({r["centre"] for r in rows}); centre = [r for r in rows if r["leg"].startswith("C_")]
m = lambda k, rr: st.mean(r[k] for r in rr)
print(f"{len(cs)} centres finished, refinement applied in {sum(r['applied'] for r in centre)}")
for lab, rr in (("all legs", rows), ("centre legs", centre)):
    print(f"{lab:12s}: train {m('train0', rr):.3f} -> {m('train1', rr):.3f} | held-out 2000-05 {m('held0', rr):+.3f} -> {m('held1', rr):+.3f}"
          f" | held-out better in {sum(r['held1'] > r['held0'] for r in rr)}/{len(rr)}")
