#!/usr/bin/env python3
"""Stage 1 (basin-seeded quads) vs the same quads in Oct2026_quad5: final besk by basin, training and held-out
2000-2005 (leg_scores2 lines in run.out), dominant winding, and the pre-Bessel manifold vs the 6-leg fit.
    python3 compare_stage1.py        (writes compare_stage1.csv)"""
import csv, json, re, statistics as st
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OLD = HERE.parent / "Oct2026_quad5"
SIX = HERE.parent / "Oct2026_generator/legs/amo_pdo_nao_nino4_feb_overnight.wow.recover.leg6.tna.extrapolate"
sample = json.loads((HERE / "sample.json").read_text())
pat = re.compile(r"^(\S+): \S+ \| train ([+-][\d.]+).*?held out \S+ ([+-][\d.]+)", re.M)
r6 = np.loadtxt(SIX / "lte_results.1.csv", delimiter=","); ref = {round(t, 4): v for t, v in zip(r6[:, 0], r6[:, 3])}
def scores(f):
    s = pat.findall(f.read_text()) if f.exists() else []
    return [(float(a), float(b)) for _, a, b in s]
rows = []
for q, basin in sorted(sample.items()):
    d = HERE / q
    if not (d / "legs.png").exists(): continue
    p = json.load(open(d / "result.p")); s0 = json.load(open(d / "seed.p"))
    new, old = scores(d / "run.out"), scores(OLD / q / "run.out")
    w = json.load(open(d / "lt.exe.windings.1.json")); kap = sorted(w["k_amp_phase"], key=lambda r: -r[1])
    a = np.loadtxt(d / "lte_results.1.csv", delimiter=","); m = np.array([round(t, 4) in ref for t in a[:, 0]])
    Q = a[m, 3]; R = np.array([ref[round(t, 4)] for t in a[m, 0]])
    rows.append(dict(quad=q, basin=basin, besk_seed=s0["besk"], besk=p["besk"], base_k=p["ltep"][-1], dom_k=kap[0][0], dom_amp=kap[0][1],
                     train_new=st.mean(x for x, _ in new), held_new=st.mean(y for _, y in new),
                     train_old=st.mean(x for x, _ in old), held_old=st.mean(y for _, y in old),
                     centre_held_new=new[0][1], centre_held_old=old[0][1],
                     manifold_ccdiff=float(np.corrcoef(np.diff(Q), np.diff(R))[0, 1]), manifold_offset=float(np.mean(Q - R)), manifold_std=float(Q.std() / R.std())))
if not rows: raise SystemExit("no finished quads yet")
with open(HERE / "compare_stage1.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
print(f"{len(rows)} of {len(sample)} quads finished")
print(f"{'quad':11s} {'basin':8s} {'besk seed->final':>17s} {'dominant k (amp)':>17s} | {'train old->new':>14s} {'held-out old->new':>18s} | manifold vs 6-leg: dCC, offset, std")
for r in rows:
    print(f"{r['quad']:11s} {r['basin']:8s} {r['besk_seed']:.4f}->{r['besk']:.4f}  {r['dom_k']:+.3f} ({r['dom_amp']:.3f})  | {r['train_old']:.3f}->{r['train_new']:.3f}  {r['held_old']:+.3f}->{r['held_new']:+.3f}   | {r['manifold_ccdiff']:.3f} {r['manifold_offset']:+.2f} {r['manifold_std']:.3f}")
for b in ("atlantic", "pacific", "indian"):
    v = [r for r in rows if r["basin"] == b]
    if v: print(f"{b:8s} {len(v):2d}: besk median {st.median(r['besk'] for r in v):.4f} [{min(r['besk'] for r in v):.4f}..{max(r['besk'] for r in v):.4f}] | "
                f"train {st.mean(r['train_old'] for r in v):.3f}->{st.mean(r['train_new'] for r in v):.3f} | held-out {st.mean(r['held_old'] for r in v):+.3f}->{st.mean(r['held_new'] for r in v):+.3f}"
                f" (better in {sum(r['held_new'] > r['held_old'] for r in v)}/{len(v)})")
