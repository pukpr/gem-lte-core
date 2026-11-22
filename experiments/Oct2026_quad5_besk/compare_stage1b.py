#!/usr/bin/env python3
"""Stage 1b summary: per quad the training/held-out (2000-2005) of both seeds without search, which seed fits better by
basin, and after the 10-min search from the better seed vs the same quad in Oct2026_quad5.   python3 compare_stage1b.py"""
import json, re, statistics as st
from pathlib import Path
HERE = Path(__file__).resolve().parent; S1B = HERE / "s1b"; OLD = HERE.parent / "Oct2026_quad5"
sample = json.loads((HERE / "sample.json").read_text())
pat = re.compile(r"^(\S+): \S+ \| train ([+-][\d.]+).*?held out \S+ ([+-][\d.]+)", re.M)
def sc(f):
    s = pat.findall(f.read_text()) if f.exists() else []
    return (st.mean(float(a) for _, a, _ in s), st.mean(float(b) for _, _, b in s)) if s else (float("nan"), float("nan"))
rows = []
print(f"{'quad':11s} {'basin':8s} | {'atl seed: train held':>21s} | {'pac seed: train held':>21s} | better | searched: train held | Oct2026_quad5: train held")
for q, b in sorted(sample.items()):
    a, p = sc(S1B / f"{q}_atl/scores.out"), sc(S1B / f"{q}_pac/scores.out")
    s, o = sc(S1B / f"{q}_best/run.out"), sc(OLD / q / "run.out"); better = "atl" if a[0] > p[0] else "pac"
    rows.append((q, b, a, p, better, s, o))
    print(f"{q:11s} {b:8s} | {a[0]:.3f} {a[1]:+.3f}          | {p[0]:.3f} {p[1]:+.3f}          | {better}    | {s[0]:.3f} {s[1]:+.3f}        | {o[0]:.3f} {o[1]:+.3f}")
for b in ("atlantic", "pacific", "indian"):
    v = [r for r in rows if r[1] == b]
    print(f"{b:8s} {len(v):2d}: training prefers atl in {sum(r[4] == 'atl' for r in v)}, pac in {sum(r[4] == 'pac' for r in v)} | "
          f"held-out without search: atl {st.mean(r[2][1] for r in v):+.3f}, pac {st.mean(r[3][1] for r in v):+.3f}"
          + (f" | after search {st.mean(r[5][1] for r in v if r[5][1] == r[5][1]):+.3f} vs Oct2026_quad5 {st.mean(r[6][1] for r in v):+.3f}" if any(r[5][1] == r[5][1] for r in v) else ""))
