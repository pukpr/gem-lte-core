#!/usr/bin/env python3
"""Pre-Bessel forcing manifold (lte_results column 4 with MANIFOLD_PRE) of every quad, before and after lpap_refine,
against the 6-leg [nao, nino4, amo, pdo, iode, tna] manifold (legs/...leg6.tna.extrapolate), on the quads' dates.
  before  = <quad>/refine_work/baseline/lte_results.1.csv  (the Oct2026_quad5 fit, scored with MANIFOLD_PRE)
  after   = the trust-region round whose lpap is refined_lpap.json
Per quad: CC, CC of month-to-month differences (offset and trend removed), mean offset, std ratio.
Writes compare_manifolds.csv and compare_manifolds.png.   python3 compare_manifolds.py"""
import csv, glob, json, re, statistics as st
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
REF = HERE.parent / "Oct2026_generator/legs/amo_pdo_nao_nino4_feb_overnight.wow.recover.leg6.tna.extrapolate"
r = np.loadtxt(REF / "lte_results.1.csv", delimiter=","); ref = {round(t, 4): v for t, v in zip(r[:, 0], r[:, 3])}
def col4(f):
    a = np.loadtxt(f, delimiter=","); return a[:, 0], a[:, 3]
def stats(t, F):
    m = np.array([round(x, 4) in ref for x in t]); R = np.array([ref[round(x, 4)] for x in t[m]]); Q = F[m]
    return dict(cc=float(np.corrcoef(Q, R)[0, 1]), cc_diff=float(np.corrcoef(np.diff(Q), np.diff(R))[0, 1]),
                offset=float(np.mean(Q - R)), std_ratio=float(Q.std() / R.std())), t[m], Q, R
rows, curves = [], {}
for d in sorted(p for p in HERE.iterdir() if p.is_dir() and re.match(r"^k[NS]\d{3}_[EW]\d{3}$", p.name) and (p / "DONE").exists()):
    W = d / "refine_work"; lp1 = json.load(open(W / "refined_lpap.json"))["lpap"]
    best = None
    for rd in sorted(W.glob("round*")):
        if np.allclose(np.array(json.load(open(rd / "lt.exe.p"))["lpap"])[:, 1:3], np.array(lp1)[:, 1:3], atol=1e-12): best = rd
    if best is None: print(f"{d.name}: no round matches refined_lpap.json, skipped"); continue
    s0, t, Q0, R = stats(*col4(W / "baseline/lte_results.1.csv")); s1, _, Q1, _ = stats(*col4(best / "lte_results.1.csv"))
    rows.append(dict(quad=d.name, **{f"before_{k}": v for k, v in s0.items()}, **{f"after_{k}": v for k, v in s1.items()}))
    curves[d.name] = (t, Q0, Q1, R)
with open(HERE / "compare_manifolds.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
m = lambda k: st.median(x[k] for x in rows)
print(f"{len(rows)} quads vs the 6-leg manifold (1950-2023), medians [min .. max]:")
for k in ("cc", "cc_diff", "offset", "std_ratio"):
    b, a = [x["before_" + k] for x in rows], [x["after_" + k] for x in rows]
    print(f"  {k:9s} before {st.median(b):+.4f} [{min(b):+.4f} .. {max(b):+.4f}]   after {st.median(a):+.4f} [{min(a):+.4f} .. {max(a):+.4f}]")
# the quads among themselves: how far apart are their manifolds?
Q0s = np.array([c[1] for c in curves.values()]); Q1s = np.array([c[2] for c in curves.values()])
for lab, Qs in (("before", Q0s), ("after", Q1s)):
    med = np.median(Qs, axis=0); dd = [np.corrcoef(np.diff(q), np.diff(med))[0, 1] for q in Qs]
    print(f"  quads vs their own median manifold, CC of differences ({lab}): median {st.median(dd):.4f}, min {min(dd):.4f}; spread across quads (rms of q - median) median {st.median([np.std(q - med) for q in Qs]):.3f}")
t, _, _, R = next(iter(curves.values()))
fig, ax = plt.subplots(2, 1, figsize=(14, 8), facecolor="#fcfcfb", sharex=True)
for A, Qs, lab in ((ax[0], Q0s, "before refinement (Oct2026_quad5 fits)"), (ax[1], Q1s, "after lpap_refine")):
    A.set_facecolor("#fcfcfb")
    for q in Qs: A.plot(t, q, color="#a3a29c", lw=0.4, alpha=0.5)
    A.plot(t, np.median(Qs, axis=0), color="#2a78d6", lw=1.4, label=f"median of {len(Qs)} quads")
    A.plot(t, R, color="#eb6834", lw=1.4, label="6-leg [nao, nino4, amo, pdo, iode, tna]")
    A.set_title(f"pre-Bessel manifold, {lab}", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, loc="upper right")
    A.grid(True, ls=":", lw=0.5, alpha=0.7)
    for s in ("top", "right"): A.spines[s].set_visible(False)
ax[1].set_xlabel("year"); ax[0].set_ylabel("F (column 4)"); ax[1].set_ylabel("F (column 4)")
fig.tight_layout(); fig.savefig(HERE / "compare_manifolds.png", dpi=110, facecolor="#fcfcfb"); print("plot:", HERE / "compare_manifolds.png")
