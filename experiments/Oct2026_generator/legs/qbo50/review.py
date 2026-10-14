#!/usr/bin/env python3
"""Review of the Feb2026/qbo50 fit (KIT monthly 50 hPa) with the current lt.exe, no search (2026-10-13).
 1. does the current binary reproduce the stored fit? (as is, and with TIDE_LEGACY TRUE)
 2. what was held out, and what the score is when the response amplitudes are trained on part of the record only
 3. how fast the windings turn against monthly sampling"""
import json, os, re, subprocess, numpy as np
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); ENV = dict(BZONE="TRUE", CLIMATE_INDEX="qbo50.dat", DLOD_REF="TRUE", F9="1", IDATE="1920.9", LOCKA="FALSE", LOCKT="FALSE", METRIC="CC", TREND="true", TEST_ONLY="TRUE", TIMEOUT="60", NUMBER_OF_PROCESSORS="2")
def run(extra):
    out = subprocess.run("ulimit -s 65536; ../lt.exe", shell=True, env={**os.environ, **ENV, **extra}, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300).stdout; dl = re.findall(r"([-\d.]+):dLOD:", out); return np.loadtxt("lte_results.csv", delimiter=","), (float(dl[-1]) if dl else float("nan"))
old = np.loadtxt("../../../Feb2026/qbo50/lte_results.csv", delimiter=","); print(f"stored fit (12 Sep): CC {cc(old[:, 1], old[:, 2]):+.4f} over {old[0, 0]:.1f}-{old[-1, 0]:.1f}; its settings held out 1953-1955 only")
for nm, ex in (("current binary, same settings", dict(EXCLUDE="true", TRAIN_START="1953", TRAIN_END="1955")), ("current binary, TIDE_LEGACY TRUE", dict(EXCLUDE="true", TRAIN_START="1953", TRAIN_END="1955", TIDE_LEGACY="TRUE"))):
    a, dl = run(ex); n = min(len(a), len(old)); print(f"1. {nm}: CC {cc(a[:, 1], a[:, 2]):+.4f}; model against the stored model CC {cc(a[:n, 1], old[:n, 1]):+.4f}; manifold against the stored manifold CC {cc(a[:n, 3], old[:n, 3]):+.4f}; dLOD {dl:.4f}")
    if nm.endswith("TRUE"): leg = cc(a[:, 1], a[:, 2])
best = dict(TIDE_LEGACY="TRUE") if leg > 0.7 else {}
print("2. response amplitudes trained on part of the record, manifold as stored (no search):")
for a0, a1 in ((2000, 2005), (1990, 2025), (1953, 1988), (1976, 2000)):
    a, _ = run(dict(EXCLUDE="true", TRAIN_START=str(a0), TRAIN_END=str(a1), **best)); t = a[:, 0]; g = (t >= a0) & (t < a1); print(f"   held out {a0}-{a1} ({int(g.sum())} months): train {cc(a[~g, 1], a[~g, 2]):+.3f}, held out {cc(a[g, 1], a[g, 2]):+.3f}")
a, _ = run(dict(EXCLUDE="true", TRAIN_START="1953", TRAIN_END="1955", **best)); t, F = a[:, 0], a[:, 3]; p = json.load(open("lt.exe.p")); k = abs(p["ltep"][-1]); H = sorted(int(h) for h in p["harm"]); dF = np.abs(np.diff(F))
print(f"3. base winding {k:.4f}, harmonics {H}; manifold range {F.min():.2f} to {F.max():.2f}; monthly change of the manifold: median {np.median(dF):.3f}, 90th percentile {np.percentile(dF, 90):.3f}")
for h in (1, H[len(H) // 2], H[-1]): c = h * k * dF; print(f"   harmonic {h:2d}: turns {h * k * (F.max() - F.min()):.1f} times over the manifold's range; per month median {np.median(c):.2f} cycles, above the half-cycle limit in {100 * np.mean(c > 0.5):.0f}% of months")
