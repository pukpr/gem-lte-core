#!/usr/bin/env python3
"""Two controls on the Feb2026/qbo50 fit (TIDE_LEGACY TRUE, manifold as stored, no search):
 A. the half-record test on phase-randomized surrogates of the QBO (same spectrum, so the same ~28-month quasi-period
    and the same coherence time, but no relation to the manifold's timing): amplitudes trained outside 1990-2025
 B. when in the year does the wind change? size of month-to-month changes by calendar month, data against model"""
import os, re, shutil, subprocess, numpy as np
from concurrent.futures import ThreadPoolExecutor
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); rng = np.random.default_rng(5); D = np.loadtxt("qbo50.dat")
ENV = dict(BZONE="TRUE", DLOD_REF="TRUE", F9="1", IDATE="1920.9", LOCKA="FALSE", LOCKT="FALSE", METRIC="CC", TREND="true", TEST_ONLY="TRUE", TIMEOUT="60", NUMBER_OF_PROCESSORS="1", TIDE_LEGACY="TRUE", EXCLUDE="true", TRAIN_START="1990", TRAIN_END="2025")
def surr(v):
    Z = np.fft.rfft(v - v.mean()); ph = rng.uniform(0, 2 * np.pi, len(Z)); ph[0] = 0; return np.fft.irfft(np.abs(Z) * np.exp(1j * ph), len(v)) + v.mean()
def run(i, y):
    d = f"w{i % 6}_{i}"; os.makedirs(d, exist_ok=True)
    for f in ("lt.exe.p", "lt.exe.resp", "dlod_ref.dat"): shutil.copy(f, d)
    np.savetxt(f"{d}/x.dat", np.column_stack([D[:, 0], y]), fmt=["%.6f", "%.4f"], delimiter="\t"); os.path.exists(f"{d}/../lt.exe") or None
    subprocess.run("ulimit -s 65536; ../../lt.exe", shell=True, cwd=d, env={**os.environ, **ENV, "CLIMATE_INDEX": "x.dat"}, stdin=subprocess.DEVNULL, capture_output=True, timeout=300); a = np.loadtxt(f"{d}/lte_results.csv", delimiter=","); shutil.rmtree(d); t = a[:, 0]; g = t >= 1990; return cc(a[~g, 1], a[~g, 2]), cc(a[g, 1], a[g, 2]), a
tr, te, a = run(0, D[:, 1]); S = [surr(D[:, 1]) for _ in range(40)]
with ThreadPoolExecutor(6) as ex: R = np.array([r[:2] for r in ex.map(lambda z: run(*z), list(enumerate(S, 1)))])
print(f"A. amplitudes trained on 1953-1989, tested on 1990-2024: real QBO train {tr:+.3f}, held out {te:+.3f}")
print(f"   40 phase-randomized surrogates: train {R[:, 0].mean():+.3f} +- {R[:, 0].std():.3f} (largest {R[:, 0].max():+.3f}); held out {R[:, 1].mean():+.3f} +- {R[:, 1].std():.3f} (largest {R[:, 1].max():+.3f}); real held-out is {(te - R[:, 1].mean()) / R[:, 1].std():+.1f} sd above")
t, m, y, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]; mo = np.round((t[1:] % 1) * 12).astype(int) % 12; dm, dy, dF = np.abs(np.diff(m)), np.abs(np.diff(y)), np.abs(np.diff(F)); raw = np.abs(np.diff(D[:, 1])); mor = np.round((D[1:, 0] % 1) * 12).astype(int) % 12
print("B. mean size of the change INTO each calendar month (Jan .. Dec), each scaled to its own average:")
for nm, v, mm_ in (("manifold", dF, mo), ("model", dm, mo), ("data (filtered, as fitted)", dy, mo), ("data (raw monthly)", raw, mor)): print(f"   {nm:28s} " + " ".join(f"{v[mm_ == k].mean() / v.mean():5.2f}" for k in range(12)))
prof = lambda v, mm_: np.array([v[mm_ == k].mean() for k in range(12)]); print(f"   seasonal profile of the changes, model against raw data: CC {cc(prof(dm, mo), prof(raw, mor)):+.2f}; manifold against raw data: CC {cc(prof(dF, mo), prof(raw, mor)):+.2f}")
