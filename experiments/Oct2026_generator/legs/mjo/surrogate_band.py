#!/usr/bin/env python3
"""Chance level for the MJO leg's held-out skill in the 15-100 day band (2026-10-15).
The two-leg result's manifold and windings (mjo_soi/result.p) are kept; the MJO data are replaced by phase-randomized
copies (same spectrum, random timing) and the MJO leg is re-scored alone, no search: amplitudes regressed outside
2000-2005, then scored on 2000-2005, exactly as the real leg is."""
import os, re, shutil, subprocess, numpy as np
from concurrent.futures import ThreadPoolExecutor
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); sm = lambda v, k: np.convolve(v, np.ones(k) / k, mode="same"); band = lambda v: sm(v, 3) - sm(v, 21); rng = np.random.default_rng(13)
M = np.loadtxt("mjo_140E_5.dat"); R = re.sub(r"^-?(LEGS|SELECT_TOP|SELECT_TIE)\s.*\n", "", open("mjo_soi/lt.exe.resp").read(), flags=re.M)
ENV = dict(DLOD_REF="TRUE", F9="1", IDATE="1880.0", LOCKA="FALSE", LOCKT="FALSE", METRIC="CC", TREND="false", EXCLUDE="true", TRAIN_START="2000", TRAIN_END="2005", CLIMATE_INDEX="x.dat", TEST_ONLY="TRUE", TIMEOUT="60", NUMBER_OF_PROCESSORS="1")
def surr(v):
    Z = np.fft.rfft(v - v.mean()); ph = rng.uniform(0, 2 * np.pi, len(Z)); ph[0] = 0; return np.fft.irfft(np.abs(Z) * np.exp(1j * ph), len(v)) + v.mean()
def run(i, y):
    d = f"sw_{i}"; os.makedirs(d, exist_ok=True); shutil.copy("mjo_soi/result.p", f"{d}/lt.exe.p"); shutil.copy("dlod_ref.dat", d); open(f"{d}/lt.exe.resp", "w").write(R); np.savetxt(f"{d}/x.dat", np.column_stack([M[:, 0], y]), fmt=["%.6f", "%.4f"], delimiter="\t")
    subprocess.run("ulimit -s 65536; ../../lt.exe", shell=True, cwd=d, env={**os.environ, **ENV}, stdin=subprocess.DEVNULL, capture_output=True, timeout=300); a = np.loadtxt(f"{d}/lte_results.csv", delimiter=","); shutil.rmtree(d)
    t = a[:, 0]; g = (t >= 2000) & (t < 2005); return cc(a[g, 1], a[g, 2]), cc(band(a[:, 1])[g], band(a[:, 2])[g]), cc(band(a[:, 1])[~g], band(a[:, 2])[~g])
r0 = run(0, M[:, 1]); 
with ThreadPoolExecutor(8) as ex: S = np.array(list(ex.map(lambda z: run(*z), [(i, surr(M[:, 1])) for i in range(1, 201)])))
print(f"real MJO, two-leg manifold, re-scored alone: held out {r0[0]:+.3f}; 15-100 day band held out {r0[1]:+.3f} (train {r0[2]:+.3f})")
for k, nm in ((0, "held out, all time scales"), (1, "held out, 15-100 day band"), (2, "train, 15-100 day band")):
    s = S[:, k]; real = r0[k]; print(f"   200 surrogates, {nm:28s}: {s.mean():+.3f} +- {s.std():.3f}; real is {(real - s.mean()) / s.std():+.1f} sd, exceeded by {100 * np.mean(s >= real):.0f}% of surrogates")
