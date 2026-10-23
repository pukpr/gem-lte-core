#!/usr/bin/env python3
"""Tests 1 and 2 of landreman_lte_synthesis.md section 7, model-free (2026-10-21).

Manifold: the common manifold F(t) of the AMO/PDO/NINO4/IODE four-leg fit (legs/amo_pdo_nino4_iode, lte_results.1.csv
column 4), base winding k_b = 0.2017 (ltep[-1] of that fit), slow winding 0.1173 (= 0.581 k_b).
Records: data NOT used to fit that manifold: the 89 Kaplan SST quads (Oct2026_dlod), other indices (Feb2026), the 60
tide-gauge series (Oct2026_msl5/data). Each record: dates within the manifold's span, missing months dropped,
quadratic trend removed. "Independent" subset: records whose detrended data correlate |r| <= 0.5 with each of the four
fitting series (AMO, PDO, NINO4, IODE).
Spectrum: P(M) = |sum y(t) exp(-2 pi i M F(t))|^2 for M in [0.03, 2.0] (step 0.0025), the data projected on the
windings of the manifold, no harmonic ladder imposed.
Null: 200 phase-randomized surrogates of each record (same power spectrum), projected on the same manifold, so any
structure the manifold itself imposes appears in the null too.
 Test 2 (rational iota in the data): S = power within +-0.01 of h*k_b, h = 1..9, as a fraction of the power in
   M in [0.1, 2.0]; each record's percentile among its surrogates; mean percentile over records against 50.
 Test 1 (common slow ratio): each record's strongest winding in M in [0.03, 0.19] (below the base); the spread of those
   winding values across records compared with records drawn from the surrogates (one surrogate per record, 200 draws):
   fraction of records within +-0.004 of 0.1173 (ratio 0.581) and the tightest cluster (highest kernel density).
    python3 landreman_tests.py      -> landreman_tests.json, landreman_tests.png, digest on stdout
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EX = HERE.parent / "experiments"
RUN = EX / "Oct2026_generator/legs/amo_pdo_nino4_iode"
rng = np.random.default_rng(2026)
NS = 200
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

p = json.load(open(RUN / "lt.exe.p")); KB = abs(p["ltep"][-1]); KS = abs(p["ltep"][0])
man = np.loadtxt(RUN / "lte_results.1.csv", delimiter=","); TM, FM = man[:, 0], man[:, 3]
M = np.arange(0.03, 2.0001, 0.0025)
fit_series = {}
for k, n in enumerate(("AMO", "PDO", "NINO4", "IODE"), 1):
    a = np.loadtxt(RUN / f"lte_results.{k}.csv", delimiter=","); fit_series[n] = {round(t, 3): v for t, v, ok in zip(a[:, 0], a[:, 2], a[:, 2] != 0) if ok}


def records():
    out = {}
    for d in sorted((EX / "Oct2026_dlod").glob("k[NS]*_[EW]*")):
        out[f"quad {d.name}"] = d / f"{d.name}.dat"
    for i in ("pna", "sam", "nao", "baltic", "soi", "tna", "tpi", "emi", "ao", "darwin1880", "iod", "iodw", "nino34", "k3_amo_full"):
        out[f"index {i}"] = EX / "Feb2026" / i / f"{i}.dat"
    for f in sorted((EX / "Oct2026_msl5/data").glob("msl*.dat")):
        out[f"tide gauge {f.stem[3:]}"] = f
    return out


def prep(path):
    a = np.loadtxt(path); a = a[np.isfinite(a[:, 1]) & (a[:, 1] != 0)]
    a = a[(a[:, 0] >= TM[0]) & (a[:, 0] <= TM[-1])]
    if len(a) < 240:
        return None
    t, y = a[:, 0], a[:, 1]; y = y - np.polyval(np.polyfit(t, y, 2), t)
    return t, y


def surrogates(y, n):
    Z = np.fft.rfft(y); ph = rng.uniform(0, 2 * np.pi, (n, len(Z))); ph[:, 0] = 0
    return np.fft.irfft(np.abs(Z)[None, :] * np.exp(1j * ph), len(y), axis=1)


win_int = np.zeros_like(M, bool)
for h in range(1, 10):
    win_int |= np.abs(M - h * KB) <= 0.01
band = (M >= 0.1) & (M <= 2.0); sub = (M >= 0.03) & (M <= 0.19)

res = {}
for name, path in records().items():
    r = prep(path)
    if r is None:
        continue
    t, y = r
    key = [round(x, 3) for x in t]
    indep = True
    for s in fit_series.values():
        c = [(y[i], s[k]) for i, k in enumerate(key) if k in s]
        if len(c) > 120:
            cy = np.array(c); rr = np.corrcoef(cy[:, 0], cy[:, 1] - np.polyval(np.polyfit(np.arange(len(cy)), cy[:, 1], 2), np.arange(len(cy))))[0, 1]
            if abs(rr) > 0.5:
                indep = False
    E = np.exp(-2j * np.pi * M[:, None] * np.interp(t, TM, FM)[None, :])
    Y = np.vstack([y, surrogates(y, NS)])                       # row 0 = data
    P = np.abs(E @ Y.T) ** 2                                    # (nM, 1+NS)
    S = P[win_int & band].sum(axis=0) / P[band].sum(axis=0)
    slow = M[sub][np.argmax(P[sub], axis=0)]
    res[name] = dict(indep=indep, n=len(t), S=float(S[0]), S_pct=float(np.mean(S[1:] < S[0]) * 100),
                     slow=float(slow[0]), slow_surr=slow[1:].tolist(),
                     top=[float(M[band][i]) for i in np.argsort(-P[band, 0])[:5]])

names = list(res); ind = [n for n in names if res[n]["indep"]]
out = {"manifold_run": str(RUN.relative_to(EX.parent)), "k_base": KB, "k_slow": KS, "records": res}
(HERE / "landreman_tests.json").write_text(json.dumps(out, indent=1))


def digest(sel, label):
    pct = np.array([res[n]["S_pct"] for n in sel])
    se = 100 / np.sqrt(12 * len(sel))                           # sd of a mean of uniform percentiles
    print(f"\n[{label}: {len(sel)} records]")
    print(f" Test 2  power at integer multiples of the base: mean percentile vs surrogates {pct.mean():.1f} (null 50 +- {se:.1f}; z {(pct.mean() - 50) / se:+.2f}); "
          f"records above the 95th percentile {np.sum(pct > 95)} (expected {0.05 * len(sel):.1f})")
    real = np.array([res[n]["slow"] for n in sel]); surr = np.array([res[n]["slow_surr"] for n in sel])     # (nrec, NS)
    near = lambda v: np.mean(np.abs(v - KS) <= 0.004)
    nulln = np.array([near(surr[:, j]) for j in range(NS)])
    grid = np.arange(0.03, 0.19, 0.001)
    def dens(v):
        return np.max(np.exp(-0.5 * ((grid[:, None] - v[None, :]) / 0.003) ** 2).sum(axis=1)) / len(v)
    dr = dens(real); dn = np.array([dens(surr[:, j]) for j in range(NS)])
    print(f" Test 1  strongest sub-base winding: share within +-0.004 of {KS:.4f} (ratio {KS / KB:.3f}) {near(real):.2f} "
          f"(surrogates {nulln.mean():.2f}, 95th pct {np.percentile(nulln, 95):.2f}; p {np.mean(nulln >= near(real)):.3f}); "
          f"tightest cluster density {dr:.2f} at ratio {grid[np.argmax(np.exp(-0.5 * ((grid[:, None] - real[None, :]) / 0.003) ** 2).sum(axis=1))] / KB:.3f} "
          f"(surrogates 95th pct {np.percentile(dn, 95):.2f}; p {np.mean(dn >= dr):.3f})")
    return pct, real, surr


print(f"manifold {out['manifold_run']}, k_base {KB:.4f}, slow {KS:.4f} (ratio {KS / KB:.3f}); {NS} surrogates per record")
pa, ra, sa = digest(names, "all records")
pi, ri, si = digest(ind, "independent records (|r| <= 0.5 with AMO, PDO, NINO4, IODE)")
for grp in ("quad", "index", "tide gauge"):
    sel = [n for n in ind if n.startswith(grp)]
    if sel:
        digest(sel, f"independent {grp}s")

fig, axs = plt.subplots(1, 2, figsize=(15, 5), facecolor=SURF)
for a in axs:
    a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.tick_params(colors=MUTED)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
axs[0].hist(pi, bins=np.arange(0, 101, 5), color="#2a78d6", edgecolor=SURF)
axs[0].axhline(len(pi) / 20, color=INK, lw=1, ls="--", label="expected if no preference (uniform)")
axs[0].set_xlabel("percentile of the data among its surrogates", color=MUTED); axs[0].set_ylabel("records", color=MUTED)
axs[0].set_title(f"Test 2: power at integer multiples of the base winding\n({len(pi)} independent records)", loc="left", fontsize=11, color=INK)
axs[0].legend(frameon=False, fontsize=9)
bins = np.arange(0.03, 0.1901, 0.004) / KB
axs[1].hist(ri / KB, bins=bins, color="#2a78d6", edgecolor=SURF, label="data", density=True)
axs[1].hist(si.ravel() / KB, bins=bins, histtype="step", color="#eb6834", lw=1.5, label="surrogates", density=True)
axs[1].axvline(KS / KB, color=INK, lw=1, ls="--", label=f"slow winding of the fit ({KS / KB:.3f})")
axs[1].set_xlabel("strongest sub-base winding / base winding", color=MUTED); axs[1].set_ylabel("density", color=MUTED)
axs[1].set_title("Test 1: where each record's strongest slow winding falls", loc="left", fontsize=11, color=INK)
axs[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(HERE / "landreman_tests.png", dpi=100, facecolor=SURF)
print("\nplot: landreman_tests.png; per-record results: landreman_tests.json")
