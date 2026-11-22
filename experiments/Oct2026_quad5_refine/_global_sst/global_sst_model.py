#!/usr/bin/env python3
"""Global monthly SST from the last quad sweep's model fits, against the consensus of observed records (2026-11-22).

Model: every Oct2026_quad5_refine quad's refined lt.exe.p, its centre leg rescored alone TEST_ONLY with F9=-1
(nofilt/<quad>, rescore.sh; with F9=1 the same rescore reproduces the sweep's leg-1 model to 1.5e-11), so columns 2/3
are the unfiltered monthly model and data. Both are pro-rated by each quad's ocean area on the Kaplan grid
(Sep2026/global_sst_from_quads: weights renormalized each month over the quads present).

Observed global mean SST, each as monthly anomalies from its own 1961-1990 calendar-month means:
  Kaplan    cos(lat)-weighted mean of the Kaplan SST V2 grid the quads were cut from (Feb2026/k_sst)
  HadSST    Met Office HadSST 4.2.0.0 monthly GLOBE (raw/hadsst4_globe.csv)
  NOAA      NOAA NCEI Climate at a Glance, global ocean (NOAAGlobalTemp, ERSST v5; raw/noaa_ocean.csv)
  consensus the mean of the three, over the months all three have
The model is put on the same footing (its own 1961-1990 calendar-month means removed); no scale is fitted.
Scores: CC and RMS against the consensus, raw monthly, after removing each series' own linear trend over the span
(the trend is the shared nuisance), and 12-month running means; training months and the 2000-2005 held-out gap apart.
    python3 global_sst_model.py   -> global_sst_model.png, global_sst_model.txt (monthly table)"""
import sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "Sep2026")); import global_sst_from_quads as G

quads = G.load_quads(HERE / "nofilt")
nc = G.load_kaplan(EXP / "Feb2026/k_sst/sst.mean.anom.nc")
weights, covered = G.ocean_area_weights(sorted(quads), nc)
keys = np.arange(min(G.month_key(q["dates"][0]) for q in quads.values()),
                 max(G.month_key(q["dates"][-1]) for q in quads.values()) + 1)
model, frac = G.prorate(quads, weights, "model", keys)
recon, _ = G.prorate(quads, weights, "data", keys)
kk, sst, lat, _ = nc; kap = dict(zip(kk, G.global_mean(sst, lat)))
had = {int(r[0]) * 12 + int(r[1]) - 1: r[2] for r in np.genfromtxt(HERE / "raw/hadsst4_globe.csv", delimiter=",", skip_header=1)}
noaa = {}
for l in open(HERE / "raw/noaa_ocean.csv"):
    x = l.strip().split(",")
    if len(x) == 2 and x[0].isdigit(): noaa[int(x[0][:4]) * 12 + int(x[0][4:]) - 1] = float(x[1])

def rebase(d):  # remove own 1961-1990 calendar-month means
    out = dict(d)
    for m in range(12):
        base = [v for k, v in d.items() if k % 12 == m and 1961 * 12 <= k < 1991 * 12 and np.isfinite(v)]
        for k in d:
            if k % 12 == m: out[k] = d[k] - np.mean(base)
    return out
S = {"model": rebase(dict(zip(keys, model))), "quad data": rebase(dict(zip(keys, recon))),
     "Kaplan": rebase(kap), "HadSST": rebase(had), "NOAA": rebase(noaa)}
obs = ("Kaplan", "HadSST", "NOAA")
common = [k for k in keys if all(k in S[o] and np.isfinite(S[o][k]) for o in obs) and np.isfinite(S["model"][k])]
S["consensus"] = {k: np.mean([S[o][k] for o in obs]) for k in common}
t = np.array(common) / 12.0; held = (t >= 2000) & (t < 2005)
V = {n: np.array([S[n][k] for k in common]) for n in S}

def detr(y): return y - np.polyval(np.polyfit(t, y, 1), t)
def run12(y): return np.convolve(y, np.ones(12) / 12, mode="same")
def score(a, b, m):
    return np.corrcoef(a[m], b[m])[0, 1], np.sqrt(np.mean((a[m] - b[m]) ** 2))
ref = V["consensus"]
lines = [f"{len(quads)} quads, {t[0]:.0f}-{t[-1]:.2f} ({len(t)} months common to model, Kaplan, HadSST and NOAA)",
         f"quad coverage: mean {np.nanmean(frac):.3f} of the quads' total ocean weight present",
         "", f"{'vs consensus':14s} {'monthly CC / RMS':>18s} {'detrended CC':>13s} {'12-mo CC':>9s} {'detr 12-mo':>11s} {'held 2000-05 CC':>16s}"]
for n in ("model", "quad data", "Kaplan", "HadSST", "NOAA"):
    a = V[n]; c, r = score(a, ref, ~held); cd, _ = score(detr(a), detr(ref), ~held)
    cs, _ = score(run12(a), run12(ref), ~held); cds, _ = score(run12(detr(a)), run12(detr(ref)), ~held)
    ch, _ = score(a, ref, held)
    lines.append(f"{n:14s} {c:+.3f} / {r:.3f} C      {cd:+.3f}        {cs:+.3f}     {cds:+.3f}        {ch:+.3f}")
for n in ("model", "consensus"):
    lines.append(f"linear trend {n:10s}: {np.polyfit(t, V[n], 1)[0] * 100:+.2f} C/century")
lines.append(f"observed spread (mean |dataset - consensus|): {np.mean([np.abs(V[o] - ref).mean() for o in obs]):.3f} C")
print("\n".join(lines))
with open(HERE / "global_sst_model.txt", "w") as fh:
    fh.write("# " + "\n# ".join(lines) + "\n# date model consensus Kaplan HadSST NOAA quad_data\n")
    for i, k in enumerate(common): fh.write(f"{k / 12:.4f} " + " ".join(f"{V[n][i]:+.4f}" for n in ("model", "consensus", "Kaplan", "HadSST", "NOAA", "quad data")) + "\n")

INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
fig, ax = plt.subplots(3, 1, figsize=(15, 11), facecolor=SURF, gridspec_kw=dict(height_ratios=[3, 2, 2]))
for A in ax:
    A.set_facecolor(SURF); A.axvspan(2000, 2005, color=GRID, alpha=0.8, lw=0); A.grid(color=GRID, lw=0.6)
    for s in ("top", "right"): A.spines[s].set_visible(False)
A = ax[0]
A.plot(t, V["consensus"], color="#9a9893", lw=0.8, label="consensus, monthly")
A.plot(t, run12(V["consensus"]), color=INK, lw=2, label="consensus (Kaplan, HadSST, NOAA), 12-month mean")
A.plot(t, run12(V["model"]), color="#eb6834", lw=2, label="model, 12-month mean")
A.set_title("Global mean SST: quad model (89 refined Oct2026_quad5 fits, area-weighted) vs observed consensus", loc="left")
A.set_ylabel("anomaly from 1961-1990 (C)"); A.legend(frameon=False, ncol=3, fontsize=9)
A = ax[1]
A.plot(t, run12(detr(V["consensus"])), color=INK, lw=1.8, label="consensus, detrended")
A.plot(t, run12(detr(V["model"])), color="#eb6834", lw=1.8, label="model, detrended")
A.set_ylabel("detrended, 12-mo (C)"); A.legend(frameon=False, ncol=2, fontsize=9)
A = ax[2]
for o, c in zip(obs, ("#2a78d6", "#1baf7a", "#e87ba4")): A.plot(t, run12(V[o] - ref), color=c, lw=1.2, label=f"{o} - consensus")
A.plot(t, run12(V["model"] - ref), color="#eb6834", lw=1.8, label="model - consensus")
A.axhline(0, color=MUTED, lw=0.8); A.set_ylabel("difference, 12-mo (C)"); A.legend(frameon=False, ncol=4, fontsize=9)
for A in ax: A.set_xlim(t[0], t[-1])
fig.tight_layout(); fig.savefig(HERE / "global_sst_model.png", dpi=100, facecolor=SURF); print(f"plot: {HERE / 'global_sst_model.png'}")
