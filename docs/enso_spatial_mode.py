#!/usr/bin/env python3
"""Spatial structure of the fitted windings across the equatorial Pacific: standing or travelling? (2026-10-21)

Purpose: give the "rotational transform" a spatial meaning for ENSO. For each winding k = h * k_base (h = 1..9) and the
slow winding of the common manifold (AMO/PDO/NINO4/IODE four-leg fit), project each location's data on the winding,
model-free:  c(x) = sum_t y(t, x) exp(-2 pi i k F(t)),  with y detrended (quadratic), 1950.0-2022.5 common span.
Amplitude |c| and phase arg c along the equator show the spatial pattern of that winding:
  standing wave  -> one phase (or two, 180 deg apart, across a node), amplitude varying with longitude;
  travelling wave -> phase changing steadily with longitude.
Locations (SST): the 8 equatorial Kaplan quads 10S-10N from 130E to 90W (Oct2026_dlod) and the Nino regions
(4, 3.4, 3, 1+2). Pressure: Darwin and Tahiti (sign of pressure, not SST).
Significance: amplitude percentile among 200 phase-randomized surrogates of each location. Stability: phase in
1950-1986 vs 1987-2022.
    python3 enso_spatial_mode.py   -> enso_spatial_mode.png, table on stdout
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
NINO = EX / "Oct2026_generator/legs/nino4_34_3_12"
DT = EX / "Oct2026_generator/legs/darwin_tahiti"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
rng = np.random.default_rng(7)
T0, T1, MID = 1950.0, 2022.5, 1986.5

p = json.load(open(RUN / "lt.exe.p")); KB = abs(p["ltep"][-1]); KS = abs(p["ltep"][0])
man = np.loadtxt(RUN / "lte_results.1.csv", delimiter=","); TM, FM = man[:, 0], man[:, 3]
# (label, longitude east, file, kind)
LOCS = [(f"quad {q}", lon, EX / "Oct2026_dlod" / q / f"{q}.dat", "SST") for q, lon in
        (("kN000_E130", 130), ("kN000_E150", 150), ("kN000_E170", 170), ("kN000_W170", 190), ("kN000_W150", 210),
         ("kN000_W130", 230), ("kN000_W110", 250), ("kN000_W090", 270))]
LOCS += [("Nino4", 185, NINO / "nino4.dat", "SST"), ("Nino3.4", 215, NINO / "nino34.dat", "SST"),
         ("Nino3", 240, NINO / "nino3.dat", "SST"), ("Nino1+2", 275, NINO / "nino12.dat", "SST"),
         ("Darwin SLP", 131, DT / "darwin_anom.dat", "SLP"), ("Tahiti SLP", 210.4, DT / "tahiti_anom.dat", "SLP")]
WIND = [(f"{h}x", h * KB) for h in range(1, 10)] + [(f"slow ({KS / KB:.3f}x)", KS)]


def series(path):
    a = np.loadtxt(path); a = a[np.isfinite(a[:, 1]) & (a[:, 1] != 0) & (a[:, 0] >= T0) & (a[:, 0] < T1)]
    t, y = a[:, 0], a[:, 1]; return t, y - np.polyval(np.polyfit(t, y, 2), t)


def coef(t, y, k):
    return np.sum(y * np.exp(-2j * np.pi * k * np.interp(t, TM, FM))) / len(y)


def surr(y, n):
    Z = np.fft.rfft(y); ph = rng.uniform(0, 2 * np.pi, (n, len(Z))); ph[:, 0] = 0
    return np.fft.irfft(np.abs(Z) * np.exp(1j * ph), len(y), axis=1)


R = {}
for lab, lon, path, kind in LOCS:
    t, y = series(path); y = y / y.std(); S = surr(y, 200); f = np.interp(t, TM, FM)
    for wl, k in WIND:
        e = np.exp(-2j * np.pi * k * f); c = np.sum(y * e) / len(y); cs = (S @ e) / len(y)
        h1, h2 = t < MID, t >= MID
        R[(lab, wl)] = dict(lon=lon, kind=kind, amp=abs(c), phase=np.degrees(np.angle(c)), pct=float(np.mean(np.abs(cs) < abs(c)) * 100),
                            ph1=np.degrees(np.angle(np.sum(y[h1] * e[h1]))), ph2=np.degrees(np.angle(np.sum(y[h2] * e[h2]))))
wrap = lambda d: (d + 180) % 360 - 180  # noqa: E731
print(f"common manifold base {KB:.4f}; amplitude in units of each record's sd x 1e-2; pct = percentile among surrogates; "
      f"phase relative to Nino3.4 (deg); dphase = 1987-2022 minus 1950-1986")
for wl, k in WIND:
    ref = R[("Nino3.4", wl)]["phase"]
    sst = [(lab, R[(lab, wl)]) for lab, *_ in LOCS if R[(lab, wl)]["kind"] == "SST"]
    sig = np.mean([r["pct"] > 95 for _, r in sst])
    print(f"\n{wl:16s} SST locations above the 95th pct: {sig:.0%}")
    for lab, *_ in LOCS:
        r = R[(lab, wl)]
        print(f"   {lab:18s} lon {r['lon']:5.0f}  amp {100 * r['amp']:5.2f}  pct {r['pct']:5.1f}  phase {wrap(r['phase'] - ref):+5.0f}  dphase {wrap(r['ph2'] - r['ph1']):+5.0f}")
json.dump({f"{a}|{b}": v for (a, b), v in R.items()}, open(HERE / "enso_spatial_mode.json", "w"), indent=1)

fig, axs = plt.subplots(2, 1, figsize=(14, 8.5), facecolor=SURF, sharex=True)
cols = {"9x": "#2a78d6", "3x": "#eb6834", "6x": "#1baf7a", "1x": "#eda100", f"slow ({KS / KB:.3f}x)": "#e87ba4"}
for a in axs:
    a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.tick_params(colors=MUTED)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
for wl, col in cols.items():
    pts = sorted([(R[(lab, wl)]["lon"], R[(lab, wl)]) for lab, *_ in LOCS if lab.startswith("quad")])
    lon = [x for x, _ in pts]; ref = R[("Nino3.4", wl)]["phase"]
    axs[0].plot(lon, [100 * r["amp"] for _, r in pts], color=col, marker="o", lw=1.5, label=wl)
    axs[1].plot(lon, [wrap(r["phase"] - ref) for _, r in pts], color=col, marker="o", lw=1.5, label=wl)
    for lab in ("Nino4", "Nino3.4", "Nino3", "Nino1+2"):
        r = R[(lab, wl)]; axs[0].scatter(r["lon"], 100 * r["amp"], color=col, marker="s", s=40, edgecolor=INK, lw=0.5)
        axs[1].scatter(r["lon"], wrap(r["phase"] - ref), color=col, marker="s", s=40, edgecolor=INK, lw=0.5)
axs[0].set_ylabel("amplitude (sd × 10⁻²)", color=MUTED)
axs[0].set_title("Equatorial Pacific SST: amplitude of each winding of the common manifold (circles: 20° quads; squares: Niño regions)",
                 loc="left", fontsize=11, color=INK)
axs[0].legend(frameon=False, fontsize=9, ncol=5, loc="upper left")
axs[1].set_ylabel("phase relative to Niño3.4 (°)", color=MUTED); axs[1].set_ylim(-190, 190); axs[1].set_yticks([-180, -90, 0, 90, 180])
axs[1].set_title("Phase along the equator: flat = standing wave (sign flips at nodes), sloping = travelling wave", loc="left", fontsize=11, color=INK)
axs[1].set_xticks([130, 150, 170, 190, 210, 230, 250, 270])
axs[1].set_xticklabels(["130E", "150E", "170E", "170W", "150W", "130W", "110W", "90W"], color=MUTED)
fig.tight_layout(); fig.savefig(HERE / "enso_spatial_mode.png", dpi=100, facecolor=SURF)
print("\nplot: enso_spatial_mode.png")
