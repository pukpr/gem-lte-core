#!/usr/bin/env python3
"""Which altitudes does the shared-manifold QBO fit carry through the excluded 2015-2021 anomalies? (2026-10-26)
Run directory: N QBO legs on one manifold (KIT monthly and/or ERA5 5-day), 2015-2021 held out (EXCLUDE); the legs are
read from lt.exe.resp.<k> (CLIMATE_INDEX), labelled KIT<hPa> / ERA<hPa> from the file names.
Reads lte_results.<k>.csv (decimal year, model, data) and reports per leg:
  train r (outside 2015-2021), held-out r and rms error / data sd over 2015-2021, and r in three sub-windows:
  2015.5-2017.0 (the 2016 disruption), 2017.0-2019.5 (between), 2019.5-2020.75 (the 2020 disruption).
usage: python3 excluded_region.py [label]      -> excluded_region.png, excluded_region.json"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
TAG = sys.argv[1] if len(sys.argv) > 1 else "result"
import re
LEGS = []
for k in range(1, 40):
    f = HERE / f"lt.exe.resp.{k}"
    if not f.exists():
        break
    name = f.read_text().split()[1]; lev = re.search(r"u(\d+)", name).group(1)
    LEGS.append((k, ("ERA" if name.startswith("era") else "KIT") + lev))
SUB = {"2016 disruption 2015.5-2017.0": (2015.5, 2017.0), "between 2017.0-2019.5": (2017.0, 2019.5),
       "2020 disruption 2019.5-2020.75": (2019.5, 2020.75)}


def r(a, b):
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 3 else float("nan")


out = {}
fig, ax = plt.subplots(len(LEGS), 1, figsize=(13, 2.3 * len(LEGS)), sharex=True)
for k, L in LEGS:
    a = np.loadtxt(HERE / f"lte_results.{k}.csv", delimiter=","); d, m, y = a[:, 0], a[:, 1], a[:, 2]
    ok = y != 0; h = (d >= 2015) & (d <= 2021)
    row = {"train_r": r(m[ok & ~h], y[ok & ~h]), "held_r": r(m[ok & h], y[ok & h]),
           "held_rmse_over_sd": float(np.sqrt(np.mean((m[ok & h] - y[ok & h]) ** 2)) / np.std(y[ok]))}
    for name, (s0, s1) in SUB.items():
        w = ok & (d >= s0) & (d < s1); row[name] = r(m[w], y[w])
    out[L] = row
    v = (d >= 2008) & (d <= 2024.5)
    ax[k - 1].plot(d[v], y[v], color="k", lw=1.2, label=L)
    ax[k - 1].plot(d[v], m[v], color="#d94a4a", lw=1.1, label="model (shared manifold)")
    ax[k - 1].axvspan(2015, 2021, color="#f2e6c9", alpha=0.7)
    for x in (2016.15, 2020.0): ax[k - 1].axvline(x, color="#2a78d6", ls=":", lw=1)
    ax[k - 1].set_ylabel("m/s"); ax[k - 1].legend(loc="upper left", fontsize=8, ncol=2)
    ax[k - 1].set_title(f"{L}: train r {row['train_r']:+.2f}, held-out 2015–21 r {row['held_r']:+.2f}, rms/sd {row['held_rmse_over_sd']:.2f}", fontsize=9.5)
ax[-1].set_xlabel("year (shaded: excluded from fitting; dotted: disruptions)")
fig.tight_layout(); fig.savefig(HERE / f"excluded_region_{TAG}.png", dpi=100)
json.dump(out, open(HERE / f"excluded_region_{TAG}.json", "w"), indent=1)
print(f"{'leg':>6s} {'train r':>8s} {'held r':>7s} {'rms/sd':>7s} | " + " | ".join(SUB))
for L, row in out.items():
    print(f"{L:>6s} {row['train_r']:+8.3f} {row['held_r']:+7.3f} {row['held_rmse_over_sd']:7.2f} | " + " | ".join(f"{row[s]:+.2f}".center(len(s)) for s in SUB))
