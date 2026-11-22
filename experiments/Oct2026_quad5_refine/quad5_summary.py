#!/usr/bin/env python3
"""Summary of the five-leg quad sweep (Oct2026_quad5): for every finished centre (legs.png present)
  - centre leg: train / held-out (2000-2005) correlation, also with a training-fitted quadratic trend removed
  - each neighbour leg: the same, its detrended data correlation with the centre, and its response agreement with the
    centre = complex correlation of the fitted winding terms (k_amp_phase), +1 identical response, -1 mirror image
  - reference: the centre's single-quad fit from Oct2026_dlod scored on the same held-out years
Writes quad5_summary.csv and quad5_maps.png (centre held-out skill; mean neighbour agreement), prints a digest.
    python3 quad5_summary.py
"""
import csv
import json
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DLOD = HERE.parent / "Oct2026_dlod"
H0, H1 = 2000.0, 2005.0
NAME = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1])  # noqa: E731


def centre(q):
    m = NAME.match(q)
    return (int(m.group(2)) * (1 if m.group(1) == "N" else -1), int(m.group(4)) * (1 if m.group(3) == "E" else -1))


def scores(a):
    t, m, y = a[:, 0], a[:, 1], a[:, 2]
    ok = y != 0; h = (t >= H0) & (t <= H1); tr = ok & ~h; ho = ok & h
    dq = lambda v: v - np.polyval(np.polyfit(t[tr], v[tr], 2), t)  # noqa: E731
    return cc(m[tr], y[tr]), cc(m[ho], y[ho]), cc(dq(m)[tr], dq(y)[tr]), cc(dq(m)[ho], dq(y)[ho]), dq(y)


def cvec(path):
    return np.array([x[1] * np.exp(1j * x[2]) for x in json.load(open(path))["k_amp_phase"]])


agree = lambda a, b: float(np.real(np.vdot(a, b)) / np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2)))  # noqa: E731
nb = json.load(open(HERE / "neighbours.json"))
rows = []
for q, legs in nb.items():
    d = HERE / q
    if not (d / "legs.png").exists():
        continue
    A = {k: np.loadtxt(d / f"lte_results.{k}.csv", delimiter=",") for k in range(1, 6)}
    C = {k: cvec(d / f"lt.exe.windings.{k}.json") for k in range(1, 6)}
    s1 = scores(A[1])
    ref = np.loadtxt(DLOD / q / "lte_results.csv", delimiter=",")
    rt, rm, ry = ref[:, 0], ref[:, 1], ref[:, 2]; rh = (rt >= H0) & (rt <= H1) & (ry != 0)
    row = dict(centre=q, lat=centre(q)[0], lon=centre(q)[1], train=s1[0], held=s1[1], train_detr=s1[2], held_detr=s1[3],
               single_fit_held=cc(rm[rh], ry[rh]))
    ag, dr, nh = [], [], []
    for k, n in zip(range(2, 6), legs):
        s = scores(A[k]); a = agree(C[1], C[k]); r = cc(s1[4], s[4])
        row.update({f"{n['side']}_quad": n["quad"], f"{n['side']}_role": n["role"], f"{n['side']}_held": s[1],
                    f"{n['side']}_agree": a, f"{n['side']}_data_r": r})
        ag.append(a); dr.append(r); nh.append(s[1])
    row.update(mean_agree=float(np.mean(ag)), min_agree=float(np.min(ag)), mean_data_r=float(np.mean(dr)), mean_neighbour_held=float(np.mean(nh)))
    rows.append(row)

if not rows:
    raise SystemExit("no finished centres yet")
keys = list(rows[0])
with open(HERE / "quad5_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
v = lambda k: np.array([r[k] for r in rows])  # noqa: E731
print(f"{len(rows)} of {len(nb)} centres finished")
print(f"centre leg, median: train {np.median(v('train')):+.3f}, held out {np.median(v('held')):+.3f}; trend removed: "
      f"train {np.median(v('train_detr')):+.3f}, held out {np.median(v('held_detr')):+.3f}")
print(f"same quads as single fits (Oct2026_dlod), held out median {np.median(v('single_fit_held')):+.3f}; "
      f"five-leg better in {np.mean(v('held') > v('single_fit_held')):.0%} of centres")
print(f"neighbour response agreement with the centre: median {np.median(v('mean_agree')):+.2f}; centres with any neighbour below 0: "
      f"{int(np.sum(v('min_agree') < 0))}; relation of agreement to data correlation across centres r = {cc(v('mean_agree'), v('mean_data_r')):+.2f}")
for lab, k in (("best held out", "held"), ("lowest neighbour agreement", "min_agree")):
    o = np.argsort(v(k))
    pick = o[-5:][::-1] if lab.startswith("best") else o[:5]
    print(f"{lab}: " + ", ".join(f"{rows[i]['centre']} {rows[i][k]:+.2f}" for i in pick))

fig, axs = plt.subplots(1, 2, figsize=(16, 5.2), facecolor=SURF)
for ax, k, title, cmap, lim in ((axs[0], "held", "centre leg, held-out correlation 2000-2005", "RdBu_r", 0.8),
                                (axs[1], "mean_agree", "mean response agreement of the 4 neighbours with the centre", "RdBu_r", 1.0)):
    ax.set_facecolor(SURF)
    sc = ax.scatter(v("lon"), v("lat"), c=v(k), cmap=cmap, vmin=-lim, vmax=lim, s=420, marker="s", edgecolor=SURF, lw=1)
    for r in rows:
        ax.text(r["lon"], r["lat"], f"{r[k]:+.2f}", ha="center", va="center", fontsize=6.5, color=INK)
    ax.set_xlim(-185, 185); ax.set_ylim(-75, 92); ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-60, 81, 20))
    ax.set_xlabel("longitude of quad centre", color=MUTED); ax.set_ylabel("latitude", color=MUTED)
    ax.set_title(title, loc="left", fontsize=11, color=INK); ax.tick_params(colors=MUTED); ax.grid(True, color=GRID, lw=0.5)
    for s in ax.spines.values():
        s.set_color(GRID)
    fig.colorbar(sc, ax=ax, shrink=0.8)
fig.suptitle(f"Five-leg fits, every Kaplan quad as centre with its four neighbours ({len(rows)} of {len(nb)} centres)",
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout()
fig.savefig(HERE / "quad5_maps.png", dpi=100, facecolor=SURF)
print("plot: quad5_maps.png; table: quad5_summary.csv")
