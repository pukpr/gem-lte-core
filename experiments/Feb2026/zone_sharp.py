#!/usr/bin/env python3
"""Sharpened BZONE (Brillouin-zone) fold of the LTE modulation, from lte_results.csv and lt.exe.p only.

plot.py's BZONE option (formerly ZONE) folds the latent forcing F (column 4) by the base winding k (the last ltep
entry), phase = (k F) mod 1, and plots the data against that phase. Only the base and its integer
harmonics are periodic in that phase. Everything else in the fit smears the fold:
    level + k0 F + trend t + accel (t - t0)^2      (the secular part and the slow ~60-yr arch in F)
    the low winding(s), a sin(2 pi k_low F + phi)  (not a harmonic of the base)
This script removes those before folding.

IR. The model is fitted to DR(t) = Data(t) + IR Data(t-12), and the saved model is
M(t) - IR M(t-12). As in plot.py, the data side uses DR and the model side the recursive inverse,
so the lag-12 term is already out of both; nothing further is needed for it.

Method (no lt.exe.windings.json needed; the regression is re-solved here):
  1. DR is regressed on  1, F, t, t^2, sin/cos(2 pi k_low F) for each low winding,
     and sin/cos(2 pi h k F) for the base (h = 1) and every harmonic h in the .p "harm" list.
  2. compensated data = DR - (fitted non-periodic terms); the periodic model is the fitted
     base + harmonics, a single-valued curve in the phase.
  3. Optional refinements, each kept only if the fold gets sharper:
       --tune     scan the fold factor around k (the phase is very sensitive to k over a wide F range)
       --seasonal add annual and semi-annual terms in t to the terms removed
Sharpness = share of the plotted data's variance explained by the phase alone (eta^2 from phase bins).

usage: zone_sharp.py [run_dir] [--exclude 2000 2005] [--tune] [--seasonal] [--bins 40] [--out file.png]
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("run_dir", nargs="?", default=".")
ap.add_argument("--exclude", nargs=2, type=float, default=None, metavar=("START", "END"),
                help="interval left out of the regression (the run's held-out gap)")
ap.add_argument("--tune", action="store_true", help="refine the fold factor around the base winding")
ap.add_argument("--seasonal", action="store_true", help="also remove annual and semi-annual terms")
ap.add_argument("--bins", type=int, default=40)
ap.add_argument("--out", default=None)
args = ap.parse_args()

d = args.run_dir
a = np.loadtxt(os.path.join(d, "lte_results.csv"), delimiter=",")
t, model, data, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
p = json.load(open(os.path.join(d, "lt.exe.p")))
ir = float(p.get("IR", p.get("ir", 0.0)))
ltep = [float(x) for x in p["ltep"]]
kbase, klow = ltep[-1], ltep[:-1]
harm = [1.0] + [float(h) for h in (p.get("harm") or [])]

# --- IR, exactly as plot.py does it -----------------------------------------------------------
DR = data.copy()
DR[12:] = data[12:] + ir * data[:-12]
Mo = model.copy()
if abs(ir) < 1.0:
    for i in range(12, len(model)):
        Mo[i] = model[i] + ir * Mo[i - 12]

fit = np.ones(len(t), bool)
if args.exclude:
    fit &= ~((t >= args.exclude[0]) & (t < args.exclude[1]))
fit[:12] = False                                   # DR has no lag-12 partner in the first year


def design(k, seasonal):
    """columns removed (non-periodic in the base phase) and columns kept (base + harmonics)"""
    tt = t - t[0]
    rem = [np.ones_like(t), F, tt, tt ** 2]
    for kl in klow:
        rem += [np.sin(2 * np.pi * kl * F), np.cos(2 * np.pi * kl * F)]
    if seasonal:
        for n in (1, 2):
            rem += [np.sin(2 * np.pi * n * t), np.cos(2 * np.pi * n * t)]
    keep = []
    for h in harm:
        keep += [np.sin(2 * np.pi * h * k * F), np.cos(2 * np.pi * h * k * F)]
    return np.column_stack(rem), np.column_stack(keep)


def eta2(y, phase, nb):
    """share of the variance of y explained by its phase-bin means"""
    b = np.minimum((phase * nb).astype(int), nb - 1)
    tot = np.sum((y - y.mean()) ** 2)
    within = sum(np.sum((y[b == j] - y[b == j].mean()) ** 2) for j in range(nb) if np.any(b == j))
    return 1.0 - within / tot


def solve(k, seasonal):
    R, K = design(k, seasonal)
    X = np.column_stack([R, K])
    c = np.linalg.lstsq(X[fit], DR[fit], rcond=None)[0]
    nonper = R @ c[:R.shape[1]]
    per = K @ c[R.shape[1]:]
    comp = DR - nonper
    ph = np.mod(k * F, 1.0)
    return dict(k=k, seasonal=seasonal, c=c, nR=R.shape[1], nonper=nonper, per=per, comp=comp, phase=ph,
                eta=eta2(comp[fit], ph[fit], args.bins),
                r_per=float(np.corrcoef(per[fit], comp[fit])[0, 1]),
                r_fit=float(np.corrcoef((nonper + per)[fit], DR[fit])[0, 1]))


ph0 = np.mod(kbase * F, 1.0)
eta_raw = eta2(DR[fit], ph0[fit], args.bins)
best = solve(kbase, False)
steps = [("as plot.py folds it (DR vs phase)", eta_raw, kbase),
         ("non-periodic terms removed", best["eta"], kbase)]
if args.seasonal:
    s = solve(kbase, True)
    steps.append(("+ annual and semi-annual removed", s["eta"], kbase))
    if s["eta"] > best["eta"]:
        best = s
if args.tune:
    # one cycle of phase slip across the range of F is dk = 1 / range(F); scan a few of those
    span = 3.0 / (F.max() - F.min())
    ks = kbase + np.linspace(-span, span, 241)
    cand = max((solve(k, best["seasonal"]) for k in ks), key=lambda s: s["eta"])
    steps.append((f"+ fold factor tuned ({kbase:.6f} -> {cand['k']:.6f})", cand["eta"], cand["k"]))
    if cand["eta"] > best["eta"]:
        best = cand

print(f"{os.path.abspath(d)}")
print(f"  IR {ir:+.4f} | base winding k {kbase:.6f} (zone = {1 / abs(kbase):.3f} in F) | low winding(s) {klow} | harmonics {harm[1:]}")
print(f"  manifold F range {F.min():+.3f} .. {F.max():+.3f} = {(F.max() - F.min()) * abs(kbase):.1f} zones")
print(f"  re-solved regression vs DR: CC {best['r_fit']:.4f}; saved model (IR removed) vs DR: CC {np.corrcoef(Mo[fit], DR[fit])[0, 1]:.4f};"
      f" re-solved vs saved model: CC {np.corrcoef((best['nonper'] + best['per'])[fit], Mo[fit])[0, 1]:.4f}")
print("  sharpness (share of plotted-data variance explained by the phase alone):")
for lab, e, k in steps:
    print(f"    {lab:55s} eta^2 {e:.3f}")
print(f"  periodic model (base + harmonics) vs compensated data: CC {best['r_per']:.3f}")
v = lambda x: float(np.var(x[fit]))
print(f"  variance shares of DR: non-periodic terms {v(best['nonper']) / v(DR):.2f}, periodic terms {v(best['per']) / v(DR):.2f},"
      f" residual {v(DR - best['nonper'] - best['per']) / v(DR):.2f}")

# --- figure: before / after, phase on the vertical axis as in the GUI's LTE Modulation panel -----
nb = args.bins
edges = np.linspace(0, 1, nb + 1)
ctr = 0.5 * (edges[1:] + edges[:-1])
def binned(y, ph):
    b = np.minimum((ph * nb).astype(int), nb - 1)
    return np.array([np.median(y[b == j]) if np.any(b == j) else np.nan for j in range(nb)])
pg = np.linspace(0, 1, 400)
c = best["c"][best["nR"]:]
curve = sum(c[2 * i] * np.sin(2 * np.pi * h * pg) + c[2 * i + 1] * np.cos(2 * np.pi * h * pg) for i, h in enumerate(harm))

fig, ax = plt.subplots(1, 3, figsize=(17, 6.5), facecolor="#e6f2ff")
ax[0].plot(DR, ph0, ls="None", marker="_", color="blue", alpha=0.6, label="Data (IR-matched)")
ax[0].plot(Mo, ph0, ls="None", marker="_", color="red", alpha=0.6, label="Model (IR removed)")
ax[0].set_title(f"Zone fold as plot.py draws it\n$\\eta^2$ = {eta_raw:.2f}")
ax[1].plot(best["comp"], best["phase"], ls="None", marker="_", color="blue", alpha=0.6, label="Data, non-periodic terms removed")
ax[1].plot(binned(best["comp"], best["phase"]), ctr, color="navy", lw=1.5, label="Data, median per phase bin")
ax[1].plot(curve, pg, color="red", lw=2.2, label="Model: base + harmonics")
ax[1].set_title(f"Compensated fold (k = {best['k']:.6f})\n$\\eta^2$ = {best['eta']:.2f}, model vs data CC {best['r_per']:.2f}")
for x in ax[:2]:
    x.set_ylim(0, 1); x.set_ylabel(f"Folded value (k F mod 1)"); x.set_xlabel("Level modulation"); x.legend(loc="upper left", fontsize=8)
# what was removed, against time, so the 60-yr part is visible
ax[2].plot(t, DR, color="blue", lw=0.6, alpha=0.6, label="Data (IR-matched)")
ax[2].plot(t, best["nonper"], color="black", lw=1.6, label="removed: secular + low winding" + (" + seasonal" if best["seasonal"] else ""))
ax[2].plot(t, best["per"], color="red", lw=0.7, alpha=0.8, label="kept: base + harmonics")
ax[2].set_title("What is removed before folding"); ax[2].set_xlabel("Year"); ax[2].legend(loc="upper left", fontsize=8)
fig.tight_layout()
out = args.out or os.path.join(d, "zone_sharp.png")
fig.savefig(out, dpi=110)
print(f"  plot: {out}")
