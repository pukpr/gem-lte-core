"""Vestige cancellation and cubic cross-terms across the Nino regions.

For each leg of a multi-leg run (one shared manifold F, per-leg windings), fit the slow-winding
response y = k0*F + sum A sin(2 pi k F + phi) (slow windings only) as a cubic in u = standardised F,
split the cubic term into its in-phase part (projection onto u: the vestige gain) and the residual
cross-terms, and compare peak periods of model and data in 2-8 yr.

  python3 regions_crossterms.py [--run DIR] [--labels NINO4 NINO34 NINO3 NINO12]
"""
import argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import peak_estimate as pe

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="../../experiments/Oct2026_generator/legs/nino_regions_seed_nino4")
ap.add_argument("--labels", nargs="+", default=["NINO4", "NINO34", "NINO3", "NINO12"])
ap.add_argument("--bench", default="../../experiments/Oct2026_generator/legs/nino4_slow_then_harm/benchmark_enso_2026-11-01")
ap.add_argument("--out", default="regions_crossterms.png")
a = ap.parse_args()


def spec(x):
    f, Pk = pe.multitaper(pe.prep(x)); return f, Pk.mean(0)


def pk(x, n=3):
    f, P = spec(x); return [1 / v for v in pe.band_peaks(f, P, n)]


def band_corr(x, y, lo=2.0, hi=8.0):
    f, P = spec(x); _, Q = spec(y); m = (f >= 1 / hi) & (f <= 1 / lo)
    return np.corrcoef(P[m], Q[m])[0, 1]


def band_power(x, lo, hi):
    f, P = spec(x); m = (f >= 1 / hi) & (f <= 1 / lo); return P[m].sum()


sets = [("NINO4 benchmark (1 leg)", f"{a.bench}/lte_results.csv", f"{a.bench}/lt.exe.windings.json")]
sets += [(f"{lab} (4-leg)", f"{a.run}/lte_results.{i}.csv", f"{a.run}/lt.exe.windings.{i}.json") for i, lab in enumerate(a.labels, 1)]

fmt = lambda v: ", ".join(f"{p:.2f}" for p in v)
print(f"{'series':26s} {'data peaks':20s} {'model peaks':20s} {'cross-term peaks':20s} {'spec r':>6s} {'c1':>6s} {'c3':>6s} {'net vestige':>11s} {'3.5-4.1 drop':>12s} {'cubic r':>7s}")
rows = []
for name, res, wf in sets:
    d = np.genfromtxt(res, delimiter=","); t, model, data, F = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    ok = data != 0  # leading zero padding
    first = np.argmax(ok); t, model, data, F = t[first:], model[first:], data[first:], F[first:]
    w = json.load(open(wf)); kmin = min(k for k, _, _ in w["k_amp_phase"])
    slow = [x for x in w["k_amp_phase"] if x[0] < 3 * kmin]
    y = w["k0"] * F + sum(A * np.sin(2 * np.pi * k * F + ph) for k, A, ph in slow)
    u = (F - F.mean()) / F.std()
    c = np.polyfit(u, y, 3)[::-1]
    lin, cub = c[1] * u, c[3] * u ** 3
    beta = (cub @ u) / (u @ u); resid = cub - beta * u
    r3 = np.corrcoef(np.polyval(c[::-1], u), y)[0, 1]
    drop = band_power(lin, 3.5, 4.1) / band_power(lin + cub, 3.5, 4.1)
    sr = band_corr(model, data)
    rows.append((name, t, model, data))
    print(f"{name:26s} {fmt(pk(data)):20s} {fmt(pk(model)):20s} {fmt(pk(resid)):20s} {sr:6.2f} {c[1]:6.2f} {c[3]:6.2f} {c[1] + beta:+11.2f} {drop:11.1f}x {r3:7.3f}")

fig, axes = plt.subplots(len(rows), 1, figsize=(10, 2.3 * len(rows)), sharex=True)
for ax, (name, t, model, data) in zip(axes, rows):
    for x, col, lab in [(data, "#555555", "data"), (model, "#1b4f72", "model")]:
        f, P = spec(x); m = (f > 1 / 12) & (f < 1)
        ax.plot(1 / f[m], P[m] / P[m].max(), color=col, lw=1.5, label=lab)
    for v, ls in [(3.796, "-"), (4.762, "--")]:
        ax.axvline(v, color="#c0392b", lw=1, ls=ls)
    ax.set_title(name, fontsize=10, loc="left"); ax.grid(alpha=.3); ax.set_xlim(1.5, 10)
axes[0].legend(fontsize=8, loc="upper right")
axes[0].plot([], [], color="#c0392b", label="Mf 3.796 / Mf' 4.762 yr aliases")
axes[-1].set_xlabel("period (yr)")
fig.tight_layout(); fig.savefig(a.out, dpi=120); print(f"\nplot: {a.out}")
