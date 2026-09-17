#!/usr/bin/env python3
"""Cross-wavelet scalograms: Forcing (col 4) vs Data (col 3), and
Forcing (col 4) vs Model (col 2), from each index's lte_results.csv.

lte_results.csv columns (1-indexed, matching lte_gui.py's own usage):
  1 = Year (date)   2 = Model   3 = Data   4 = Forcing

Year is the natural, already-uniform time axis (monthly steps) — no
resampling is needed for these datasets, but the script re-grids onto a
uniform axis defensively in case a run has any gaps.

Method: a from-scratch Morlet continuous wavelet transform (Torrence &
Compo 1998), no extra dependencies beyond numpy/matplotlib — consistent
with this project's "modest compute" recipe (see RECIPE_FOR_RESEARCHERS.md).
For each pair (Forcing, X), the cross-wavelet transform Wxy = W_forcing *
conj(W_x) gives a scalogram of shared power at each (year, period), plus a
phase difference at each point (visualized as arrows): pointing right means
Forcing and X are in phase at that time/period; left means anti-phase;
up/down means a quarter-cycle lead/lag. The hatched region is the cone of
influence (COI) — periods there are long enough, relative to distance from
the record's edges, that edge effects make the result unreliable.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DEFAULT_INDICES = ["nino4", "pdo", "amo", "tna", "nao", "baltic", "iode", "emi", "nino34", "88", "darwin1880", "111", "qbo30"]
# DEFAULT_INDICES = ["pdo_iaaft_detuned_surrogate", "qbo30", "qbo50", "npi", "pna", "noi", "tpi", "kap10-10-20-30"]
# DEFAULT_INDICES = ["155", "brestexcl"]


def morlet_cwt(x, dt, dj=0.125, s0=None, w0=6.0):
    """Continuous Morlet wavelet transform via FFT (Torrence & Compo 1998).

    Returns (W, periods, coi_period) where W has shape (n_scales, len(x)),
    periods is the Fourier period (same time units as dt) for each scale,
    and coi_period(t) is the cone-of-influence period threshold at each
    time index.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    x = x - x.mean()

    n_pad = int(2 ** np.ceil(np.log2(n)))
    xpad = np.zeros(n_pad)
    xpad[:n] = x
    fx = np.fft.fft(xpad)
    k = 2.0 * np.pi * np.fft.fftfreq(n_pad, d=dt)

    if s0 is None:
        s0 = 2.0 * dt
    j = int(np.floor(np.log2(n * dt / s0) / dj))
    scales = s0 * 2.0 ** (np.arange(j + 1) * dj)

    W = np.empty((len(scales), n), dtype=complex)
    for i, s in enumerate(scales):
        norm = (np.pi ** -0.25) * np.sqrt(2.0 * np.pi * s / dt)
        daughter = norm * np.exp(-0.5 * (s * k - w0) ** 2)
        daughter[k <= 0] = 0.0
        W[i] = np.fft.ifft(fx * daughter)[:n]

    fourier_factor = 4.0 * np.pi / (w0 + np.sqrt(2.0 + w0 ** 2))
    periods = fourier_factor * scales

    t_idx = np.arange(n)
    dist_from_edge = np.minimum(t_idx, n - 1 - t_idx) * dt
    coi_period = fourier_factor * np.sqrt(2.0) * dist_from_edge

    return W, periods, coi_period


def load_columns(idx: str):
    csv_path = ROOT / idx / "lte_results.csv"
    data = np.loadtxt(csv_path, delimiter=",", usecols=(0, 1, 2, 3))
    year, model, obs, forcing = data.T

    dt = float(np.median(np.diff(year)))
    grid = np.arange(year[0], year[-1] + dt / 2, dt)
    model_g = np.interp(grid, year, model)
    obs_g = np.interp(grid, year, obs)
    forcing_g = np.interp(grid, year, forcing)
    return grid, dt, model_g, obs_g, forcing_g


def standardize(x):
    x = np.asarray(x, dtype=float)
    s = x.std()
    return (x - x.mean()) / s if s > 0 else x - x.mean()


def plot_pair(ax, year, dt, w_forcing, periods, coi_period, w_other,
              title, w0):
    cross = w_forcing * np.conj(w_other)
    power = np.abs(cross) ** 2
    power = np.maximum(power, 1e-12)
    log_power = np.log2(power)

    Y, X = np.meshgrid(periods, year)
    pcm = ax.pcolormesh(X.T, Y.T, log_power, shading="auto", cmap="viridis")

    ax.fill_between(year, coi_period, periods.max(), color="white",
                     alpha=0.35, hatch="//", edgecolor="black", linewidth=0.0)
    ax.plot(year, np.clip(coi_period, periods.min(), periods.max()),
            color="black", linewidth=1.0, linestyle="--")

    # Sparse phase-difference arrows, skipped inside the COI to reduce clutter.
    phase = np.angle(cross)
    nt, ns = len(year), len(periods)
    t_step = max(nt // 40, 1)
    s_step = max(ns // 14, 1)
    for si in range(0, ns, s_step):
        for ti in range(0, nt, t_step):
            if periods[si] > coi_period[ti]:
                continue
            ang = phase[si, ti]
            ax.annotate("", xy=(year[ti] + 0.35 * np.cos(ang) * dt * t_step,
                                 periods[si] * (2 ** (0.12 * np.sin(ang)))),
                        xytext=(year[ti], periods[si]),
                        arrowprops=dict(arrowstyle="->", color="white",
                                        linewidth=0.6, alpha=0.8))

    ax.set_yscale("log", base=2)
    ax.set_ylim(periods.min(), periods.max())
    ax.set_ylabel("period (years)")
    ax.set_title(title, fontweight="bold")
    yticks = [2 ** p for p in range(int(np.floor(np.log2(periods.min()))),
                                     int(np.ceil(np.log2(periods.max()))) + 1)]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{y:g}" for y in yticks])
    return pcm


def make_scalogram(idx: str, w0: float, dj: float, outdir: Path | None):
    csv_path = ROOT / idx / "lte_results.csv"
    if not csv_path.exists():
        print(f"  (skipping {idx}: no lte_results.csv)", file=sys.stderr)
        return

    year, dt, model, obs, forcing = load_columns(idx)

    w_forcing, periods, coi_period = morlet_cwt(standardize(forcing), dt,
                                                 dj=dj, w0=w0)
    w_obs, _, _ = morlet_cwt(standardize(obs), dt, dj=dj, w0=w0)
    w_model, _, _ = morlet_cwt(standardize(model), dt, dj=dj, w0=w0)

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
    pcm1 = plot_pair(axes[0], year, dt, w_forcing, periods, coi_period,
                      w_obs, f"{idx}: Forcing (col 4) x Data (col 3)", w0)
    pcm2 = plot_pair(axes[1], year, dt, w_forcing, periods, coi_period,
                      w_model, f"{idx}: Forcing (col 4) x Model (col 2)", w0)
    axes[1].set_xlabel("year")

    for pcm, ax in ((pcm1, axes[0]), (pcm2, axes[1])):
        cb = fig.colorbar(pcm, ax=ax, pad=0.01)
        cb.set_label("log2 cross-wavelet power")

    fig.suptitle(f"Cross-wavelet scalograms — {idx}  (Morlet w0={w0:g}, "
                 f"dt={dt:.4f}yr; hatched = cone of influence; arrows = "
                 f"phase, right=in-phase, left=anti-phase)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out_dir = outdir if outdir is not None else (ROOT / idx)
    out_path = out_dir / "wavelet_scalogram.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("indices", nargs="*", default=DEFAULT_INDICES,
                     help="index subdirectory names (default: all seven "
                          "survey indices)")
    ap.add_argument("--w0", type=float, default=6.0,
                     help="Morlet wavelet parameter (default 6.0)")
    ap.add_argument("--dj", type=float, default=0.125,
                     help="scale resolution, octaves per step (default 0.125)")
    ap.add_argument("--outdir", type=Path, default=None,
                     help="write all PNGs here instead of each index's own "
                          "subdirectory")
    args = ap.parse_args()

    for idx in args.indices:
        print(f"-- {idx} --")
        make_scalogram(idx, args.w0, args.dj, args.outdir)


if __name__ == "__main__":
    main()
