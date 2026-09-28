#!/usr/bin/env python3
"""resolved_winding_scalogram.py — winding_scalogram.py's own
transform (winding_transform/noise_floor/plot_panel, unmodified, no new
machinery), pointed at a RESOLVED SHALLOW-WATER MODEL's own gauge output
instead of the production Ada regression's Data/Model columns.

Direct response to a specific methodological question: across
baltic/nao/brestexcl/pna, the resolved-PDE gauge output's own visual
character looks like it's mostly just PASSING THE FORCING'S OWN
already-"squared-off" (Impulse_Delta comb + near-integrator IIR) shape
straight through -- i.e. a largely non-autonomous, undamped transmission
-- except for AMO, where genuine PDE-driven shaping (the moving-gauge
frequency-doubling) is visible. The 12mo delayed-difference is a
genuinely autonomous (self-referential) operation, but since it's
applied to an already-squared-off input, its output stays squared-off
too, rather than becoming smooth. And PNA's own real high-order
harmonics (M=4.73 especially) don't show up anywhere in this project's
own ridge-scan checks of the resolved models, raising the question of
whether the PDE is even CAPABLE of exciting/carrying that content, or is
acting as an inadvertent low-pass filter.

The direct, no-new-degrees-of-freedom way to check any of this: run the
exact same scalogram transform used throughout this project on real
data against the resolved model's OWN gauge output against its OWN
driving Forcing, side by side with the real data panel -- if the model
is just passing the forcing through, its scalogram should look like a
smeared/scaled version of the SAME comb-and-IIR structure at every M
(no real M-selective structure of its own); if it's genuinely creating
winding structure via the shallow-water dynamics, the scalogram should
show real ridges, including possibly at the higher M values the model's
own ridge-scan checks have been missing.

Usage
-----
    ./resolved_winding_scalogram.py --index pna --variant fixed
    ./resolved_winding_scalogram.py --index pna --variant moving --m-max 5
    ./resolved_winding_scalogram.py --index nao --variant diff
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from winding_scalogram import (standardize, winding_transform, noise_floor,  # noqa: E402
                               compute_log_power, robust_scale, plot_panel)
from signal_operators import delayed_difference                           # noqa: E402

ROOT = Path(__file__).resolve().parent

# each already-built resolved model's own cache file + established real
# winding_rank ridges, gathered in one place for reuse
CACHES = {
    "pna": ("pna/pna_resolved_check_cache.npz", [1.460, 4.730, 2.970]),
    "brestexcl": ("brestexcl/brestexcl_resolved_check_cache.npz",
                 [1.870, 0.955, 2.490, 0.410, 0.060]),
    "nao": ("nao/nao_baltic_check_cache.npz", [0.8307]),
    "baltic": ("baltic/baltic_atmospheric_check_cache.npz", [1.2454]),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, choices=list(CACHES))
    ap.add_argument("--variant", default="fixed",
                     choices=["fixed", "moving", "diff"])
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.02)
    ap.add_argument("--sigma", type=float, default=15.0)
    ap.add_argument("--t0-step", type=float, default=5.0)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    cache_path, fitted_m = CACHES[args.index]
    d = np.load(ROOT / cache_path)
    gauge_dates, gauge_forcing, fixed, moving = (
        d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"])

    if args.variant == "fixed":
        series, t, forcing = fixed, gauge_dates, gauge_forcing
    elif args.variant == "moving":
        series, t, forcing = moving, gauge_dates, gauge_forcing
    else:
        series = delayed_difference(fixed, 12)
        t = gauge_dates[12:]
        forcing = gauge_forcing[12:]

    dat = np.loadtxt(ROOT / args.index / f"{args.index}.dat")
    obs = np.interp(t, dat[:, 0], dat[:, 1])

    m_grid = np.arange(0.0, args.m_max + args.dm / 2, args.dm)
    t0_grid = np.arange(t[0] + args.sigma / 2, t[-1] - args.sigma / 2 + 1e-9,
                        args.t0_step)
    if len(t0_grid) < 2:
        raise SystemExit(f"record too short for sigma={args.sigma}yr windows")

    obs_s = standardize(obs)
    series_s = standardize(series)

    print(f"-- {args.index} ({args.variant}): resolved-model winding "
          f"scalogram vs real data, both against the model's own driving "
          f"Forcing --")
    G_obs, edge_mask, _ = winding_transform(t, obs_s, forcing, m_grid,
                                            t0_grid, args.sigma)
    G_model, _, _ = winding_transform(t, series_s, forcing, m_grid, t0_grid,
                                      args.sigma)
    floor = noise_floor(t, forcing, m_grid, t0_grid, args.sigma)

    log_power_obs = compute_log_power(G_obs, floor)
    log_power_model = compute_log_power(G_model, floor)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True,
                             layout="constrained")
    vmin, vmax = robust_scale(np.concatenate([log_power_obs.ravel(),
                                              log_power_model.ravel()]))
    pcm1 = plot_panel(axes[0], t0_grid, m_grid, log_power_obs, edge_mask,
                      fitted_m, f"{args.index}: REAL DATA winding against "
                      f"the resolved model's own Forcing", vmin, vmax)
    pcm2 = plot_panel(axes[1], t0_grid, m_grid, log_power_model, edge_mask,
                      fitted_m, f"{args.index}: RESOLVED-PDE MODEL "
                      f"({args.variant}) winding against its own Forcing",
                      vmin, vmax)
    axes[0].set_xlabel("year (window center)")
    axes[1].set_xlabel("year (window center)")
    axes[1].set_ylabel("")
    cb = fig.colorbar(pcm2, ax=axes, pad=0.01, shrink=0.9)
    cb.set_label("log2 power / matched-noise floor")
    fig.suptitle(f"Resolved-model winding scalogram -- {args.index} "
                f"({args.variant}) (window sigma={args.sigma:g}yr; dashed "
                f"red = this index's own established real ridges; hatched "
                f"= within one sigma of record edge)", fontsize=9)

    out_dir = args.outdir if args.outdir is not None else ROOT / args.index
    out_path = out_dir / f"resolved_winding_scalogram_{args.variant}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
