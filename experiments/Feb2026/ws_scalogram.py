#!/usr/bin/env python3
"""ws_scalogram.py — winding-power scalogram for ws.py's LOCAL-KERNEL FIT.

winding_scalogram.py's scalogram (G(t0,M) windowed correlation against the
Forcing basis, normalized by a matched white-noise floor) is built around
the PRODUCTION Ada model's saved output (<idx>/lte_results.csv, columns
Data/Model/Forcing). This script computes the exact same transform via
winding_scalogram.py's own library functions (winding_transform,
noise_floor, compute_log_power, robust_scale, plot_panel -- imported, not
copied), but feeds it ws.py's own (obs, local-kernel-fit) pair instead of
(Data, Model) -- so you can see whether ws.py's local fit is riding a
genuine, time-persistent winding ridge or just locally-flexible texture,
the same visual diagnostic winding_scalogram.py gives the production
model, applied to the OTHER model this session built.

Usage
-----
    ./ws_scalogram.py --index nino4
    ./ws_scalogram.py --index nino4 --sigma 10 --m-max 5 --scalogram-sigma 15
    ./ws_scalogram.py --index brestexcl --root /some/path --cmap hot

Panels: OBS winding-power (left) vs ws.py LOCAL FIT winding-power (right),
sharing one color scale, with ws.py's own fitted winding numbers (the
same Ms the local fit's design matrix actually used) drawn as reference
lines -- directly analogous to winding_scalogram.py's Data/Model panels.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import ws  # noqa: E402 -- load_index, find_index_root, local_fit
from winding_scalogram import (winding_transform, noise_floor,  # noqa: E402
                                compute_log_power, robust_scale, plot_panel)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True)
    ap.add_argument("--root", help="see ws.py --root (default: auto-discover)")
    ap.add_argument("--sigma", type=float, default=10.0,
                     help="ws.py LOCAL FIT kernel width, years (default 10)")
    ap.add_argument("--nmax", type=int, default=3,
                     help="ws.py LOCAL FIT max harmonic order per winding")
    ap.add_argument("--m-max", type=float, default=5.0,
                     help="scalogram: max winding number shown (default 5.0)")
    ap.add_argument("--dm", type=float, default=0.01,
                     help="scalogram: winding-number grid resolution")
    ap.add_argument("--scalogram-sigma", type=float, default=15.0,
                     help="scalogram: time-window half-width, years (default 15) "
                          "-- distinct from --sigma, which is ws.py's OWN local-fit "
                          "kernel width; the two need not match")
    ap.add_argument("--t0-step", type=float, default=5.0,
                     help="scalogram: spacing between window centers, years")
    ap.add_argument("--n-reps", type=int, default=32,
                     help="noise-floor white-noise surrogate reps")
    ap.add_argument("--cmap", choices=["viridis", "hot"], default="viridis")
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    root = ws.find_index_root(args.index, args.root)
    prep = ws.load_index(root, args.index)

    nz = prep["data"] != 0.0
    t = prep["dates"][nz]
    x = prep["data"][nz]
    x = (x - x.mean()) / x.std()
    F = prep["forcing"][nz]
    Ms = sorted({round(abs(mm), 4) for mm in prep["m"] if 0.002 < abs(mm) < 10})

    print(f"[{args.index}] root={root}  local-fit sigma={args.sigma}yr  "
          f"nmax={args.nmax}  Ms={Ms}")
    y = ws.local_fit(F, t, x, args.sigma, Ms, nmax=args.nmax,
                      trend_on=prep["trend_on"])
    r_raw = float(np.corrcoef(x, y)[0, 1])
    print(f"[{args.index}] local fit r_raw={r_raw:+.4f} (see ws.py --mode local "
          f"for the surrogate-honesty number -- not recomputed here)")

    m_grid = np.arange(0.0, args.m_max + args.dm / 2, args.dm)
    t0_grid = np.arange(t[0] + args.scalogram_sigma / 2,
                         t[-1] - args.scalogram_sigma / 2 + 1e-9, args.t0_step)
    if len(t0_grid) < 2:
        print(f"[{args.index}] record too short for --scalogram-sigma "
              f"{args.scalogram_sigma} -- reduce it", file=sys.stderr)
        return 1

    G_obs, edge_mask, _ = winding_transform(t, x, F, m_grid, t0_grid,
                                             args.scalogram_sigma)
    G_fit, _, _ = winding_transform(t, y, F, m_grid, t0_grid,
                                     args.scalogram_sigma)
    floor = noise_floor(t, F, m_grid, t0_grid, args.scalogram_sigma,
                         n_reps=args.n_reps)

    log_power_obs = compute_log_power(G_obs, floor)
    log_power_fit = compute_log_power(G_fit, floor)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True,
                              layout="constrained")
    vmin, vmax = robust_scale(np.concatenate([log_power_obs.ravel(),
                                               log_power_fit.ravel()]))
    pcm1 = plot_panel(axes[0], t0_grid, m_grid, log_power_obs, edge_mask,
                       Ms, f"{args.index}: OBS winding against Forcing",
                       vmin, vmax, cmap=args.cmap)
    pcm2 = plot_panel(axes[1], t0_grid, m_grid, log_power_fit, edge_mask,
                       Ms, f"{args.index}: ws.py LOCAL FIT (sigma={args.sigma:g}yr, "
                           f"n<={args.nmax}, r_raw={r_raw:+.3f}) winding against Forcing",
                       vmin, vmax, cmap=args.cmap)
    axes[0].set_xlabel("year (window center)")
    axes[1].set_xlabel("year (window center)")
    axes[1].set_ylabel("")
    cb = fig.colorbar(pcm2, ax=axes, pad=0.01, shrink=0.9)
    cb.set_label("log2 power / matched-noise floor")
    fig.suptitle(
        f"ws.py local-fit winding scalogram — {args.index}  (scalogram window "
        f"sigma={args.scalogram_sigma:g}yr; red lines = ws.py's own local-fit "
        f"windings; hatched = within one scalogram-sigma of record edge; "
        f"color = power relative to a white-noise floor through the same "
        f"Forcing basis)", fontsize=9)

    out_dir = args.outdir if args.outdir is not None else (root / args.index)
    out_path = out_dir / f"{args.index}_ws_scalogram.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[saved] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
