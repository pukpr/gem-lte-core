#!/usr/bin/env python3
"""ws_params.py — plots ws.py's LOCAL-KERNEL FIT's own per-winding
amplitude(t) and phase(t) trajectories, not just the resulting curve.

Motivation (direct user critique of ws.py --mode local): the fitted
CURVE's r-value is mostly self-referential flexibility (see ws.py's
--n-surrogates honesty check and the blind-holdout collapse it also
revealed) -- not, by itself, informative about the global model. What
IS a genuinely different question local mode can answer: does any
winding's own amplitude or phase drift systematically over the record?
A real, structural drift (not noise) is an actionable clue for
IMPROVING the global (production-matching) fit -- e.g. a winding whose
amplitude visibly strengthens in one multi-decade regime and fades in
another suggests either a genuine slow envelope/beat the global model's
single fixed amplitude can't capture, or (if two windings' amplitudes
move in lockstep) that they should be merged/reconsidered as one
near-degenerate pair rather than two independent terms (the same
mechanism already found this session for PDO/NAO/brestexcl's beat
pairs). A winding whose amplitude is small and noisy everywhere is a
candidate to drop from the global model instead.

Usage
-----
    ./ws_params.py --index nino4
    ./ws_params.py --index brestexcl --sigma 15 --root /some/path
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

import ws  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True)
    ap.add_argument("--root")
    ap.add_argument("--sigma", type=float, default=10.0,
                     help="local-fit kernel width, years")
    ap.add_argument("--nmax", type=int, default=3)
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

    print(f"[{args.index}] root={root}  sigma={args.sigma}yr  nmax={args.nmax}  Ms={Ms}")
    fit, params = ws.local_fit(F, t, x, args.sigma, Ms, nmax=args.nmax,
                                trend_on=prep["trend_on"], return_params=True)

    # EDGE MASK -- within one local-fit sigma of either end of the record,
    # the Gaussian window has support on only one side, so the regression
    # is under-constrained and prone to blowing up on noise (found
    # directly: nino4 showed every winding's amplitude spike upward in
    # its last ~sigma years with no plausible shared physical cause).
    # Matches winding_scalogram.py's own edge_mask convention (hatched,
    # not deleted -- the values are real model output, just low-
    # confidence) so the two tools read consistently.
    t_min, t_max = t.min(), t.max()
    edge_mask_pt = (t - t_min < args.sigma) | (t_max - t < args.sigma)
    interior = ~edge_mask_pt
    if interior.sum() < 10:
        print(f"[{args.index}] WARNING: record barely longer than 2*sigma -- "
              f"almost everything is edge-flagged; consider a smaller --sigma")

    n_M = len(Ms)
    fig, axes = plt.subplots(n_M, 1, figsize=(11, 1.8 * n_M), sharex=True,
                              squeeze=False)
    axes = axes[:, 0]
    for ax, M in zip(axes, Ms):
        amp = params[M]["amp"]
        ax2 = ax.twinx()
        ax.plot(t, amp, color="tab:blue", lw=1.1)
        ax.set_ylabel(f"M={M:g}\namplitude", fontsize=8, color="tab:blue")
        ax.tick_params(axis="y", labelcolor="tab:blue")
        # phase plotted unwrapped-ish via sin(phase) so it reads as a
        # smooth trace rather than a sawtooth at the +/-pi wrap
        ax2.plot(t, params[M]["phase"], color="0.7", lw=0.8, alpha=0.7)
        ax2.set_ylabel("phase (rad)", fontsize=7, color="0.5")
        ax2.tick_params(axis="y", labelcolor="0.5", labelsize=7)
        if t_min <= t_min + args.sigma:
            ax.axvspan(t_min, t_min + args.sigma, color="white", alpha=0.55,
                       hatch="//", linewidth=0, zorder=5)
        if t_max - args.sigma <= t_max:
            ax.axvspan(t_max - args.sigma, t_max, color="white", alpha=0.55,
                       hatch="//", linewidth=0, zorder=5)
        amp_int = amp[interior] if interior.sum() > 0 else amp
        amp_cv = float(np.std(amp_int) / np.mean(amp_int)) if np.mean(amp_int) > 0 else float("nan")
        ax.text(0.995, 0.88,
                f"amp mean={amp_int.mean():.3f}  CV={amp_cv:.2f}  "
                f"(interior only, hatched edges excluded)",
                transform=ax.transAxes, ha="right", fontsize=7, color="tab:blue")
    axes[-1].set_xlabel("year")
    fig.suptitle(f"{args.index}: local-fit per-winding amplitude(t)/phase(t), "
                 f"sigma={args.sigma:g}yr  (blue=amplitude, left axis; "
                 f"grey=phase, right axis; hatched = within one sigma of "
                 f"record edge, low-confidence -- see edge-bias note)",
                 fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out_dir = args.outdir if args.outdir is not None else (root / args.index)
    out_path = out_dir / f"{args.index}_ws_params.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"[saved] {out_path}")

    print(f"\n[{args.index}] per-winding amplitude coefficient of variation, "
          f"INTERIOR ONLY (excludes the hatched edge regions, "
          f"{t_min:.1f}-{t_min+args.sigma:.1f} and "
          f"{t_max-args.sigma:.1f}-{t_max:.1f}, where the local fit is "
          f"under-constrained and prone to edge-bias blowup -- CV = "
          f"std/mean; high CV means this winding's strength genuinely "
          f"drifts across the record; near-constant amplitude, low CV, "
          f"means the global model's single fixed amplitude is already a "
          f"good summary for that winding):")
    for M in Ms:
        amp = params[M]["amp"][interior] if interior.sum() > 0 else params[M]["amp"]
        cv = float(np.std(amp) / np.mean(amp)) if np.mean(amp) > 0 else float("nan")
        edge_amp = params[M]["amp"][edge_mask_pt]
        edge_flag = ""
        if edge_mask_pt.sum() > 0 and amp.mean() > 0 and edge_amp.max() > 3 * amp.mean():
            edge_flag = f"  <-- edge amplitude peaks at {edge_amp.max():.3f}, " \
                        f">3x the interior mean: likely edge-bias artifact, not real"
        print(f"    M={M:8.4f}   interior mean amp={amp.mean():.4f}   "
              f"CV={cv:.3f}{edge_flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
