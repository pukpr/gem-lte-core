#!/usr/bin/env python3
"""Two things WINDING_NARRATIVE.md was missing, called out directly: (1) no
time-series graphs, so a claimed cycle (e.g. AMO's dominant M=-0.0134 "60
year" mode) can't be checked by eye against the real data; (2) no
demonstration that the search actually covers/resolves the claimed optimum,
so there's no way to tell a genuine, well-defined fit from wherever the
optimizer happened to land.

This script builds both, per index:
  1. Isolates the single largest-amplitude fitted mode (base or harmonic),
     reconstructs it alone as A*sin(2*pi*M*Forcing(t)+phase) using the
     REAL fitted M/amp/phase, and overlays it directly against the
     standardized Data and Model time series -- the reader can see whether
     the claimed cycle is actually there, not just read an M value.
  2. Sweeps candidate M across a range spanning the claimed mode, doing a
     direct OLS regression (offset + Forcing + sin(2*pi*M*Forcing) +
     cos(2*pi*M*Forcing)) against Data at each candidate M, and plots the
     resulting fit-quality (correlation) landscape. A genuine, well-found
     optimum should show a clear, well-defined peak at the reported M --
     not a flat landscape where many candidates fit equally well (which
     would mean the reported M is just wherever the optimizer happened to
     land, not something the search actually resolved).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns, standardize
from param_survey import compute_winding

ROOT = Path(__file__).resolve().parent


def fit_at_M(forcing, data, M):
    """Direct OLS regression at a fixed candidate M: offset + k0*Forcing +
    A*sin(2*pi*M*Forcing) + B*cos(2*pi*M*Forcing). Returns correlation
    between the fit and data."""
    design = np.column_stack([
        np.ones_like(forcing), forcing,
        np.sin(2 * np.pi * M * forcing), np.cos(2 * np.pi * M * forcing),
    ])
    coef, _, _, _ = np.linalg.lstsq(design, data, rcond=None)
    fitted = design @ coef
    return np.corrcoef(fitted, data)[0, 1]


def sweep_M(forcing, data, m_range, n_points=400):
    m_grid = np.linspace(m_range[0], m_range[1], n_points)
    cc = np.array([fit_at_M(forcing, data, m) for m in m_grid])
    return m_grid, cc


def main(idx: str, m_sweep_halfwidth: float = 0.5, n_sweep: int = 400):
    year, dt, model, obs, forcing = load_columns(idx)
    w = compute_winding(idx)
    if w is None:
        print(f"[{idx}] no fit available, skipping")
        return
    m = np.asarray(w["m"]); amp = np.asarray(w["amp"]); ph = np.asarray(w["phase"])
    i_dom = np.argmax(np.abs(amp))
    M0, A0, P0 = m[i_dom], amp[i_dom], ph[i_dom]
    print(f"[{idx}] dominant mode: M={M0:.4f}  amp={A0:.4f}  phase={P0:.4f}")

    isolated = A0 * np.sin(2 * np.pi * M0 * forcing + P0)

    data_s = standardize(obs)
    model_s = standardize(model)
    isolated_s = standardize(isolated)

    # search-coverage check: does a direct M-sweep resolve a clear peak at
    # M0, or is the landscape flat/ambiguous there?
    lo, hi = M0 - m_sweep_halfwidth, M0 + m_sweep_halfwidth
    m_grid, cc_curve = sweep_M(forcing, obs, (lo, hi), n_sweep)
    cc_at_M0 = fit_at_M(forcing, obs, M0)
    best_idx = np.argmax(cc_curve)
    print(f"[{idx}] CC at reported M0: {cc_at_M0:.4f};  "
          f"best CC in sweep: {cc_curve[best_idx]:.4f} at M={m_grid[best_idx]:.4f}  "
          f"(offset from M0: {m_grid[best_idx]-M0:+.4f})")

    fig, axes = plt.subplots(2, 1, figsize=(12, 9))

    axes[0].plot(year, data_s, color="black", linewidth=0.8, label="Data (standardized)")
    axes[0].plot(year, model_s, color="tab:blue", linewidth=0.7, alpha=0.6,
                 label="Model (standardized)")
    axes[0].plot(year, isolated_s, color="tab:red", linewidth=1.6,
                 label=f"isolated dominant mode alone: M={M0:.4f} "
                       f"(amp={A0:.3f}, {'strongest' if i_dom==np.argmax(np.abs(amp)) else ''} term)")
    axes[0].set_xlabel("year")
    axes[0].set_ylabel("standardized value")
    axes[0].set_title(f"{idx}: real Data/Model vs. the isolated dominant "
                       f"winding mode alone -- does the claimed cycle "
                       f"actually show up?", fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].plot(m_grid, cc_curve, color="tab:green", linewidth=1.2)
    axes[1].axvline(M0, color="red", linestyle="--", linewidth=1.0,
                     label=f"reported M0={M0:.4f} (CC={cc_at_M0:.3f})")
    axes[1].axvline(m_grid[best_idx], color="gray", linestyle=":",
                     linewidth=1.0,
                     label=f"sweep best M={m_grid[best_idx]:.4f} "
                           f"(CC={cc_curve[best_idx]:.3f})")
    axes[1].set_xlabel("candidate winding number M")
    axes[1].set_ylabel("correlation (direct OLS fit vs. Data)")
    axes[1].set_title(f"search-coverage check: does a direct sweep near M0 "
                       f"find a well-defined peak there, or is it "
                       f"arbitrary/flat?", fontweight="bold")
    axes[1].legend(loc="best", fontsize=8)

    fig.tight_layout()
    out = ROOT / idx / "mode_evidence.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    idx = sys.argv[1] if len(sys.argv) > 1 else "amo"
    halfwidth = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    main(idx, halfwidth)
