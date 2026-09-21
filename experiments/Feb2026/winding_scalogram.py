#!/usr/bin/env python3
"""Winding scalograms: how strongly Data (col 3) and Model (col 2) contain a
component "winding" at rate M around the Forcing signal (col 4), localized
in time (year).

This is deliberately NOT a conventional time-frequency wavelet scalogram.
In this model (see `LTE` in src/gem-lte-primitives.adb), a mode with
wave-number M contributes a term

    Amplitude * Sin(2*pi * M * Forcing(t) + Phase)

The winding argument is `2*pi * M * Forcing(t)` — M winds around the
*forcing signal's own value*, not around calendar time. A standard
wavelet's scale axis is conjugate to time (cycles per year); M is
conjugate to Forcing itself. These are different clocks, so a standard
scalogram's period axis can't just be relabeled into M — a different
transform is needed.

What this script computes instead, for each candidate M and each time
window centered at t0:

    G_X(t0, M) = < window(t - t0) * X(t) * exp(-i*2*pi*M*Forcing(t)) >
                 -------------------------------------------------------
                          < window(t - t0) >

for X in {Data, Model} — a windowed (Gabor-style) correlation of X against
the model's own native winding basis, localized in time. |G_X(t0,M)| is
the winding-power scalogram: how strongly X locally resembles a winding at
rate M around Forcing, as a function of when in the record. Dashed
horizontal lines mark the M values this index's own regression search
actually found (from `compute_winding` in param_survey.py), so you can see
directly whether that power is persistent across the whole record or
concentrated in a particular era.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns, standardize, DEFAULT_INDICES
from param_survey import compute_winding

ROOT = Path(__file__).resolve().parent


def load_ir(idx: str) -> float:
    p = json.loads((ROOT / idx / "lt.exe.p").read_text())
    return float(p.get("IR", p.get("ir", 0.0)))


def uncompensate_yearly_feedback(values, ir):
    """Invert the Ada model's post-LTE lag-12 IR adjustment
    (gem-lte-primitives-solution.adb:1062-1066): Model_final(I) =
    Model_orig(I) - ir*Model_orig(I-12), both taps reading the ORIGINAL
    (pre-adjustment) Model. Recovering Model_orig needs the recursive
    inverse built up in increasing-index order using the ALREADY-RECOVERED
    lag-12 value, not the raw one -- see plot.py's copy of this function
    for the round-trip verification. Apply only to Model; Data is never
    IR-adjusted on the Ada side (Save() writes the raw Data_Records, not
    the separate IR-adjusted DR used only internally as a regression
    target), so decompensating it would introduce a spurious transform.

    Unrolled, this recursion is a geometric series in ir and only stays
    bounded for |ir| < 1 (regardless of sign) -- for |ir| >= 1 (e.g.
    tpi's ir=-1.21) it diverges across the record. Guarded: skip and
    return the series unmodified rather than returning garbage."""
    values = np.asarray(values, dtype=float)
    uncompensated = values.copy()
    if abs(ir) >= 1.0:
        print(f"WARNING: |IR|={abs(ir):.4f} >= 1 -- decompensation is "
              f"unstable at this magnitude. Skipping for this series.",
              file=sys.stderr)
        return uncompensated
    if ir != 0.0:
        for i in range(12, len(values)):
            uncompensated[i] = values[i] + ir * uncompensated[i - 12]
    return uncompensated


def detect_gaps(idx: str, threshold: float = 0.5):
    """Real gaps in the underlying lte_results.csv time column (not the
    interpolated grid load_columns produces) -- e.g. brestexcl has a
    1944.33-1954.33 gap where no observations exist at all. Returns a list
    of (gap_start, gap_end) in years, for any consecutive-row jump larger
    than `threshold` years."""
    year = np.loadtxt(ROOT / idx / "lte_results.csv", delimiter=",",
                       usecols=(0,))
    diffs = np.diff(year)
    gaps = []
    for i in np.where(diffs > threshold)[0]:
        gaps.append((year[i], year[i + 1]))
    return gaps


def valid_mask_from_gaps(t, gaps):
    """True at grid points t that fall outside every detected gap -- i.e.
    NOT the fabricated linear-interpolation filler load_columns produces
    across a real gap (see detect_gaps)."""
    valid = np.ones(len(t), dtype=bool)
    for start, end in gaps:
        valid &= ~((t > start) & (t < end))
    return valid


def winding_transform(t, x, forcing, m_grid, t0_grid, sigma, valid=None):
    """Returns complex G, shape (len(m_grid), len(t0_grid)), a boolean
    edge_mask (True where the window at that t0 sits within one sigma of
    either end of the record -- analogous to a wavelet's cone of
    influence), and a boolean gap_mask (True where t0 itself falls inside a
    real data gap -- e.g. brestexcl's 1944-1954 -- i.e. there is no genuine
    sample anywhere near that window center, only fabricated
    linear-interpolation filler; unlike edge_mask this is "no real data
    here at all" rather than "reduced confidence," so it's blanked outright
    rather than just hatched). `valid`: boolean array over t, True at
    genuine (non-fabricated) samples; None means all points are genuine.
    `valid` also excludes fabricated points from the weighted sum itself,
    so a window straddling a gap is computed from real data only -- the
    transform is a local, non-autonomous weighted sum, not a global
    transform requiring uniform sampling, so this doesn't need any other
    special handling."""
    n_m, n_t0 = len(m_grid), len(t0_grid)
    G = np.zeros((n_m, n_t0), dtype=complex)
    edge_mask = np.zeros(n_t0, dtype=bool)
    gap_mask = np.zeros(n_t0, dtype=bool)
    t_min, t_max = t.min(), t.max()
    if valid is None:
        valid = np.ones(len(t), dtype=bool)
    for j, t0 in enumerate(t0_grid):
        nearest = np.argmin(np.abs(t - t0))
        gap_mask[j] = not valid[nearest]
        w = np.exp(-0.5 * ((t - t0) / sigma) ** 2) * valid
        wsum = w.sum()
        if wsum <= 0:
            continue
        wx = w * x
        basis = np.exp(-1j * 2.0 * np.pi * np.outer(m_grid, forcing))
        G[:, j] = (basis @ wx) / wsum
        edge_mask[j] = (t0 - t_min < sigma) or (t_max - t0 < sigma)
    return G, edge_mask, gap_mask


def noise_floor(year, forcing, m_grid, t0_grid, sigma, n_reps=32, seed=0,
                 valid=None, rho=0.97):
    """The winding basis exp(-i*2*pi*M*Forcing(t)) is not a clean orthogonal
    basis in M — Forcing's own recurrence structure makes some M values
    intrinsically more "resonant" than others for *any* input, independent
    of signal content. This estimates that per-M leakage floor from
    matched-length surrogates, so it can be divided out, leaving only power
    that exceeds what noise alone would produce at that M. `valid` is
    threaded through so the floor is estimated the same way the real
    transform is -- excluding fabricated filler across a real data gap, if
    any.

    Surrogates are AR(1) red noise (rho=0.97 by default), matching
    winding_rank.ar1_floor's convention, NOT white noise. Real climate/
    geophysical series are strongly autocorrelated month-to-month, so they
    carry far more low-frequency power than white noise does with zero
    forced signal required -- a white-noise floor is anti-conservative for
    this kind of data (confirmed directly: pure AR1 red noise with no real
    signal produced a false "bright ridge" of ~4 bits against a white-noise
    floor, vs. ~1.2 bits -- correctly unremarkable -- against this one, at
    the same M). Use a lower rho (or plug in the series' own measured lag-1
    autocorrelation) if the target index is known to be less persistent
    than typical monthly SST."""
    rng = np.random.default_rng(seed)
    floor = np.zeros(len(m_grid))
    s = np.sqrt(max(1e-12, 1.0 - rho * rho))
    n = len(year)
    for _ in range(n_reps):
        e = rng.standard_normal(n)
        z = np.zeros(n)
        for i in range(1, n):
            z[i] = rho * z[i - 1] + s * e[i]
        noise = standardize(z)
        Gn, _, _ = winding_transform(year, noise, forcing, m_grid, t0_grid,
                                      sigma, valid=valid)
        floor += np.mean(np.abs(Gn) ** 2, axis=1)
    return floor / n_reps


def compute_log_power(G, floor=None):
    """Raw log2(power) by default -- floor is no longer divided out here.

    Dividing power by a per-M floor before display reshapes the very
    thing being looked at: the floor itself has real M-dependent
    structure (the winding basis's own non-orthogonality), so dividing
    by it can distort relative-amplitude comparisons ACROSS M, which is
    exactly what this tool is for (see spectrum_panel for the intended
    role of the floor -- a reference curve alongside the signal, the way
    an AR1 line sits alongside a power spectrum, not a divisor baked into
    it). `floor` is accepted and ignored if passed, kept only so old call
    sites don't break."""
    power = np.abs(G) ** 2
    return np.log2(np.maximum(power, 1e-12))


def robust_scale(log_power):
    # Robust (percentile-based) color scale: a handful of extreme outlier
    # bins otherwise wash out the visible contrast in the bulk of the map.
    vmin, vmax = np.percentile(log_power, [3, 99.5])
    return vmin, vmax


def spectrum_panel(ax, m_grid, log_power_obs, log_power_model, floor,
                    fitted_m):
    """The floor's actual role: a 1D reference curve against M (like an
    AR1 dashed line on a power spectrum), plotted alongside the
    time-averaged raw power -- not divided into the 2D heatmaps above."""
    mean_obs = np.log2(np.maximum(np.mean(2.0 ** log_power_obs, axis=1), 1e-12))
    mean_model = np.log2(np.maximum(np.mean(2.0 ** log_power_model, axis=1), 1e-12))
    log_floor = np.log2(np.maximum(floor, 1e-12))
    ax.plot(mean_obs, m_grid, color="tab:blue", lw=1.0, label="Data (time-mean)")
    ax.plot(mean_model, m_grid, color="tab:orange", lw=1.0, label="Model (time-mean)")
    ax.plot(log_floor, m_grid, color="0.4", lw=1.0, linestyle="--",
             label="AR1 red-noise floor")
    if fitted_m is not None:
        for m in fitted_m:
            ax.axhline(abs(m), color="red", linestyle=":", linewidth=0.6,
                       alpha=0.6)
    ax.set_xlabel("log2 power")
    ax.legend(fontsize=6, loc="upper right")
    ax.set_title("mean spectrum vs. AR1 floor", fontsize=8)


def plot_panel(ax, t0_grid, m_grid, log_power, edge_mask, fitted_m, title,
               vmin, vmax, cmap="viridis", gap_mask=None):
    X, Y = np.meshgrid(t0_grid, m_grid)
    pcm = ax.pcolormesh(X, Y, log_power, shading="auto", cmap=cmap,
                         vmin=vmin, vmax=vmax)

    for j in np.where(edge_mask)[0]:
        ax.axvspan(t0_grid[j] - (t0_grid[1] - t0_grid[0]) / 2,
                   t0_grid[j] + (t0_grid[1] - t0_grid[0]) / 2,
                   color="white", alpha=0.35, hatch="//", linewidth=0)

    if gap_mask is not None:
        # Solid (opaque) blank, distinct from edge_mask's translucent
        # diagonal hatch -- this is a real absence of data (e.g.
        # brestexcl's 1944-1954 gap), not just reduced-confidence
        # near-edge coverage, so it reads as genuinely blank, not "partly
        # trustworthy."
        for j in np.where(gap_mask)[0]:
            ax.axvspan(t0_grid[j] - (t0_grid[1] - t0_grid[0]) / 2,
                       t0_grid[j] + (t0_grid[1] - t0_grid[0]) / 2,
                       color="white", alpha=1.0, linewidth=0, zorder=5)

    if fitted_m is not None:
        for m in fitted_m:
            ax.axhline(abs(m), color="red", linestyle=":", linewidth=0.8,  # "--"
                       alpha=0.8)
            ax.text(t0_grid[-1], abs(m), f" {abs(m):.3f}", color="red",  # "white"
                    fontsize=6, va="center", ha="left")

    ax.set_title(title, fontweight="bold")
    ax.set_ylabel("winding number  M  (dimensionless)")
    return pcm


def make_winding_scalogram(idx: str, m_max: float, dm: float, sigma: float,
                            t0_step: float, outdir: Path | None,
                            cmap: str = "viridis", stacked: bool = False):
    csv_path = ROOT / idx / "lte_results.csv"
    if not csv_path.exists():
        print(f"  (skipping {idx}: no lte_results.csv)", file=sys.stderr)
        return

    year, dt, model, obs, forcing = load_columns(idx)
    ir = load_ir(idx)
    model = uncompensate_yearly_feedback(model, ir)
    obs_s = standardize(obs)
    model_s = standardize(model)

    m_grid = np.arange(0.0, m_max + dm / 2, dm)
    t0_grid = np.arange(year[0] + sigma / 2, year[-1] - sigma / 2 + 1e-9,
                         t0_step)
    if len(t0_grid) < 2:
        print(f"  (skipping {idx}: record too short for sigma={sigma})",
              file=sys.stderr)
        return

    gaps = detect_gaps(idx)
    valid = valid_mask_from_gaps(year, gaps)
    if gaps:
        gap_desc = ", ".join(f"{s:.2f}-{e:.2f}" for s, e in gaps)
        print(f"  [{idx}] excluding fabricated filler across real gap(s): "
              f"{gap_desc}")

    G_obs, edge_mask, gap_mask_obs = winding_transform(
        year, obs_s, forcing, m_grid, t0_grid, sigma, valid=valid)
    G_model, _, gap_mask_model = winding_transform(
        year, model_s, forcing, m_grid, t0_grid, sigma, valid=valid)
    gap_mask = gap_mask_obs | gap_mask_model
    floor = noise_floor(year, forcing, m_grid, t0_grid, sigma, valid=valid)

    w = compute_winding(idx)
    fitted_m = np.abs(w["m"]) if w is not None else None

    log_power_obs = compute_log_power(G_obs)
    log_power_model = compute_log_power(G_model)

    title_obs = f"{idx}: Data (col 3) winding against Forcing (col 4)"
    title_model = (f"{idx}: Model (col 2, IR-corrected) winding against "
                    f"Forcing (col 4)")

    if stacked:
        fig = plt.figure(figsize=(13, 9))
        gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], hspace=0.3, wspace=0.15)
        ax_spec = fig.add_subplot(gs[:, 0])
        axes = (fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1]))
        axes[0].sharex(axes[1])
        vmin_obs, vmax_obs = robust_scale(log_power_obs)
        vmin_model, vmax_model = robust_scale(log_power_model)
        pcm1 = plot_panel(axes[0], t0_grid, m_grid, log_power_obs, edge_mask,
                           fitted_m, title_obs, vmin_obs, vmax_obs, cmap=cmap,
                           gap_mask=gap_mask)
        pcm2 = plot_panel(axes[1], t0_grid, m_grid, log_power_model,
                           edge_mask, fitted_m, title_model, vmin_model,
                           vmax_model, cmap=cmap, gap_mask=gap_mask)
        spectrum_panel(ax_spec, m_grid, log_power_obs, log_power_model,
                       floor, fitted_m)
        ax_spec.set_ylabel("winding number  M  (dimensionless)")
        axes[1].set_xlabel("year (window center)")
        for pcm, ax in ((pcm1, axes[0]), (pcm2, axes[1])):
            cb = fig.colorbar(pcm, ax=ax, pad=0.01)
            cb.set_label("log2 power (raw)")
    else:
        fig, (ax_spec, axes0, axes1) = plt.subplots(
            1, 3, figsize=(19, 6.5), sharey=True, layout="constrained",
            gridspec_kw={"width_ratios": [1, 2.5, 2.5]})
        axes = (axes0, axes1)
        vmin, vmax = robust_scale(np.concatenate([log_power_obs.ravel(),
                                                   log_power_model.ravel()]))
        pcm1 = plot_panel(axes[0], t0_grid, m_grid, log_power_obs, edge_mask,
                           fitted_m, title_obs, vmin, vmax, cmap=cmap,
                           gap_mask=gap_mask)
        pcm2 = plot_panel(axes[1], t0_grid, m_grid, log_power_model,
                           edge_mask, fitted_m, title_model, vmin, vmax,
                           cmap=cmap, gap_mask=gap_mask)
        spectrum_panel(ax_spec, m_grid, log_power_obs, log_power_model,
                       floor, fitted_m)
        ax_spec.set_ylabel("winding number  M  (dimensionless)")
        axes[0].set_xlabel("year (window center)")
        axes[1].set_xlabel("year (window center)")
        axes[1].set_ylabel("")
        cb = fig.colorbar(pcm2, ax=list(axes), pad=0.01, shrink=0.9)
        cb.set_label("log2 power (raw)")

    fig.suptitle(
        f"Winding scalograms — {idx}  (window sigma={sigma:g}yr; dashed "
        f"lines = this index's own fitted winding numbers; hatched = "
        f"within one sigma of record edge; color = RAW log2 power, "
        f"undivided -- see left panel for the AR1 red-noise floor as a "
        f"reference curve, not a divisor)",
        fontsize=9)
    if stacked:
        fig.tight_layout(rect=(0, 0, 1, 0.96))

    out_dir = outdir if outdir is not None else (ROOT / idx)
    suffix = "" if cmap == "viridis" else f"_{cmap}"
    out_path = out_dir / f"winding_scalogram{suffix}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("indices", nargs="*", default=DEFAULT_INDICES)
    ap.add_argument("--m-max", type=float, default=5.0,
                     help="max winding number shown (default 5.0)")
    ap.add_argument("--dm", type=float, default=0.01,
                     help="winding-number grid resolution (default 0.01)")
    ap.add_argument("--sigma", type=float, default=15.0,
                     help="time-window half-width, years (default 15)")
    ap.add_argument("--t0-step", type=float, default=5.0,
                     help="spacing between window centers, years (default 5)")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--cmap", choices=["viridis", "hot"], default="viridis",
                     help="'hot' = black (weak) through red to white-hot "
                          "(strong) instead of viridis (default viridis)")
    ap.add_argument("--stacked", action="store_true",
                     help="original layout: Data/Model stacked vertically, "
                          "each with its own colorbar (default: side by "
                          "side, sharing one colorbar)")
    args = ap.parse_args()

    for idx in args.indices:
        print(f"-- {idx} --")
        make_winding_scalogram(idx, args.m_max, args.dm, args.sigma,
                                args.t0_step, args.outdir, cmap=args.cmap,
                                stacked=args.stacked)


if __name__ == "__main__":
    main()
