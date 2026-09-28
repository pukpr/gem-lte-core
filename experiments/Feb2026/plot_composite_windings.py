#!/usr/bin/env python3
"""plot_composite_windings.py -- plot a few cycles of the composite sinusoid
built from a region's fitted k_amp_phase triplets (as saved by Save_Windings
into lt.exe.windings.json), for two regions/fits side by side.

Each entry in "k_amp_phase" is a [k, amp, phase] triplet describing one term
    amp * sin(k * t + phase)
and the composite signal is just their sum. k is a rate in the model's own
forcing-manifold coordinate, not calendar time -- units don't matter here,
only the *shape* of the two composites relative to each other, so time is
plotted in units of the shared ~0.2075 "backbone" cycle (the one term nearly
every GEM-LTE index converges on) unless overridden.

Usage:
    plot_composite_windings.py <windings.json> <windings.json>
        [--cycles N] [--labels A B] [--out FILE.png]

Example (this project's Red Sea/Persian Gulf post-1950 vs full-record fits):
    plot_composite_windings.py kN020_E050/lt.exe.windings.json \\
        kN020_E050_/lt.exe.windings.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BACKBONE_K = 0.2075  # universal constant shared by nearly every index this project has fit


def load_k_amp_phase(path: Path) -> list[tuple[float, float, float]]:
    data = json.loads(path.read_text())
    return [(float(k), float(a), float(p)) for k, a, p in data["k_amp_phase"]]


def composite(t: np.ndarray, terms: list[tuple[float, float, float]]) -> np.ndarray:
    y = np.zeros_like(t)
    for k, amp, phase in terms:
        y += amp * np.sin(2.0 * np.pi * k * t + phase)
    return y


def backbone_k(terms: list[tuple[float, float, float]]) -> float:
    """The term nearest the shared ~0.2075 backbone frequency, used only to
    set a sensible, comparable time window -- falls back to the median k of
    the fit's own terms if nothing is close."""
    ks = [k for k, _, _ in terms]
    best = min(ks, key=lambda k: abs(k - BACKBONE_K))
    if abs(best - BACKBONE_K) < 0.05:
        return best
    return float(np.median(ks))


def find_best_shift(t: np.ndarray, y_a: np.ndarray,
                     terms_b: list[tuple[float, float, float]],
                     search_range: float, n_grid: int = 801) -> tuple[float, float]:
    """Coarse-then-fine grid search for the time shift tau (applied to B's
    OWN time coordinate, i.e. B is re-evaluated at t+tau) that maximizes
    corrcoef(y_a, composite(t+tau, terms_b)) over +/- search_range. B is a
    pure sum of sinusoids, defined for any t, so shifting just means
    re-evaluating it -- no resampling/edge-trimming needed."""

    def cc_at(tau: float) -> float:
        y_b_shifted = composite(t + tau, terms_b)
        return float(np.corrcoef(y_a, y_b_shifted)[0, 1])

    def refine(center: float, half_width: float) -> tuple[float, float]:
        taus = np.linspace(center - half_width, center + half_width, n_grid)
        ccs = [cc_at(tau) for tau in taus]
        i = int(np.argmax(ccs))
        return float(taus[i]), float(ccs[i])

    best_tau, best_cc = refine(0.0, search_range)
    # One refinement pass, zoomed into the coarse winner's own neighborhood.
    best_tau, best_cc = refine(best_tau, 2.0 * search_range / n_grid)
    return best_tau, best_cc


def label_for(path: Path) -> str:
    name = path.parent.name
    if name.endswith("_"):
        return f"{name} (full record)"
    return f"{name} (production)"


def run(path_a: Path, path_b: Path, cycles: float, labels: tuple[str, str] | None,
        points_per_cycle: int, out: Path, secondary: bool = False,
        flip: bool = False, align: bool = False) -> None:
    terms_a = load_k_amp_phase(path_a)
    terms_b = load_k_amp_phase(path_b)

    k_ref = 0.5 * (backbone_k(terms_a) + backbone_k(terms_b))
    period = 1.0 / k_ref
    t = np.linspace(0.0, cycles * period, int(cycles * points_per_cycle))

    y_a = composite(t, terms_a)
    y_b = composite(t, terms_b)

    if align and not secondary:
        print("--align only applies with --secondary; ignoring.")
        align = False

    tau = 0.0
    if align:
        cc_before = float(np.corrcoef(y_a, y_b)[0, 1])
        # "slide slightly": search within +/- one backbone period, not the
        # whole shown window.
        tau, cc_after = find_best_shift(t, y_a, terms_b, search_range=period)
        print(f"\n--align: shifted {path_b.parent.name} by tau={tau:+.4f} "
              f"time-units  (CC {cc_before:+.4f} -> {cc_after:+.4f})")
        y_b = composite(t + tau, terms_b)

    label_a, label_b = labels if labels else (label_for(path_a), label_for(path_b))

    print(f"reference backbone k = {k_ref:.5f}  (period = {period:.3f} time-units, "
          f"{cycles:g} cycles shown)")
    for name, terms in [(label_a, terms_a), (label_b, terms_b)]:
        top = sorted(terms, key=lambda t: -abs(t[1]))[:3]
        print(f"\n[{name}] {len(terms)} terms, top 3 by |amp|:")
        for k, amp, phase in top:
            print(f"  k={k:9.4f}  amp={amp:8.4f}  phase={phase:7.3f}")

    fig, (ax1) = plt.subplots(1, 1, figsize=(11, 3))

    line_a, = ax1.plot(t, y_a, color="tab:blue", lw=0.5, label=label_a)
    if secondary:
        # Independent y-scales: each series is plotted against its own
        # amplitude range, so relative SHAPE/phase can be compared even
        # when the two composites' absolute amplitudes differ a lot
        # (e.g. an SST index vs an MSL index).
        ax2 = ax1.twinx()
        line_b, = ax2.plot(t, y_b, color="tab:red", lw=0.5, label=label_b)
        ax1.set_ylabel(f"{label_a} amplitude (a.u.)", color="tab:blue")
        ax2.set_ylabel(f"{label_b} amplitude (a.u.)", color="tab:red")
        ax1.tick_params(axis="y", labelcolor="tab:blue")
        ax2.tick_params(axis="y", labelcolor="tab:red")
        ax1.legend(handles=[line_a, line_b], fontsize=9)
        ax2.axhline(0.0, color="0.7", lw=0.6)
        if flip:
            # Mirror the secondary axis through y=0 (not just reverse
            # top/bottom): a point drawn at +v now sits where -v used to
            # be, so an inverse/anti-correlated relationship between the
            # two series can be checked visually without negating B's
            # actual data.
            lo, hi = ax2.get_ylim()
            ax2.set_ylim(-hi, -lo)
        if align:
            # Top x-axis shows B's OWN native time coordinate: B is
            # plotted (at each bottom-axis position t) using its value at
            # t+tau, so the top axis -- aligned pixel-for-pixel with the
            # bottom one -- must read t+tau there too.
            ax_top = ax1.twiny()
            x_lo, x_hi = ax1.get_xlim()
            ax_top.set_xlim(x_lo + tau, x_hi + tau)
            ax_top.set_xlabel(f"{label_b} own manifold (shifted by {tau:+.4f})",
                               color="tab:red")
            ax_top.tick_params(axis="x", labelcolor="tab:red")
    else:
        ax1.plot(t, y_b, color="tab:red", lw=0.5, label=label_b)
        ax1.set_ylabel("composite amplitude (a.u.)")
        ax1.legend(fontsize=9)
        ax1.axhline(0.0, color="0.7", lw=0.6)
    ax1.set_title("Composite sinusoid comparison: sum of fitted $A_i sin(2\pi k_i M + \phi_i)$ terms over manifold range")

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\nsaved {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windings_a", type=Path, help="first lt.exe.windings.json")
    ap.add_argument("windings_b", type=Path, help="second lt.exe.windings.json")
    ap.add_argument("--cycles", type=float, default=10.0,
                     help="number of backbone cycles to show (default 6)")
    ap.add_argument("--labels", nargs=2, metavar=("A", "B"), default=None,
                     help="override the two legend labels (default: derived from directory names)")
    ap.add_argument("--points-per-cycle", type=int, default=300,
                     help="sample density (default 300/cycle)")
    ap.add_argument("--out", type=Path, default=Path("composite_windings_comparison.png"),
                     help="output plot path")
    ap.add_argument("--secondary", action="store_true",
                     help="plot the second series (B) on an independent "
                          "secondary y-axis instead of sharing A's scale "
                          "-- use when the two composites' amplitudes "
                          "aren't directly comparable")
    ap.add_argument("--flip", action="store_true",
                     help="flip the secondary y-axis about y=0 (mirrors B's "
                          "axis, not its data) -- only meaningful with "
                          "--secondary; use to visually check an "
                          "inverse/anti-correlated relationship")
    ap.add_argument("--align", action="store_true",
                     help="slide B in time (within +/- one backbone period) "
                          "to maximize its correlation with A, and add a "
                          "second x-axis on top showing B's own native time "
                          "at the shifted position -- only used with "
                          "--secondary")
    args = ap.parse_args()
    run(args.windings_a, args.windings_b, args.cycles, args.labels,
        args.points_per_cycle, args.out, args.secondary, args.flip, args.align)


if __name__ == "__main__":
    main()
