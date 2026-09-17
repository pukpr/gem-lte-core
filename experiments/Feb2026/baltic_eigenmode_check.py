#!/usr/bin/env python3
"""baltic_eigenmode_check.py — the cheap check before committing to a full
fluid-dynamics run: (1) translate the empirical ridge M* found by
cv_ridge_transfer.py into an actual calendar period for Baltic's own
forcing trajectory, then (2) compute the Baltic basin's own free linear
oscillation (seiche) periods -- barotropic AND baroclinic -- from a
lumped multi-basin "box network" model, and see whether any of them land
near the translated period.

Part 1 -- translation
----------------------
winding_scalogram.py's own docstring is explicit that M is conjugate to
Forcing itself, not calendar time: a term amp*sin(2*pi*M*Forcing(t)+phase)
winds at instantaneous phase rate d/dt[2*pi*M*Forcing(t)] = 2*pi*M*F'(t),
which is NOT constant, because Forcing(t) is a nonlinear, IIR-integrated
composite of many tidal terms, not a single sinusoid. So "M*=1.25" has no
single calendar period -- it has a *distribution* of instantaneous periods
1/(M*|F'(t)|), which this script reports two ways:
  - pointwise percentiles of the local instantaneous period (dominated by
    the broad, slowly-varying stretches between bursts of fast winding)
  - a total-variation-weighted average period: record span / (M * total
    unsigned variation of Forcing / 2*pi) -- how long, ON AVERAGE, a full
    2*pi revolution takes if you weight by how much winding actually
    happens (dominated instead by the bursts, since that's where most of
    the phase gets used up)
These can differ by an order of magnitude or more for a bursty series --
that gap is itself informative about how "unclocklike" this manifold is,
and both numbers matter for judging what an eigenmode comparison should
even be compared against.

Part 2 -- box-network eigenmodes
----------------------------------
A lumped-parameter ("hydraulic network") model: each major Baltic
sub-basin is a box holding a spatially-uniform sea-level (or interface)
anomaly zeta_i over area A_i; each connecting strait carries a volume flux
Q_ij with its own inertia, linearized as

    A_i * d(zeta_i)/dt = sum_j Q_ij                  (mass balance)
    d(Q_ij)/dt = (g * a_ij / L_ij) * (zeta_i - zeta_j)   (momentum, "inductor")

Combining gives A_i * d2(zeta_i)/dt2 = sum_j (g*a_ij/L_ij)*(zeta_j-zeta_i),
a generalized eigenvalue problem  K v = omega^2 M v  with mass matrix
M=diag(A_i) and a graph-Laplacian-like stiffness matrix K built from each
strait's conductance g*a_ij/L_ij. This is the same reduction Wubber &
Krauss-style seiche analyses use, just with hand-lumped boxes instead of a
full 2D grid -- appropriate for a first-pass "is this even in the right
ballpark" check, not a survey-grade computation.

The barotropic run uses full gravity g and total depth; the baroclinic run
substitutes reduced gravity g' = g*delta_rho/rho and the upper-layer
(halocline) thickness, giving the much slower internal-seiche family that
is the more plausible candidate for a monthly-to-multiyear resonance.

Basin/strait geometry below is hand-estimated from general knowledge of
Baltic Sea geography (areas, mean depths, sill depths) -- NOT a bathymetry
survey product. Treat absolute periods as order-of-magnitude, and rerun
with --area/--depth overrides (see --list-basins) to see how sensitive the
conclusion is to these numbers before trusting it further.

Usage
-----
    ./baltic_eigenmode_check.py
    ./baltic_eigenmode_check.py --m-star 1.25 --drho-rho 0.005 --h1 50
    ./baltic_eigenmode_check.py --list-basins
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import scipy.linalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare  # noqa: E402

ROOT = Path(__file__).resolve().parent
G = 9.81

# ---------------------------------------------------------------------------
# Part 1: translate M* into a calendar period
# ---------------------------------------------------------------------------

def translate_ridge(idx: str, m_star: float) -> dict:
    prep = prepare(idx)
    dates, forcing = prep["dates"], prep["forcing"]
    dfdt = np.gradient(forcing, dates)
    with np.errstate(divide="ignore"):
        inst_period = 1.0 / np.abs(m_star * dfdt)
    finite = np.isfinite(inst_period) & (inst_period < 1000.0)
    pct = {q: float(np.percentile(inst_period[finite], q))
           for q in (5, 25, 50, 75, 95)}

    total_time = float(dates[-1] - dates[0])
    total_variation = float(np.sum(np.abs(np.diff(forcing))))
    # one full cycle of sin(2*pi*M*F(t)) is one unit of change in M*F(t)
    # (not in the 2*pi*M*F(t) argument itself), so the number of unsigned
    # windings over the record is M*total_variation, not that over 2*pi.
    tv_period = total_time / (m_star * total_variation) \
        if total_variation > 0 else float("inf")

    return dict(idx=idx, m_star=m_star, pct_years=pct,
                tv_period_years=tv_period, total_time=total_time)


def print_translation(t: dict) -> None:
    print(f"-- translating {t['idx']}'s M*={t['m_star']:.3f} into calendar time --")
    print(f"  record span: {t['total_time']:.1f} yr")
    print("  pointwise instantaneous period 1/(M*|dF/dt|), percentiles:")
    for q, v in t["pct_years"].items():
        print(f"    p{q:2d}: {v:8.3f} yr  = {v*365.25:9.1f} days")
    tv = t["tv_period_years"]
    print(f"  total-variation-weighted average period: {tv:.4f} yr "
          f"= {tv*365.25:.1f} days = {tv*12:.2f} months")
    print("  (these differ a lot because Forcing winds in bursts, not at a "
          "steady rate -- see this script's docstring)")


# ---------------------------------------------------------------------------
# Part 2: box-network eigenmodes
# ---------------------------------------------------------------------------

# name: (area_km2, depth_m)
BASINS = {
    "bothnian_bay":   (36_000, 40.0),
    "bothnian_sea":   (78_000, 60.0),
    "baltic_proper":  (211_000, 60.0),
    "gulf_of_finland": (30_000, 40.0),
    "gulf_of_riga":   (18_000, 25.0),
    "danish_straits": (8_000, 15.0),
}

# (basin_a, basin_b, width_km, sill_depth_m, length_km)
STRAITS = [
    ("bothnian_bay", "bothnian_sea", 80.0, 25.0, 100.0),      # Kvarken
    ("bothnian_sea", "baltic_proper", 100.0, 60.0, 150.0),    # Aland Sea (lumped)
    ("baltic_proper", "gulf_of_finland", 70.0, 50.0, 50.0),   # Hango/Osmussaar mouth
    ("baltic_proper", "gulf_of_riga", 27.0, 25.0, 30.0),      # Irbe Strait
    ("baltic_proper", "danish_straits", 15.0, 15.0, 60.0),    # Darss/Drogden sills
]
# open (grounded) connection to the North Sea/Kattegat, treated as an
# infinite reservoir at zeta=0 -- not a dynamical box of its own.
OPEN_CONNECTION = ("danish_straits", 30.0, 20.0, 80.0)


def build_network(basins: dict, straits: list, open_conn: tuple,
                  depth_scale: float = 1.0, gravity: float = G,
                  h1_cap: float | None = None) -> tuple[np.ndarray, np.ndarray, list]:
    """Mass matrix M=diag(area) and stiffness K (weighted graph Laplacian
    of strait conductances g*a/L), for either the barotropic case
    (depth_scale=1, gravity=g, h1_cap=None -> full depth) or the
    baroclinic case (gravity=g', h1_cap=upper-layer thickness -> each
    strait's active depth is capped at h1_cap, since only the upper layer
    participates in the internal-seiche mode, and a sill shallower than
    h1_cap is fully within the upper layer there)."""
    names = list(basins)
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    area_m2 = np.array([basins[nm][0] * 1e6 for nm in names])
    M = np.diag(area_m2)
    K = np.zeros((n, n))

    def active_depth(sill_depth_m: float) -> float:
        d = sill_depth_m * depth_scale
        return min(d, h1_cap) if h1_cap is not None else d

    for a, b, width_km, sill_depth_m, length_km in straits:
        depth = active_depth(sill_depth_m)
        cross_section = width_km * 1e3 * depth
        conductance = gravity * cross_section / (length_km * 1e3)
        i, j = idx[a], idx[b]
        K[i, i] += conductance
        K[j, j] += conductance
        K[i, j] -= conductance
        K[j, i] -= conductance

    open_basin, width_km, sill_depth_m, length_km = open_conn
    depth = active_depth(sill_depth_m)
    cross_section = width_km * 1e3 * depth
    conductance = gravity * cross_section / (length_km * 1e3)
    K[idx[open_basin], idx[open_basin]] += conductance

    return M, K, names


def eigenperiods(M: np.ndarray, K: np.ndarray) -> np.ndarray:
    """omega^2 from the generalized eigenproblem K v = omega^2 M v ->
    periods 2*pi/omega, sorted ascending (fastest mode first)."""
    w2 = scipy.linalg.eigh(K, M, eigvals_only=True)
    w2 = np.clip(w2, 0.0, None)
    omega = np.sqrt(w2)
    periods = np.full_like(omega, np.inf)
    nz = omega > 1e-12
    periods[nz] = 2 * np.pi / omega[nz]
    return np.sort(periods)


def run(idx: str, m_star: float, drho_rho: float, h1: float,
        outdir: Path | None) -> None:
    t = translate_ridge(idx, m_star)
    print_translation(t)

    M_bt, K_bt, names = build_network(BASINS, STRAITS, OPEN_CONNECTION,
                                       gravity=G)
    periods_bt_s = eigenperiods(M_bt, K_bt)
    periods_bt_hr = periods_bt_s / 3600.0

    g_reduced = G * drho_rho
    M_bc, K_bc, _ = build_network(BASINS, STRAITS, OPEN_CONNECTION,
                                  gravity=g_reduced, h1_cap=h1)
    periods_bc_s = eigenperiods(M_bc, K_bc)
    periods_bc_days = periods_bc_s / 86400.0

    n_modes = len(names)
    print(f"\n-- {idx}: box-network free-oscillation periods "
          f"({n_modes} basins -> {n_modes} modes) --")
    print(f"  barotropic (g={G}, full depth), fastest {n_modes} modes:")
    for p in periods_bt_hr[:n_modes]:
        print(f"    {p:8.2f} hours" if np.isfinite(p) else "    (zero mode)")
    print(f"  baroclinic (g'={g_reduced:.4f}, upper layer capped at "
          f"{h1:.0f} m, drho/rho={drho_rho:g}), fastest {n_modes} modes:")
    for p in periods_bc_days[:n_modes]:
        print(f"    {p:8.2f} days  ({p/30.44:5.2f} months)"
              if np.isfinite(p) else "    (zero mode)")

    tv_days = t["tv_period_years"] * 365.25
    p25_days = t["pct_years"][25] * 365.25
    p75_days = t["pct_years"][75] * 365.25
    slowest_bt_days = periods_bt_hr[-1] / 24.0
    fastest_bc_days = periods_bc_days[0]
    slowest_bc_days = periods_bc_days[-1]

    def overlaps(lo: float, hi: float, target: float) -> bool:
        return lo <= target <= hi

    print(f"\n-- verdict: two different readings of 'the' target period, "
          f"checked separately --")
    print(f"  (a) percentile-based target (bulk of pointwise instantaneous "
          f"periods): ~{p25_days:.0f}-{p75_days:.0f} days")
    print(f"      vs. barotropic (up to {slowest_bt_days:.1f} d): "
          f"{'OVERLAPS' if overlaps(0, slowest_bt_days, p25_days) else 'no overlap -- barotropic is far too fast'}")
    print(f"      vs. baroclinic ({fastest_bc_days:.1f}-{slowest_bc_days:.1f} d): "
          f"{'OVERLAPS' if (fastest_bc_days <= p75_days and slowest_bc_days >= p25_days) else 'no overlap -- even the slowest baroclinic mode is too fast'}")
    print(f"  (b) total-variation-weighted target (winds in bursts; this "
          f"is the statistic closer in spirit to what winding_scalogram.py's "
          f"own windowed transform actually measures): {tv_days:.0f} days")
    print(f"      vs. barotropic (up to {slowest_bt_days:.1f} d): "
          f"{'OVERLAPS' if overlaps(0, slowest_bt_days, tv_days) else 'no overlap -- barotropic is too fast'}")
    print(f"      vs. baroclinic ({fastest_bc_days:.1f}-{slowest_bc_days:.1f} d): "
          f"{'OVERLAPS -- a real internal-seiche mode is a live candidate' if overlaps(fastest_bc_days, slowest_bc_days, tv_days) else 'no overlap'}")
    print("  reading (a) says no free linear mode explains this; reading "
          "(b) says a slow baroclinic mode plausibly could. Which reading "
          "is right is itself an open question this script can't settle --"
          " see the docstring.")

    plot_comparison(idx, t, periods_bt_hr, periods_bc_days,
                    outdir if outdir is not None else ROOT / idx)


def plot_comparison(idx: str, t: dict, periods_bt_hr: np.ndarray,
                    periods_bc_days: np.ndarray, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.set_xscale("log")

    bt_days = periods_bt_hr[np.isfinite(periods_bt_hr)] / 24.0
    bc_days = periods_bc_days[np.isfinite(periods_bc_days)]
    ax.scatter(bt_days, np.full(len(bt_days), 2), marker="|", s=400,
               color="tab:blue", label="barotropic modes")
    ax.scatter(bc_days, np.full(len(bc_days), 1), marker="|", s=400,
               color="tab:orange", label="baroclinic modes")

    p5, p95 = t["pct_years"][5] * 365.25, t["pct_years"][95] * 365.25
    p25, p75 = t["pct_years"][25] * 365.25, t["pct_years"][75] * 365.25
    ax.axvspan(p5, p95, color="tab:green", alpha=0.10)
    ax.axvspan(p25, p75, color="tab:green", alpha=0.25,
               label="translated M* target band (p25-p75, p5-p95 shaded)")
    ax.axvline(t["tv_period_years"] * 365.25, color="darkgreen",
               linestyle="--", linewidth=1.3,
               label="total-variation average period")

    ax.set_yticks([1, 2])
    ax.set_yticklabels(["baroclinic\n(internal seiche)",
                        "barotropic\n(surface seiche)"])
    ax.set_ylim(0.5, 2.5)
    ax.set_xlabel("period (days, log scale)")
    ax.set_title(
        f"{idx}: does any free linear basin mode land near the "
        f"M*={t['m_star']:.3f}-implied period?", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_path = out_dir / "baltic_eigenmode_check.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def list_basins() -> None:
    print("basins: name  area_km2  depth_m")
    for name, (area, depth) in BASINS.items():
        print(f"  {name:16s} {area:8.0f}  {depth:5.1f}")
    print("straits: a <-> b  width_km  sill_depth_m  length_km")
    for a, b, w, d, ln in STRAITS:
        print(f"  {a:16s} <-> {b:16s} {w:6.1f} {d:6.1f} {ln:6.1f}")
    a, w, d, ln = OPEN_CONNECTION
    print(f"open: {a} <-> North Sea/Kattegat (ground)  {w:.1f}km wide "
          f"{d:.1f}m deep {ln:.1f}km long")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indices", nargs="*", default=["baltic"])
    ap.add_argument("--m-star", type=float, default=1.250,
                     help="ridge location to translate (default: the "
                          "cv_ridge_transfer.py baltic finding, 1.250)")
    ap.add_argument("--drho-rho", type=float, default=0.005,
                     help="fractional density contrast across the Baltic "
                          "halocline (default 0.005, typical order of "
                          "magnitude for this basin)")
    ap.add_argument("--h1", type=float, default=50.0,
                     help="upper-layer (halocline) thickness in meters "
                          "for the baroclinic network (default 50)")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--list-basins", action="store_true",
                     help="print the basin/strait table and exit")
    args = ap.parse_args()

    if args.list_basins:
        list_basins()
        return 0

    for idx in args.indices:
        run(idx, args.m_star, args.drho_rho, args.h1, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
