#!/usr/bin/env python3
"""geospatial_winding_analysis.py -- for every subdirectory containing a
lt.exe.windings.json (this project's own post-run winding-table output,
added earlier this session: GEM.LTE.Primitives.Shared.Save_Windings),
test whether geographically CLOSE indices share more of their winding
numbers (matching M, coherent phase) than geographically REMOTE ones.

METHODOLOGY CORRECTION (made after the first version of this script):
comparing each index's own JSON phase directly against another index's
own JSON phase is unsound. Every index here has been through its own
independent Ada optimizer run since creation (confirmed: even LPAP
Doodson tables have drifted apart from any common donor, not just the
Bessel-step impA/impB/bg/offs), so each index's manifold is only ~99.6%
correlated with any other's overall -- and that 0.4% disagreement
decorrelates CATASTROPHICALLY once multiplied by a winding number M
(verified directly: kN_baltic vs kN_northsea's own manifolds correlate
0.9966 overall, but sin(2*pi*M*F) between them collapses to +0.31 at
M~0.21 and NEGATIVE at M~1.25 -- manufacturing a spurious "anti-phase"
reading that is not present in the real data at all, R=0.974 there).

Fix: for each PAIR (A, B), use ONE shared forcing -- A's own, rebuilt from
A's lt.exe.p/resp -- for BOTH regressions. A's own JSON amp/phase are
already expressed against A's own forcing (no refit needed). B's RAW data
(F9-filtered, from B's own .dat) is freshly regressed against A's forcing
(interpolated onto B's own date grid) at the M values that match A's
JSON M's within tolerance. This is the same fix already applied and
verified for Baltic MSL vs Kaplan SST at the very start of this session's
investigation, and for kN_baltic vs kN_northsea specifically this turn.

Coordinate sources: see the previous version's docstring (unchanged) --
exact grid-box centers for kN040_E030/kN040_E130/kS040_E130, a confirmed
point for kN_brest, standard basin/box centers elsewhere, several
indices (NAO, PNA, SAM, SOI, TPI) flagged as not point-like and excluded
from the distance analysis.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lte_forward import (read_resp, Config, doodson_periods, year_length, tide_sum,
                         impulse_delta, iir, bessel, filter9point)

ROOT = Path(__file__).resolve().parent

COORDS = {
    "amo":           (40.0, -40.0, True,  "basin-wide N. Atlantic, representative center"),
    "nao":           (50.0, -22.0, False, "Azores/Iceland dipole -- not point-like"),
    "pdo":           (40.0, -180.0, True, "N. Pacific poleward of 20N, representative center"),
    "pna":           (50.0, -160.0, False, "diffuse wave train, central Pacific-N.America -- not point-like"),
    "nino4":         (0.0, -175.0, True,  "NOAA Nino4 box 5S-5N,160E-150W center"),
    "nino34":        (0.0, -145.0, True,  "NOAA Nino3.4 box 5S-5N,170W-120W center"),
    "tna":           (15.0, -36.0, True,  "Tropical N. Atlantic box center"),
    "tpi":           (0.0, -180.0, False, "tripole index -- not point-like"),
    "sam":           (-55.0, 0.0, False,  "zonally symmetric around Antarctica -- not point-like"),
    "soi":           (15.0, 160.0, False, "Tahiti/Darwin dipole -- not point-like"),
    "iode":          (-5.0, 100.0, True,  "IOD east pole box 10S-0,90-110E center"),
    "kN_baltic":     (58.0, 20.0, True,   "Baltic Sea proper, standard center"),
    "kN_northsea":   (56.0, 3.0, True,    "North Sea proper, standard center"),
    "kN_brest":      (48.4, -4.5, True,   "Brest, France -- confirmed via loc.png"),
    "kN_indostraits":(-2.0, 120.0, True,  "Indonesian throughflow/straits, standard center"),
    "kN_okhotsk":    (55.0, 150.0, True,  "Sea of Okhotsk, standard center"),
    "kN040_E030":    (40.0, 30.0, True,   "exact grid-box center"),
    "kN040_E130":    (40.0, 130.0, True,  "exact grid-box center"),
    "kS040_E130":    (-40.0, 130.0, True, "exact grid-box center"),
}


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def load_windings(idx: str) -> list[tuple[float, float, float]]:
    d = json.loads((ROOT / idx / "lt.exe.windings.json").read_text())
    return [(row[0], row[1], row[2]) for row in d["k_amp_phase"]]


def build_index(idx: str) -> dict:
    """Rebuild this index's own manifold (Forcing) and its own F9-filtered
    data, from its own lt.exe.p/resp/dat -- everything needed to either
    BE the shared reference forcing for a pair, or be regressed against
    someone else's."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    resp = read_resp(idx_dir / "lt.exe.resp")
    resp["CLIMATE_INDEX"] = str(idx_dir / f"{idx}.dat")
    overrides = {}
    if "harm" in params:
        overrides["NH"] = " ".join(str(int(float(h))) for h in params["harm"])
    cfg = Config(resp, overrides)
    dat = np.loadtxt(idx_dir / f"{idx}.dat")
    dates, data0 = dat[:, 0].copy(), dat[:, 1].copy()
    lpap = np.array(params["lpap"], dtype=float)
    lt = np.array(params["ltep"], dtype=float)
    year_startup = cfg.get("YEAR", 0.0)
    year_cand = float(params.get("year", 0.0))
    sampling = int(cfg.get("SAMPLING", 12.0))
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
          "ann1", "ann2", "IR", "ma", "mp", "shfT", "init")}
    periods = doodson_periods(year_startup, year_cand)
    ap = lpap[:, 1:3].copy()
    yl = year_length(year_startup, year_cand)
    tf = tide_sum(dates, ap, periods, yl, scaling=0.0, integ=B["shfT"])
    forcing = iir(tf * impulse_delta(dates, B["delA"], B["delB"], B["asym"], sampling),
                  lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
                  start_date=cfg.get("IDATE", 0.0), dates=dates)
    forcing = bessel(forcing, B["impA"], B["impB"], lt[-1], B["offs"], B["bg"])
    f9 = cfg.get("F9", 1)
    data = data0.copy()
    if f9 > 0:
        for _ in range(f9):
            data = filter9point(data)
    return dict(dates=dates, forcing=forcing, data=data)


def fit_at_M(F: np.ndarray, y: np.ndarray, M: float) -> tuple[float, float, float]:
    design = np.column_stack([np.ones_like(F), np.sin(2 * np.pi * M * F), np.cos(2 * np.pi * M * F)])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    amp = math.hypot(coef[1], coef[2])
    phase = math.atan2(coef[1], coef[2])
    r = float(np.corrcoef(design @ coef, y)[0, 1]) if y.std() > 0 else 0.0
    return amp, phase, r


def match_by_M(wa: list, wb: list, tol: float = 0.05) -> list[tuple[int, int]]:
    """Greedy nearest-neighbor index-pairs (i into wa, j into wb) with
    |Ma-Mb|<tol."""
    used_b = set()
    pairs = []
    for i, (ma, _, _) in enumerate(wa):
        best_j, best_d = None, tol
        for j, (mb, _, _) in enumerate(wb):
            if j in used_b:
                continue
            d = abs(ma - mb)
            if d < best_d:
                best_d, best_j = d, j
        if best_j is not None:
            used_b.add(best_j)
            pairs.append((i, best_j))
    return pairs


def corrected_pair_score(a: str, b: str, built: dict, windings: dict, tol: float = 0.05) -> dict:
    """Refit BOTH A's and B's raw data freshly, via the SAME single-mode
    OLS, against ONE shared forcing (A's own, interpolated onto each
    index's own date grid) -- an internally consistent apples-to-apples
    comparison. (An earlier version of this function reused A's JSON
    amp/phase -- from the real, multi-term Ada regression -- directly for
    A while only single-mode-refitting B; that mixes two different
    fitting procedures and was verified to give a materially different,
    less trustworthy answer than refitting both sides the same simple
    way, which is what a direct hand check of kN_baltic/kN_northsea
    confirmed as correct.)"""
    wa, wb = windings[a], windings[b]
    pairs = match_by_M(wa, wb, tol)
    if not pairs:
        return dict(n_shared=0, shared_frac=0.0, mean_cos_dphi=np.nan, mean_amp_log_ratio=np.nan)

    ref = built[a]
    F_ref_on_a = ref["forcing"]
    F_ref_on_b = np.interp(built[b]["dates"], ref["dates"], ref["forcing"])
    data_a, data_b = built[a]["data"], built[b]["data"]

    cos_dphi, log_ratio = [], []
    for i, j in pairs:
        ma, _, _ = wa[i]
        amp_a, phase_a, _ = fit_at_M(F_ref_on_a, data_a, ma)
        amp_b, phase_b, _ = fit_at_M(F_ref_on_b, data_b, ma)  # A's own M exactly, for both
        dphi = (phase_a - phase_b + math.pi) % (2 * math.pi) - math.pi
        cos_dphi.append(math.cos(dphi))
        if amp_a > 0 and amp_b > 0:
            log_ratio.append(abs(math.log10(amp_a / amp_b)))
    denom = min(len(wa), len(wb))
    return dict(n_shared=len(pairs), shared_frac=len(pairs) / denom if denom else 0.0,
                mean_cos_dphi=float(np.mean(cos_dphi)),
                mean_amp_log_ratio=float(np.mean(log_ratio)) if log_ratio else np.nan)


def main() -> None:
    indices = sorted(d.name for d in ROOT.iterdir()
                     if d.is_dir() and (d / "lt.exe.windings.json").exists())
    point_indices = [i for i in indices if i in COORDS and COORDS[i][2]]
    print(f"[geospatial_winding_analysis] {len(point_indices)} point-like indices: {point_indices}")

    windings = {i: load_windings(i) for i in indices}
    print("[geospatial_winding_analysis] rebuilding each index's own manifold + F9 data "
          "(needed as the shared reference forcing for its pairs)...")
    built = {}
    for i in point_indices:
        try:
            built[i] = build_index(i)
        except Exception as e:
            print(f"  (skipping {i}: {e})")
    point_indices = [i for i in point_indices if i in built]

    rows = []
    for i, a in enumerate(point_indices):
        for b in point_indices[i + 1:]:
            lat_a, lon_a, _, _ = COORDS[a]
            lat_b, lon_b, _, _ = COORDS[b]
            dist = haversine_km(lat_a, lon_a, lat_b, lon_b)
            sc = corrected_pair_score(a, b, built, windings)
            rows.append(dict(a=a, b=b, dist_km=dist, **sc))

    rows_sorted = sorted(rows, key=lambda r: r["dist_km"])
    print(f"\n[geospatial_winding_analysis] {len(rows)} pairs, CORRECTED (shared-forcing) phase comparison, tol=0.05 in M")
    print(f"{'pair':35s} {'dist(km)':>9s} {'n_shared':>9s} {'shared_frac':>12s} {'mean_cos(dphi)':>15s} {'amp_log_ratio':>14s}")
    for r in rows_sorted:
        print(f"{r['a']+'-'+r['b']:35s} {r['dist_km']:9.0f} {r['n_shared']:9d} "
              f"{r['shared_frac']:12.3f} {r['mean_cos_dphi']:15.3f} {r['mean_amp_log_ratio']:14.3f}")

    dists = np.array([r["dist_km"] for r in rows])
    fracs = np.array([r["shared_frac"] for r in rows])
    cosd = np.array([r["mean_cos_dphi"] for r in rows])
    valid = ~np.isnan(cosd)

    r_dist_frac = np.corrcoef(dists, fracs)[0, 1]
    r_dist_cos = np.corrcoef(dists[valid], cosd[valid])[0, 1] if valid.sum() > 2 else np.nan
    print(f"\n[geospatial_winding_analysis] correlation(distance, shared_frac)          = {r_dist_frac:+.3f}")
    print(f"[geospatial_winding_analysis] correlation(distance, CORRECTED mean_cos_dphi) = {r_dist_cos:+.3f}")
    print(f"[geospatial_winding_analysis] mean CORRECTED cos(dphi) overall = {cosd[valid].mean():+.3f}  "
          f"(vs the UNCORRECTED version's, which was dominated by forcing-mismatch artifacts)")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].scatter(dists, fracs, alpha=0.7)
    axes[0].set_xlabel("great-circle distance (km)")
    axes[0].set_ylabel("shared-winding fraction")
    axes[0].set_title(f"Distance vs. shared windings  (r={r_dist_frac:+.3f})", fontsize=10)
    axes[1].scatter(dists[valid], cosd[valid], alpha=0.7, color="tab:green")
    axes[1].axhline(0, color="0.7", lw=0.8)
    axes[1].set_xlabel("great-circle distance (km)")
    axes[1].set_ylabel("mean cos(phase difference), CORRECTED (shared forcing)")
    axes[1].set_title(f"Distance vs. CORRECTED phase coherence  (r={r_dist_cos:+.3f})", fontsize=10)
    fig.suptitle("Geospatial structure of shared winding numbers (shared-forcing corrected)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = ROOT / "geospatial_winding_analysis_corrected.png"
    fig.savefig(out, dpi=140)
    print(f"\n[geospatial_winding_analysis] saved {out}")


if __name__ == "__main__":
    main()
