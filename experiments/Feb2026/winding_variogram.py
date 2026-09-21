#!/usr/bin/env python3
"""winding_variogram.py -- systematic classification of windings across
every index with a lt.exe.windings.json, and an empirical VARIOGRAM
(coherence vs. distance) split by harmonic-family membership.

Motivation (per this project's own framing): the winding triplets
(M, Amplitude, Phase) are the ONLY thing that distinguishes one climate
index from another in this model -- the forcing manifold itself is
shared. If winding coherence between sites varies smoothly with distance
(the way any spatial field needs to for kriging/interpolation to make
sense), that's the empirical variogram a real kriging scheme would be
built on. This script computes that directly, using the shared-forcing
-corrected phase methodology already validated in
geospatial_winding_analysis.py (naive same-file phase comparison was
shown this session to manufacture spurious antiphase readings).

Classification: a winding is "harmonic-family" if its M is within
tolerance of an integer multiple of the shared backbone M(NM)~0.2075
(the same backbone independently confirmed across 7+ indices in
WINDING_NARRATIVE.md) -- "non-harmonic" otherwise. The hypothesis under
test: harmonic-family windings behave like a real spatial field
(coherence decaying smoothly with distance); non-harmonic windings do
not (arbitrary local/noisy relationships regardless of distance).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geospatial_winding_analysis import (build_index, load_windings, match_by_M,
                                          fit_at_M, haversine_km, COORDS, ROOT)

BACKBONE = 0.2075
HARMONIC_TOL = 0.05


def is_harmonic(m: float) -> tuple[bool, int]:
    n = m / BACKBONE
    rn = round(n)
    return abs(n - rn) < HARMONIC_TOL and rn != 0, rn


def main() -> None:
    indices = sorted(d.name for d in ROOT.iterdir()
                     if d.is_dir() and (d / "lt.exe.windings.json").exists())
    point_indices = [i for i in indices if i in COORDS and COORDS[i][2]]
    windings = {i: load_windings(i) for i in indices}

    print(f"[winding_variogram] {len(point_indices)} point-like indices, building manifolds...")
    built = {}
    for i in point_indices:
        try:
            built[i] = build_index(i)
        except Exception as e:
            print(f"  (skipping {i}: {e})")
    point_indices = [i for i in point_indices if i in built]

    # -------------------------------------------------------------
    # Step 1: classify every DISTINCT winding across the whole dataset
    # -------------------------------------------------------------
    print("\n[winding_variogram] classifying all windings across all indices "
          f"(harmonic tolerance {HARMONIC_TOL} of backbone {BACKBONE}):")
    all_m = []
    for idx, w in windings.items():
        for m, amp, ph in w:
            harm, n = is_harmonic(m)
            all_m.append((idx, m, amp, harm, n))
    harm_count = sum(1 for *_, h, _ in all_m if h)
    print(f"  {len(all_m)} total winding instances across {len(windings)} indices: "
          f"{harm_count} harmonic-family ({harm_count/len(all_m):.0%}), "
          f"{len(all_m)-harm_count} non-harmonic")

    # how many DISTINCT indices does each harmonic order appear in?
    from collections import defaultdict
    by_order = defaultdict(set)
    for idx, m, amp, harm, n in all_m:
        if harm:
            by_order[n].add(idx)
    print("\n  harmonic order -> n indices carrying it (of", len(windings), "):")
    for n in sorted(by_order):
        print(f"    order {n:+3d} (M~{n*BACKBONE:+.3f}): {len(by_order[n]):2d} indices")

    # -------------------------------------------------------------
    # Step 2: pairwise, per-matched-winding coherence vs distance
    # -------------------------------------------------------------
    records = []  # (dist, M, is_harm, cos_dphi, amp_a, amp_b)
    for i, a in enumerate(point_indices):
        for b in point_indices[i + 1:]:
            wa, wb = windings[a], windings[b]
            pairs = match_by_M(wa, wb, tol=0.05)
            if not pairs:
                continue
            dist = haversine_km(*COORDS[a][:2], *COORDS[b][:2])
            F_ref_a = built[a]["forcing"]
            F_ref_b = np.interp(built[b]["dates"], built[a]["dates"], built[a]["forcing"])
            for pi, pj in pairs:
                ma, _, _ = wa[pi]
                amp_a, phase_a, _ = fit_at_M(F_ref_a, built[a]["data"], ma)
                amp_b, phase_b, _ = fit_at_M(F_ref_b, built[b]["data"], ma)
                dphi = (phase_a - phase_b + math.pi) % (2 * math.pi) - math.pi
                harm, n = is_harmonic(ma)
                records.append((dist, ma, harm, math.cos(dphi), amp_a, amp_b))

    print(f"\n[winding_variogram] {len(records)} matched-winding instances across "
          f"{len(point_indices)*(len(point_indices)-1)//2} pairs")

    dist = np.array([r[0] for r in records])
    cosd = np.array([r[3] for r in records])
    harm_mask = np.array([r[2] for r in records])

    r_harm = np.corrcoef(dist[harm_mask], cosd[harm_mask])[0, 1] if harm_mask.sum() > 2 else np.nan
    r_nonharm = np.corrcoef(dist[~harm_mask], cosd[~harm_mask])[0, 1] if (~harm_mask).sum() > 2 else np.nan
    print(f"\n[winding_variogram] correlation(distance, cos_dphi):")
    print(f"  harmonic-family windings   (n={harm_mask.sum():3d}): r={r_harm:+.3f}   mean cos(dphi)={cosd[harm_mask].mean():+.3f}")
    print(f"  non-harmonic windings      (n={(~harm_mask).sum():3d}): r={r_nonharm:+.3f}   mean cos(dphi)={cosd[~harm_mask].mean():+.3f}")

    # binned empirical variogram
    bins = [0, 2000, 5000, 8000, 12000, 16000, 20000]
    print(f"\n[winding_variogram] binned empirical variogram (mean cos(dphi) per distance bin):")
    print(f"  {'bin (km)':>16s}  {'harmonic (n, mean)':>22s}  {'non-harmonic (n, mean)':>24s}")
    bin_centers, harm_means, nonharm_means = [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m_bin = (dist >= lo) & (dist < hi)
        h = m_bin & harm_mask
        nh = m_bin & ~harm_mask
        h_mean = cosd[h].mean() if h.sum() > 0 else np.nan
        nh_mean = cosd[nh].mean() if nh.sum() > 0 else np.nan
        bin_centers.append((lo + hi) / 2)
        harm_means.append(h_mean)
        nonharm_means.append(nh_mean)
        print(f"  {f'{lo}-{hi}':>16s}  {f'n={h.sum():3d}  {h_mean:+.3f}':>22s}  {f'n={nh.sum():3d}  {nh_mean:+.3f}':>24s}")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.scatter(dist[harm_mask], cosd[harm_mask], alpha=0.5, color="tab:blue", label="harmonic-family", s=25)
    ax.scatter(dist[~harm_mask], cosd[~harm_mask], alpha=0.5, color="tab:red", label="non-harmonic", s=25, marker="x")
    ax.plot(bin_centers, harm_means, color="tab:blue", lw=2.5, marker="o", label="harmonic bin mean")
    ax.plot(bin_centers, nonharm_means, color="tab:red", lw=2.5, marker="s", label="non-harmonic bin mean")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("great-circle distance (km)")
    ax.set_ylabel("cos(phase difference)  [shared-forcing corrected]")
    ax.set_title("Winding coherence vs. distance: empirical variogram,\n"
                 "harmonic-family vs. non-harmonic windings", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = ROOT / "winding_variogram.png"
    fig.savefig(out, dpi=140)
    print(f"\n[winding_variogram] saved {out}")


if __name__ == "__main__":
    main()
