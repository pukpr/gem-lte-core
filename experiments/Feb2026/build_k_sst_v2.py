#!/usr/bin/env python3
"""build_k_sst_v2.py — re-runs the k_sst regional sweep using properly
GATED winding discovery (winding_rank.rank_series: AR1-floor + FWHM +
continuity triple test) instead of build_k_sst.py's crude, ungated
top-3-by-excess-R2 method. Each region gets however many windings
actually PASS the gate (0 to N), not a forced count of 3.

Reuses build_k_sst.py's region extraction/donor-manifold machinery
unchanged; only the discovery step and the lt.exe.p write-out differ.
Regions with zero passing ridges still get a lt.exe.p (backbone-only,
flagged n_passing=0 in the summary) so they stay in the comparison
table rather than silently disappearing.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ws  # noqa: E402
import build_k_sst as bk  # noqa: E402
import winding_rank as wr  # noqa: E402

OUT_JSON = bk.K_SST_DIR / "k_sst_summary_v2.json"


def main() -> int:
    donor = bk.load_donor_params()
    print(f"[donor] borrowing shared manifold shape from {bk.DONOR_INDEX}, "
          f"backbone_k1={donor['backbone_k1']:.6f}")
    regions = bk.extract_regions()
    print(f"[extract] {len(regions)} regions (reusing build_k_sst.py's "
          f"extraction -- same regions, .dat files already on disk)")

    yl = ws.year_length(0.0, float(donor["year"]))
    results = []
    t_start = time.time()
    for i, r in enumerate(regions):
        idx = r["name"]
        idx_dir = bk.K_SST_DIR / idx
        fake_prep = dict(params={**donor, "ltep": [donor["backbone_k1"]]},
                          yl=yl, dates=r["dates"], idate=r["dates"][0], nm=1,
                          year_startup=0.0, year_cand=float(donor["year"]))
        F = ws.build_forcing_at(fake_prep, r["dates"])

        # winding_rank.rank_series's default window (sigma=15yr, t0_step=5yr)
        # was calibrated against the full 1856-2023 (167yr) Kaplan record.
        # Direct AR1-noise calibration (20 trials) showed it is NOT
        # length-invariant: on a 73yr (post-1950) record it lets through a
        # mean of 3.25 spurious "significant" ridges/trial (vs 0.15/trial
        # at 167yr) because the shorter span leaves too few independent
        # sigma-windows to power the continuity>=0.70 test. Scaling sigma
        # and t0_step proportionally to the ACTUAL record span (same
        # window-count density as the 167yr baseline) restored calibration
        # to 0.3/trial in the same test -- close enough to the 0.15
        # baseline to trust. Applied generically here (not hardcoded to one
        # record length) so this stays correct if MIN_YEAR or the dataset
        # ever changes.
        span_years = float(r["dates"][-1] - r["dates"][0])
        sigma_scaled = span_years * (15.0 / 167.0)
        t0_step_scaled = span_years * (5.0 / 167.0)
        ridges = wr.rank_series(r["dates"], r["values"], F,
                                 sigma=sigma_scaled, t0_step=t0_step_scaled)
        # M~0 is a legitimate winding_rank ridge (a real, very-low-frequency
        # mode -- the same kind of thing as AMO's own 0.0134 winding) but
        # cos(2*pi*0*F)=1 for every t, an EXACT duplicate of the model's
        # own intercept column -- crashed the very first region carrying
        # one with a singular-matrix error. Excluded here, matching the
        # 0.002 < |M| convention ws.py's own local mode already uses.
        passing = [x for x in ridges if x["pass"] and abs(x["M"]) > 0.01]
        top_M = [x["M"] for x in passing]

        # de-duplicate against each other AND against the appended backbone
        # (near-identical windings make the design matrix ill-conditioned
        # even short of exact singularity) -- keep the higher peak_bits one
        final_M = []
        for M in sorted(top_M, key=lambda m: -next(
                x["peak_bits"] for x in passing if x["M"] == m)):
            if all(abs(M - f) > 0.03 for f in final_M) and abs(M - donor["backbone_k1"]) > 0.03:
                final_M.append(M)
            elif abs(M - donor["backbone_k1"]) <= 0.03:
                pass  # already effectively the shared backbone, skip the duplicate
        top_M = final_M

        params = {k: v for k, v in donor.items() if k != "backbone_k1"}
        params["ltep"] = top_M + [donor["backbone_k1"]]
        params["IR"] = 0.0
        idx_dir.mkdir(parents=True, exist_ok=True)
        (idx_dir / "lt.exe.p").write_text(json.dumps(params, indent=1))

        elapsed = time.time() - t_start
        eta = elapsed / (i + 1) * (len(regions) - i - 1)
        try:
            prep = ws.load_index(bk.K_SST_DIR, idx)
            g = ws.fit_global(prep)
            train_r, val_r = g["train_r"], g["val_r"]
        except Exception as exc:  # noqa: BLE001 -- one bad region must not kill an ~28min run
            print(f"  [{i+1:3d}/{len(regions)}] {idx:12s} SKIPPED (fit error: {exc})")
            train_r = val_r = float("nan")
        results.append(dict(name=idx, lat=r["lat_c"], lon=r["lon_c"],
                             n_passing=len(passing), top_M=top_M,
                             all_ridges=ridges,
                             train_r=train_r, val_r=val_r))
        print(f"  [{i+1:3d}/{len(regions)}] {idx:12s} lat={r['lat_c']:+6.1f} "
              f"lon={r['lon_c']:+6.1f}  n_pass={len(passing)}  "
              f"M={np.round(top_M,3).tolist()}  train_r={train_r:+.3f}  "
              f"val_r={val_r:+.3f}   (elapsed={elapsed/60:.1f}min, "
              f"ETA={eta/60:.1f}min)")

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\n[written] {OUT_JSON}")

    n_pass = np.array([r["n_passing"] for r in results])
    val_r = np.array([r["val_r"] for r in results])
    print(f"\n[summary] regions with 0 passing ridges: {np.sum(n_pass==0)}/{len(results)}")
    print(f"[summary] val_r: mean={val_r.mean():+.3f} median={np.median(val_r):+.3f} "
          f"min={val_r.min():+.3f} max={val_r.max():+.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
