#!/usr/bin/env python3
"""Test the user's subharmonic-quadrature-metastability hypothesis
(2026-09-29): Ada's per-quad Regression_Factors fits each quad's own
manifold independently via OLS on [sin(k*F), cos(k*F)] basis columns for
each winding k -- for the SLOWEST winding (the ~BACKBONE/11 subharmonic
responsible for the 60-120yr AMO-type modulation), the period is
comparable to the record length, so sin and cos of that argument are
nearly collinear over the observed record: which quadrature (cos-
dominant vs sin-dominant) the regression lands on is effectively
arbitrary/metastable, decided independently per quad with no
communication between quads. Direct amo Feb2026-vs-Sep2026 comparison
showed a real but partial (r~0.26) sign-flip relationship between two
independent fits of the identical data, not a clean quadrature identity
-- consistent with "related but not the same local optimum".

This script tests the COLLECTIVE consequence across all 89 quads: if
the quadrature choice is genuinely arbitrary per-quad (not physically
determined by geography), the fitted (amp*cos(phase), amp*sin(phase))
vector for the subharmonic winding should look spatially INCOHERENT --
poorly explained by a smooth function of (lat, lon) -- compared to
other, faster windings whose period is short relative to the record
(and so are NOT quadrature-metastable, since sin/cos are well-separated
there). A canonicalization (force amp>=0, fold phase into a fixed half-
plane) is also tested: if coherence improves substantially after
canonicalizing, the residual "incoherence" was mostly this sign/
quadrature bookkeeping artifact, not real geographic disorder --
relevant to any per-quad "long-term trend" that shares a time axis with
this slow term, since a near-degenerate sin/cos split there can leak
into the fitted trend/level terms too.
"""
import json
import sys

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep
import plot_winding_regions as W
import mlr_shared_windings as M

import numpy as np

TARGET_K = M.SUBHARMONIC  # BACKBONE/11 ~= 0.0188
K_TOLERANCE = 0.010        # accept the nearest winding within +-this of TARGET_K
CONTROL_KS = [M.BACKBONE, 2 * M.BACKBONE, 4 * M.BACKBONE]  # fast windings, control group


def load_field(target_k: float, tol: float):
    """For every quad, find the k_amp_phase triple nearest target_k
    (within tol); return (lat, lon, A, B) where A=amp*cos(phase),
    B=amp*sin(phase) -- the raw Cartesian (non-canonical) components."""
    cells = sorted((p.name for p in sweep.FEB.iterdir()
                     if p.is_dir() and W.NAME_RE.match(p.name)))
    rows = []
    for cell in cells:
        path = W.freshest_windings(cell)
        if path is None:
            continue
        d = json.loads(path.read_text())
        triples = d.get("k_amp_phase", [])
        if not triples:
            continue
        best = min(triples, key=lambda t: abs(abs(t[0]) - target_k))
        if abs(abs(best[0]) - target_k) > tol:
            continue
        k, amp, phase = best
        lat, lon = W.decode(cell)
        rows.append(dict(cell=cell, lat=lat, lon=lon, k=k, amp=amp, phase=phase,
                         A=amp * np.cos(phase), B=amp * np.sin(phase)))
    return rows


def spatial_r2(rows, canonicalize: bool) -> tuple[float, float, int]:
    """Fit A(lat,lon) and B(lat,lon) each to the shared degree-1 spatial
    basis (same basis as mlr_shared_windings.py) and report R^2 for
    each -- the measure of "how spatially smooth is this field". With
    canonicalize=True, first force amp>=0 and fold phase into
    [-pi/2, pi/2) (a fixed half-plane) by flipping sign where needed --
    removes the arbitrary +-(amp,phase+pi) <-> (-amp,phase) degeneracy
    that free regression can land on for BOTH quadrature-metastable and
    ordinary windings alike, isolating whether the SPECIFIC quadrature
    ambiguity (not just an amp-sign flip) is the dominant effect."""
    lat = np.array([r["lat"] for r in rows])
    lon = np.array([r["lon"] for r in rows])
    A = np.array([r["A"] for r in rows])
    B = np.array([r["B"] for r in rows])

    if canonicalize:
        # amp*cos(phase) < 0  =>  flip sign of (A,B) together (equivalent
        # to phase += pi, amp unchanged in magnitude) so every row's A
        # component is >= 0 -- a fixed, arbitrary-free convention.
        flip = A < 0
        A = np.where(flip, -A, A)
        B = np.where(flip, -B, B)

    X = np.array([M.spatial_basis(la, lo) for la, lo in zip(lat, lon)])
    coeffs_A, *_ = np.linalg.lstsq(X, A, rcond=None)
    coeffs_B, *_ = np.linalg.lstsq(X, B, rcond=None)
    pred_A = X @ coeffs_A
    pred_B = X @ coeffs_B
    ss_res = np.sum((A - pred_A) ** 2) + np.sum((B - pred_B) ** 2)
    ss_tot = np.sum((A - A.mean()) ** 2) + np.sum((B - B.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    return r2, ss_res, len(rows)


def main() -> None:
    print(f"Subharmonic target k = {TARGET_K:.5f} (+-{K_TOLERANCE})\n")

    rows = load_field(TARGET_K, K_TOLERANCE)
    print(f"Subharmonic: {len(rows)} quads matched")
    r2_raw, _, n = spatial_r2(rows, canonicalize=False)
    r2_canon, _, _ = spatial_r2(rows, canonicalize=True)
    print(f"  spatial R^2 (raw, as-fit)        : {r2_raw:.4f}")
    print(f"  spatial R^2 (sign-canonicalized) : {r2_canon:.4f}")
    print(f"  delta: {r2_canon - r2_raw:+.4f}\n")

    print("Control group (fast windings, period << record length, should NOT")
    print("be quadrature-metastable -- same test applied for comparison):")
    for ck in CONTROL_KS:
        crows = load_field(ck, tol=0.03)
        if len(crows) < 10:
            print(f"  k={ck:.4f}: only {len(crows)} matched quads, skipping")
            continue
        cr2_raw, _, cn = spatial_r2(crows, canonicalize=False)
        cr2_canon, _, _ = spatial_r2(crows, canonicalize=True)
        print(f"  k={ck:.4f} ({cn} quads): raw R^2={cr2_raw:.4f}  "
              f"canon R^2={cr2_canon:.4f}  delta={cr2_canon-cr2_raw:+.4f}")

    print("\nInterpretation: if the subharmonic's raw->canonical R^2 jump is")
    print("much larger than the control windings' own jump, that's direct")
    print("evidence the quadrature choice (not just geography) is what's")
    print("breaking spatial coherence for this specific term.")


if __name__ == "__main__":
    main()
