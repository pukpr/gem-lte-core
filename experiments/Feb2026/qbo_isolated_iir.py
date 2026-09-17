#!/usr/bin/env python3
"""Push the isolated 9.1207d+9.1085d Doodson pair through the ACTUAL
Impulse_Amplify + IIR(mP) pipeline, instead of comparing the raw cosine
directly (qbo_compensating_dtw.py's simplified version).

This never touches or zeroes the other 40 lpap entries -- Tide_Sum is a
linear sum over constituents, so isolating two terms just means summing
only those two; nothing needs to be edited in lt.exe.p. The other 40 terms
are only read (not modified) for one purpose: validating this Python port
of Amplify+IIR against the REAL Forcing column in lte_results.csv, by
reconstructing Forcing from the FULL 42-term sum and checking it actually
matches. If that validation fails, the isolated-term result below isn't
trustworthy either -- so it has to be checked first.

Ported directly from the Ada source, matching each named parameter:
  - Amplify (gem-lte-primitives.adb) with Offset=Ramp=0 (the values
    Calc_Forcing actually passes): Impulses(t) = TideSum(t) * Impulse(t)
  - Impulse_Delta (gem-lte-primitives-solution.adb): a twice-yearly Dirac
    comb -- zero every month except one (value delA) and its opposite
    6-month mark (value asym). IMPULSE (Impulse_Only) is unset in
    qbo30/lt.exe.resp -> default True -> Impulse_Delta, not the smeared
    variant.
  - IIR (gem-lte-primitives.adb): Res(I) = Impulses(I) + Mem*Res(I-1) -
    sign(Res(I-1))*mP, Mem = 1 - mA, Res(start) = init.

Important limitation to keep in view: because IIR is nonlinear (the drag
term's sign depends on the running state), running the isolated 2-term
input through this SAME recursion is NOT the same as extracting "this
pair's contribution" out of the real 42-term Forcing -- IIR(sum of terms)
!= sum of IIR(each term). What this computes is a well-posed, different
question: "what would Forcing look like if these two terms were the ONLY
tidal input" -- informative on its own, not a decomposition of the real
signal.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns
from qbo_phase_lock import locked_M, phase_residual
from qbo_compensating_dtw import (standardize, dtw_distance, null_baseline,
                                   lpap_term)

ROOT = Path(__file__).resolve().parent
YEAR_LEN = 365.242484


def tide_sum(year, terms, integ=0.0):
    """terms: list of (period_days, amp, phase). Linear -- isolating a
    subset just means summing fewer terms, nothing else changes. Matches
    Tide_Sum_Diff (gem-lte-primitives.adb) with Scaling=0 (no TF1*TF2
    cross-term, as Calc_Forcing actually calls it): each constituent
    contributes Amplitude*Cos(2*pi*Freq*t+Phase) - Integ*Freq*Sin(...),
    Integ = D.B.ShiftT ("shfT" in lt.exe.p)."""
    s = np.zeros_like(year)
    for period_d, amp, phase in terms:
        freq = YEAR_LEN / period_d
        angle = 2.0 * np.pi * freq * year + phase
        s += amp * np.cos(angle) - integ * freq * np.sin(angle)
    return s


def impulse_delta(year, delA, delB, asym, sampling=12.0):
    frac = year - np.floor(year)
    trunc = np.floor(frac * sampling).astype(int)
    dpos = int(abs(delB) * sampling)
    dpos_opp = (dpos + int(sampling) // 2) % 12
    out = np.zeros_like(year)
    out[trunc == dpos] = delA
    out[trunc == dpos_opp] = asym
    return out


def amplify(raw, year, delA, delB, asym, sampling=12.0):
    return raw * impulse_delta(year, delA, delB, asym, sampling)


def iir(raw, mA, mP, init):
    mem = min(max(1.0 - mA, 0.0), 1.0)
    res = np.empty_like(raw)
    res[0] = init
    for i in range(1, len(raw)):
        ramp = mP if res[i - 1] >= 0 else -mP
        res[i] = raw[i] + mem * res[i - 1] - ramp
    return res


def load_params(idx: str):
    p = json.loads((ROOT / idx / "lt.exe.p").read_text())
    return p


def main():
    idx = "qbo30"
    year, dt, model, obs, forcing = load_columns(idx)
    p = load_params(idx)
    delA, delB, asym = p["delA"], p["delB"], p["asym"]
    init, mA, mP = p["init"], p["ma"], p["mp"]
    integ = p.get("shfT", 0.0)

    # --- Step 1: validate the port against the REAL Forcing column ---
    all_terms = [tuple(e) for e in p["lpap"]]
    raw_full = tide_sum(year, all_terms, integ)
    imp_full = amplify(raw_full, year, delA, delB, asym)
    forcing_full = iir(imp_full, mA, mP, init)

    corr = np.corrcoef(forcing_full, forcing)[0, 1]
    rmse = np.sqrt(np.mean((forcing_full - forcing) ** 2))
    print(f"[{idx}] port validation vs real Forcing column: "
          f"corr={corr:.5f}  rmse={rmse:.4f}  "
          f"(real Forcing std={forcing.std():.4f})")

    # --- Step 2: isolated 2-term reconstruction, same pipeline ---
    t9a = lpap_term(idx, 9.1207)
    t9b = lpap_term(idx, 9.1085)
    raw_iso = tide_sum(year, [t9a, t9b], integ)
    imp_iso = amplify(raw_iso, year, delA, delB, asym)
    forcing_iso = iir(imp_iso, mA, mP, init)

    M, decay, mP_locked = locked_M(idx)
    residual, slope = phase_residual(year, forcing, M)
    residual_s = standardize(residual)
    iso_s = standardize(forcing_iso)

    d, path = dtw_distance(residual_s, iso_s)
    rng = np.random.default_rng(0)
    null_mean, null_std = null_baseline(rng, residual_s, forcing_iso)
    z = (null_mean - d) / null_std if null_std > 0 else float("nan")
    print(f"[{idx}] isolated-pair-through-real-IIR vs phase residual: "
          f"DTW dist={d:.4f}  null mean={null_mean:.4f}  "
          f"null std={null_std:.4f}  z={z:.2f}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    axes[0].plot(year, forcing, color="black", linewidth=0.6,
                 label="real Forcing (lte_results.csv)")
    axes[0].plot(year, forcing_full, color="tab:green", linewidth=0.6,
                 alpha=0.7, label=f"ported full-42-term reconstruction "
                                   f"(corr={corr:.4f})")
    axes[0].set_title(f"{idx}: port validation -- reconstructed Forcing "
                       f"(all 42 terms) vs. the real column",
                       fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].plot(year, residual_s, color="tab:blue",
                 label="phase residual (standardized)")
    axes[1].plot(year, iso_s, color="tab:red", alpha=0.8,
                 label="isolated 9.1207d+9.1085d through real Amplify+IIR")
    axes[1].set_title(f"isolated pair, run through the actual nonlinear "
                       f"pipeline (DTW z={z:.2f} vs. null)",
                       fontweight="bold")
    axes[1].legend(loc="upper right", fontsize=8)

    path_arr = np.array(path)
    axes[2].plot(year[path_arr[:, 0]], year[path_arr[:, 1]], color="black",
                 linewidth=0.7)
    axes[2].plot([year.min(), year.max()], [year.min(), year.max()],
                 color="gray", linestyle="--", linewidth=0.7,
                 label="no-warp diagonal")
    axes[2].set_xlabel("residual time (yr)")
    axes[2].set_ylabel("candidate time (yr)")
    axes[2].set_title("DTW warping path", fontweight="bold")
    axes[2].legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    out = ROOT / idx / "qbo_isolated_iir.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
