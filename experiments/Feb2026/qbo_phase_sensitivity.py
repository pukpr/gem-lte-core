#!/usr/bin/env python3
"""'What if' sensitivity analysis: how much does the reconstructed QBO
reversal timing near 2015-2018 (the real, documented QBO disruption) move
if the two dominant lunar periods (27.2122d, 9.1207d/9.1085d) are nudged by
a small fractional amount?

Motivation: a monthly-sampled series is a genuinely marginal resolution
for a mechanism keyed to a ~29.4-day cycle -- barely coarser than the
driving period itself. If the model's regime-transition timing is highly
sensitive to tiny (much smaller than any real ephemeris uncertainty)
perturbations of the assumed lunar period, that's directly informative
about whether the 2015-2018 anomaly needs a special explanation, or falls
naturally out of ordinary phase wobble in a marginally-resolved system.

Important interpretive caveat, stated up front rather than left implicit:
extreme sensitivity to infinitesimal perturbations is NOT unambiguously
good news for the hypothesis. It cuts both ways -- it could mean a real
anomaly like 2015-2018 is unremarkable "ordinary" phase noise (supporting
the user's reading), or it could mean this specific numerical construction
is operating near an ill-conditioned / bifurcation-like regime, where
extreme sensitivity is itself a red flag about the model's reliability,
independent of whether real lunar forcing exists. Both readings are
consistent with the same sensitivity result; this script measures the
sensitivity, it doesn't resolve which reading is correct.

Reuses the full-42-term port from qbo_isolated_iir.py (corr=0.409 against
the real Forcing column -- an imperfect but structurally faithful stand-in,
not a certified replica; sensitivity *shape* is more trustworthy here than
any single absolute date).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns
from qbo_isolated_iir import tide_sum, amplify, iir, load_params

ROOT = Path(__file__).resolve().parent
TARGET_PERIODS = (27.2122, 9.1207, 9.1085)  # the terms being perturbed
DISRUPTION_WINDOW = (2015.0, 2018.5)  # the real, documented QBO anomaly


def perturbed_forcing(year, all_terms, delta, delA, delB, asym, integ, mA, mP,
                       init, mode="period"):
    """mode='period': fractional change to the target periods themselves
    (a claim about ephemeris/astronomical uncertainty -- essentially none
    exists in reality). mode='phase': additive radians shift to the target
    terms' phase only, period unchanged (a claim about modeling/epoch
    slop -- much more defensible as genuinely uncertain)."""
    terms = []
    for period_d, amp, phase in all_terms:
        if any(abs(period_d - t) < 0.01 for t in TARGET_PERIODS):
            if mode == "period":
                period_d = period_d * (1.0 + delta)
            else:
                phase = phase + delta
        terms.append((period_d, amp, phase))
    raw = tide_sum(year, terms, integ)
    imp = amplify(raw, year, delA, delB, asym)
    return iir(imp, mA, mP, init)


def reversal_dates(year, x):
    xc = x - x.mean()
    zc_idx = np.where(np.diff(np.sign(xc)) != 0)[0]
    return year[zc_idx]


def nearest_reversal(dates, target=2016.5):
    if len(dates) == 0:
        return np.nan
    return dates[np.argmin(np.abs(dates - target))]


def main():
    idx = "qbo30"
    year, dt, model, obs, forcing = load_columns(idx)
    p = load_params(idx)
    all_terms = [tuple(e) for e in p["lpap"]]
    delA, delB, asym = p["delA"], p["delB"], p["asym"]
    init, mA, mP = p["init"], p["ma"], p["mp"]
    integ = p.get("shfT", 0.0)

    baseline = perturbed_forcing(year, all_terms, 0.0, delA, delB, asym,
                                  integ, mA, mP, init)
    base_reversals = reversal_dates(year, baseline)
    base_nearest = nearest_reversal(base_reversals)
    print(f"[{idx}] baseline (delta=0) nearest reversal to 2016.5: "
          f"{base_nearest:.3f}")

    deltas = np.linspace(-0.01, 0.01, 41)  # +/- 1% fractional period change
    nearest_dates = []
    for d in deltas:
        f = perturbed_forcing(year, all_terms, d, delA, delB, asym, integ,
                               mA, mP, init, mode="period")
        rv = reversal_dates(year, f)
        nearest_dates.append(nearest_reversal(rv))
    nearest_dates = np.array(nearest_dates)
    shift_months = (nearest_dates - base_nearest) * 12.0

    phase_deltas = np.linspace(-0.5, 0.5, 41)  # +/- 0.5 rad (~4.6% of a cycle)
    phase_nearest_dates = []
    for d in phase_deltas:
        f = perturbed_forcing(year, all_terms, d, delA, delB, asym, integ,
                               mA, mP, init, mode="phase")
        rv = reversal_dates(year, f)
        phase_nearest_dates.append(nearest_reversal(rv))
    phase_nearest_dates = np.array(phase_nearest_dates)
    phase_shift_months = (phase_nearest_dates - base_nearest) * 12.0

    def report(name, xs, shifts, to_physical):
        valid = ~np.isnan(shifts)
        if valid.sum() <= 2:
            return
        i6 = np.argmin(np.abs(np.abs(shifts[valid]) - 6.0))
        x6 = xs[valid][i6]
        print(f"{name}: value needed for a ~6-month reversal-timing shift: "
              f"{x6:+.5g}  ({to_physical(x6)})")

    report("period perturbation", deltas, shift_months,
            lambda d: f"a {abs(d)*27.2122:.4f}-day shift in the 27.2122d period")
    report("phase perturbation", phase_deltas, phase_shift_months,
            lambda d: f"{abs(d)/(2*np.pi)*100:.2f}% of a full cycle, or "
                      f"{abs(d)/(2*np.pi)*27.2122:.3f} days of phase, at "
                      f"the 27.2122d rate")

    fig, axes = plt.subplots(4, 1, figsize=(11, 14))

    axes[0].plot(deltas * 100, shift_months, "o-", color="tab:blue")
    axes[0].axhline(0, color="black", linewidth=0.5)
    axes[0].axhline(6, color="tab:red", linestyle="--", linewidth=0.8,
                     label="+/-6mo (proxy for the real 2015-2018 anomaly size)")
    axes[0].axhline(-6, color="tab:red", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("fractional PERIOD perturbation on 27.2122d/9.12xxd "
                        "terms (%) -- essentially zero real uncertainty here")
    axes[0].set_ylabel("shift in nearest reversal date (months)")
    axes[0].set_title(f"{idx}: sensitivity to PERIOD perturbations "
                       f"(ephemeris-uncertainty framing)", fontweight="bold")
    axes[0].legend(loc="upper left", fontsize=8)

    axes[1].plot(phase_deltas / (2 * np.pi) * 100, phase_shift_months,
                 "o-", color="tab:purple")
    axes[1].axhline(0, color="black", linewidth=0.5)
    axes[1].axhline(6, color="tab:red", linestyle="--", linewidth=0.8,
                     label="+/-6mo (proxy for the real 2015-2018 anomaly size)")
    axes[1].axhline(-6, color="tab:red", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("PHASE perturbation on the same terms (% of a full "
                        "cycle) -- genuinely plausible modeling slop")
    axes[1].set_ylabel("shift in nearest reversal date (months)")
    axes[1].set_title("sensitivity to PHASE perturbations (epoch/modeling-"
                       "slop framing)", fontweight="bold")
    axes[1].legend(loc="upper left", fontsize=8)

    for d, c, lbl in [(0.0, "black", "baseline (delta=0)"),
                       (0.003, "tab:orange", "period delta=+0.3%"),
                       (-0.003, "tab:green", "period delta=-0.3%")]:
        f = perturbed_forcing(year, all_terms, d, delA, delB, asym, integ,
                               mA, mP, init, mode="period")
        mask = (year > 2008) & (year < 2022)
        axes[2].plot(year[mask], f[mask], color=c, alpha=0.8, label=lbl,
                     linewidth=1.0)
    axes[2].axvspan(*DISRUPTION_WINDOW, color="gray", alpha=0.15,
                     label="real QBO disruption window")
    axes[2].set_xlabel("year")
    axes[2].set_ylabel("reconstructed Forcing (full 42-term port)")
    axes[2].set_title("reconstructed signal, 2008-2022, for baseline vs. "
                       "+/-0.3% period perturbations", fontweight="bold")
    axes[2].legend(loc="upper right", fontsize=8)

    axes[3].plot(year, forcing, color="black", linewidth=0.6,
                 label="real Forcing column (lte_results.csv)")
    axes[3].axvspan(*DISRUPTION_WINDOW, color="gray", alpha=0.15,
                     label="real QBO disruption window")
    axes[3].set_xlabel("year")
    axes[3].set_ylabel("real Forcing")
    axes[3].set_title("for reference: the real Forcing column around the "
                       "same window", fontweight="bold")
    axes[3].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    out = ROOT / idx / "qbo_phase_sensitivity.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
