#!/usr/bin/env python3
"""Test whether phase-resolved LTE features improve daily Singapore QVBO skill.

The baseline is the fixed monthly LTE column-3 model, interpolated to daily
times and calibrated to Singapore wind. The phase features represent only
intramonth deviations, and are exactly zero at LTE month-start anchors.

Example:
    python3 evaluate_singapore_qvbo_multiscale.py
"""

import argparse
import csv
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fit_singapore_qvbo_by_lte_hidden import (
    LTE_COLUMN3_COEFFICIENTS,
    LTE_COLUMN3_FREQUENCY,
    dates_to_lte_time,
    load_qvbo,
    predict,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lte-input", type=Path, default=Path("lte_results.csv"))
    parser.add_argument(
        "--qvbo-input",
        type=Path,
        default=Path("singapore_qvbo_30hpa_daily.csv"),
    )
    parser.add_argument(
        "--holdout-start",
        type=date.fromisoformat,
        help="First untouched test date, YYYY-MM-DD.",
    )
    parser.add_argument(
        "--holdout-end",
        type=date.fromisoformat,
        help="Last untouched test date, YYYY-MM-DD.",
    )
    parser.add_argument(
        "--holdout-years",
        type=int,
        default=3,
        help="Latest calendar years used when holdout dates are omitted (default: 3).",
    )
    parser.add_argument(
        "--validation-years",
        type=int,
        default=2,
        help="Calendar years before holdout used to select ridge strength (default: 2).",
    )
    parser.add_argument(
        "--ridge-strengths",
        default="0,0.1,1,10,100,1000,10000",
        help="Comma-separated ridge strengths for phase terms.",
    )
    parser.add_argument(
        "--predictions-output",
        type=Path,
        default=Path("singapore_qvbo_multiscale_predictions.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("singapore_qvbo_multiscale_evaluation.png"),
    )
    return parser.parse_args()


def metrics(observed, predicted):
    return np.corrcoef(observed, predicted)[0, 1], np.sqrt(
        np.mean((observed - predicted) ** 2)
    )


def fit_linear(features, observed, mask, ridge_strength=0.0):
    """Fit intercept and monthly baseline freely; ridge-penalize phase terms."""
    phase_mean = features[mask, 2:].mean(axis=0)
    phase_scale = features[mask, 2:].std(axis=0)
    phase_scale[phase_scale == 0] = 1.0
    design = features.copy()
    design[:, 2:] = (design[:, 2:] - phase_mean) / phase_scale
    penalty = np.diag((0.0, 0.0, ridge_strength, ridge_strength))
    coefficients = np.linalg.solve(
        design[mask].T @ design[mask] + penalty,
        design[mask].T @ observed[mask],
    )
    return coefficients, phase_mean, phase_scale


def apply_linear(features, coefficients, phase_mean, phase_scale):
    design = features.copy()
    design[:, 2:] = (design[:, 2:] - phase_mean) / phase_scale
    return design @ coefficients


def parse_strengths(value):
    strengths = [float(item) for item in value.split(",")]
    if not strengths or any(item < 0 for item in strengths):
        raise ValueError("--ridge-strengths must contain non-negative values.")
    return strengths


def main():
    args = parse_args()
    if args.holdout_years < 1 or args.validation_years < 1:
        raise ValueError("--holdout-years and --validation-years must be positive.")
    if (args.holdout_start is None) != (args.holdout_end is None):
        raise ValueError("--holdout-start and --holdout-end must be supplied together.")
    if args.holdout_start and args.holdout_start > args.holdout_end:
        raise ValueError("--holdout-start must not be later than --holdout-end.")

    lte = np.loadtxt(args.lte_input, delimiter=",", ndmin=2)
    if lte.shape[1] < 4:
        raise ValueError(f"{args.lte_input} must contain at least four columns.")
    lte = lte[np.argsort(lte[:, 0])]
    lte_time, lte_hidden = lte[:, 0], lte[:, 3]
    dates, observed = load_qvbo(args.qvbo_input)
    daily_time = dates_to_lte_time(dates, lte_time)

    if args.holdout_start is None:
        last_date = dates[-1]
        holdout_start = date(last_date.year - args.holdout_years + 1, 1, 1)
        holdout_end = last_date
    else:
        holdout_start, holdout_end = args.holdout_start, args.holdout_end
    validation_start = date(holdout_start.year - args.validation_years, 1, 1)

    holdout = np.array([holdout_start <= value <= holdout_end for value in dates])
    validation = np.array(
        [validation_start <= value < holdout_start for value in dates]
    )
    training = np.array([value < validation_start for value in dates])
    if min(training.sum(), validation.sum(), holdout.sum()) < 30:
        raise ValueError("Training, validation, and holdout periods each need 30 rows.")

    monthly_model = predict(
        lte_hidden, LTE_COLUMN3_COEFFICIENTS, LTE_COLUMN3_FREQUENCY
    )
    baseline = np.interp(daily_time, lte_time, monthly_model)
    phase = LTE_COLUMN3_FREQUENCY * np.interp(daily_time, lte_time, lte_hidden)
    monthly_sine = np.sin(LTE_COLUMN3_FREQUENCY * lte_hidden)
    monthly_cosine = np.cos(LTE_COLUMN3_FREQUENCY * lte_hidden)
    # These phase deviations vanish at each stored monthly LTE timestamp.
    phase_sine = np.sin(phase) - np.interp(daily_time, lte_time, monthly_sine)
    phase_cosine = np.cos(phase) - np.interp(daily_time, lte_time, monthly_cosine)
    features = np.column_stack((np.ones_like(baseline), baseline, phase_sine, phase_cosine))

    selection_baseline_coefficients = np.linalg.lstsq(
        features[training, :2], observed[training], rcond=None
    )[0]
    selection_baseline_prediction = features[:, :2] @ selection_baseline_coefficients

    strengths = parse_strengths(args.ridge_strengths)
    validation_scores = {}
    for strength in strengths:
        coefficients, phase_mean, phase_scale = fit_linear(
            features, observed, training, strength
        )
        validation_prediction = apply_linear(
            features, coefficients, phase_mean, phase_scale
        )
        validation_scores[strength] = metrics(
            observed[validation], validation_prediction[validation]
        )[1]
    selected_strength = min(validation_scores, key=validation_scores.get)
    selection_coefficients, selection_phase_mean, selection_phase_scale = fit_linear(
        features, observed, training, selected_strength
    )
    selection_multiscale_prediction = apply_linear(
        features,
        selection_coefficients,
        selection_phase_mean,
        selection_phase_scale,
    )

    pre_holdout = training | validation
    coefficients, phase_mean, phase_scale = fit_linear(
        features, observed, pre_holdout, selected_strength
    )
    multiscale_prediction = apply_linear(
        features, coefficients, phase_mean, phase_scale
    )
    baseline_coefficients = np.linalg.lstsq(
        features[pre_holdout, :2], observed[pre_holdout], rcond=None
    )[0]
    baseline_prediction = features[:, :2] @ baseline_coefficients

    report_sets = (
        ("Training", training, baseline_prediction, multiscale_prediction),
        (
            "Validation",
            validation,
            selection_baseline_prediction,
            selection_multiscale_prediction,
        ),
        ("Holdout", holdout, baseline_prediction, multiscale_prediction),
    )
    for label, mask, baseline_values, multiscale_values in report_sets:
        baseline_correlation, baseline_rmse = metrics(
            observed[mask], baseline_values[mask]
        )
        multiscale_correlation, multiscale_rmse = metrics(
            observed[mask], multiscale_values[mask]
        )
        print(
            f"{label}: baseline r={baseline_correlation:.6f}, RMSE={baseline_rmse:.6f}; "
            f"multiscale r={multiscale_correlation:.6f}, RMSE={multiscale_rmse:.6f}"
        )
    print(f"Selected phase ridge strength = {selected_strength:g}")
    print(
        "Phase coefficients (standardized sine, cosine) = "
        f"{coefficients[2]:.6f}, {coefficients[3]:.6f}"
    )

    with args.predictions_output.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination)
        writer.writerow(
            (
                "date",
                "observed_u_30hpa_ms",
                "monthly_baseline",
                "phase_sine",
                "phase_cosine",
                "baseline_prediction",
                "multiscale_prediction",
                "split",
            )
        )
        for index, value in enumerate(dates):
            split = "holdout" if holdout[index] else "validation" if validation[index] else "training"
            writer.writerow(
                (
                    value.isoformat(),
                    f"{observed[index]:.6f}",
                    f"{baseline[index]:.12f}",
                    f"{phase_sine[index]:.12f}",
                    f"{phase_cosine[index]:.12f}",
                    f"{baseline_prediction[index]:.6f}",
                    f"{multiscale_prediction[index]:.6f}",
                    split,
                )
            )

    plot_dates = np.array(dates, dtype="datetime64[D]")
    figure, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="#e6f2ff")
    axes[0].plot(plot_dates, observed, color="steelblue", alpha=0.5, linewidth=0.6, label="Singapore daily")
    axes[0].plot(plot_dates, baseline_prediction, color="black", linewidth=1.1, label="Monthly baseline")
    axes[0].plot(plot_dates, multiscale_prediction, color="crimson", linewidth=1.0, label="Multiscale")
    axes[0].axvspan(
        np.datetime64(holdout_start),
        np.datetime64(holdout_end),
        color="gold",
        alpha=0.15,
        label="Untouched holdout",
    )
    axes[0].set_title("Singapore QVBO: monthly baseline and phase correction")
    axes[0].set_xlabel("Date")
    axes[0].set_ylabel("Zonal wind at 30 hPa (m/s)")
    axes[0].grid(alpha=0.2)
    axes[0].legend(loc="lower right")

    axes[1].scatter(
        observed[holdout], baseline_prediction[holdout],
        s=8, alpha=0.4, color="black", label="Baseline"
    )
    axes[1].scatter(
        observed[holdout], multiscale_prediction[holdout],
        s=8, alpha=0.4, color="crimson", label="Multiscale"
    )
    lower = min(observed[holdout].min(), baseline_prediction[holdout].min(), multiscale_prediction[holdout].min())
    upper = max(observed[holdout].max(), baseline_prediction[holdout].max(), multiscale_prediction[holdout].max())
    axes[1].plot((lower, upper), (lower, upper), "k--", linewidth=1)
    axes[1].set_title("Untouched holdout")
    axes[1].set_xlabel("Observed Singapore wind (m/s)")
    axes[1].set_ylabel("Predicted wind (m/s)")
    axes[1].grid(alpha=0.2)
    axes[1].legend(loc="upper left")
    figure.tight_layout()
    figure.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"Wrote {args.predictions_output}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
