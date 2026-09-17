#!/usr/bin/env python3
"""Fit lte_results.csv column 3 from column 4 with a sinusoidal trend.

Example:
    python3 fit_column3_by_column4.py

By default, D is derived from the column-4 monthly plateau slope:
D = 2*pi / abs(slope). Use --frequency-mode search to search D instead.

To assess the fitted relationship out of sample, withhold a time window:
    python3 fit_column3_by_column4.py --holdout-start-year 1980 --holdout-stop-year 1990
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


def fit_at_frequency(x: np.ndarray, y: np.ndarray, frequency: float):
    """Return mean squared error, linear coefficients, and fitted values."""
    design = np.column_stack(
        (np.ones_like(x), x, np.sin(frequency * x), np.cos(frequency * x))
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coefficients
    return np.mean((y - fitted) ** 2), coefficients, fitted


def predict(x: np.ndarray, coefficients: np.ndarray, frequency: float):
    """Apply fitted coefficients at a fixed frequency."""
    design = np.column_stack(
        (np.ones_like(x), x, np.sin(frequency * x), np.cos(frequency * x))
    )
    return design @ coefficients


def metrics(observed: np.ndarray, fitted: np.ndarray):
    """Return correlation coefficient and root mean squared error."""
    return np.corrcoef(observed, fitted)[0, 1], np.sqrt(
        np.mean((observed - fitted) ** 2)
    )


def plateau_frequency(time: np.ndarray, x: np.ndarray):
    """Derive D from the median low-magnitude monthly column-4 increment."""
    order = np.argsort(time)
    time = time[order]
    x = x[order]
    monthly_spacing = 1.0 / 12.0
    monthly_pairs = np.isclose(
        np.diff(time), monthly_spacing, rtol=0.0, atol=2e-6
    )
    increments = np.abs(np.diff(x)[monthly_pairs])
    if increments.size == 0:
        raise ValueError("No adjacent monthly column-4 values are available.")

    plateau_limit = np.quantile(increments, 0.25)
    plateau_increments = increments[increments <= plateau_limit]
    slope = np.median(plateau_increments)
    if slope <= 0:
        raise ValueError("The derived column-4 plateau slope must be non-zero.")
    return slope, 2.0 * np.pi / slope


def find_best_fit(
    x: np.ndarray, y: np.ndarray, maximum_frequency: float, grid_size: int
):
    """Profile the linear terms over frequency, then refine every grid minimum."""
    frequencies = np.linspace(1e-6, maximum_frequency, grid_size)
    errors = np.array([fit_at_frequency(x, y, frequency)[0] for frequency in frequencies])
    candidates = [0, len(frequencies) - 1]
    candidates.extend(
        index
        for index in range(1, len(frequencies) - 1)
        if errors[index] <= errors[index - 1] and errors[index] <= errors[index + 1]
    )

    best = None
    for index in candidates:
        if index in (0, len(frequencies) - 1):
            frequency = frequencies[index]
        else:
            result = minimize_scalar(
                lambda value: fit_at_frequency(x, y, value)[0],
                bounds=(frequencies[index - 1], frequencies[index + 1]),
                method="bounded",
                options={"xatol": 1e-12},
            )
            frequency = result.x

        mse, coefficients, fitted = fit_at_frequency(x, y, frequency)
        if best is None or mse < best[0]:
            best = mse, frequency, coefficients, fitted

    return best


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit column 3 as A + B*column4 + C*sin(D*column4 + E)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("lte_results.csv"),
        help="Headerless LTE CSV input file (default: lte_results.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("qbo30_column4_sinusoidal_fit.png"),
        help="Output plot path (default: qbo30_column4_sinusoidal_fit.png).",
    )
    parser.add_argument(
        "--max-frequency",
        type=float,
        default=50.0,
        help="Largest D for --frequency-mode search; use with --grid-size (default: 50).",
    )
    parser.add_argument(
        "--frequency-mode",
        choices=("plateau", "search"),
        default="plateau",
        help="Derive D from monthly plateau steps or search D (default: plateau).",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=5001,
        help="D samples for --frequency-mode search (default: 5001).",
    )
    parser.add_argument(
        "--start-year",
        type=float,
        default=1980.0,
        help="Start of the highlighted reference window (default: 1980).",
    )
    parser.add_argument(
        "--stop-year",
        type=float,
        default=1990.0,
        help="End of the highlighted reference window (default: 1990).",
    )
    parser.add_argument(
        "--holdout-start-year",
        type=float,
        help="First year withheld from fitting for out-of-sample validation.",
    )
    parser.add_argument(
        "--holdout-stop-year",
        type=float,
        help="Last year withheld from fitting for out-of-sample validation.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.max_frequency <= 0:
        raise ValueError("--max-frequency must be greater than zero.")
    if args.grid_size < 3:
        raise ValueError("--grid-size must be at least 3.")
    if args.stop_year < args.start_year:
        raise ValueError("--stop-year must not be earlier than --start-year.")
    has_holdout = args.holdout_start_year is not None
    if has_holdout != (args.holdout_stop_year is not None):
        raise ValueError(
            "--holdout-start-year and --holdout-stop-year must be supplied together."
        )
    if has_holdout and args.holdout_stop_year < args.holdout_start_year:
        raise ValueError("--holdout-stop-year must not be earlier than --holdout-start-year.")

    values = np.loadtxt(args.input, delimiter=",", ndmin=2)
    if values.shape[1] < 4:
        raise ValueError(f"{args.input} must contain at least four columns.")

    time = values[:, 0]
    x = values[:, 3]  # Column 4, using 1-based column numbering.
    y = values[:, 2]  # Column 3, using 1-based column numbering.
    holdout_mask = np.zeros(time.shape, dtype=bool)
    if has_holdout:
        holdout_mask = (time >= args.holdout_start_year) & (
            time <= args.holdout_stop_year
        )
        if holdout_mask.sum() < 2:
            raise ValueError("The holdout window must contain at least two rows.")
        if (~holdout_mask).sum() < 4:
            raise ValueError("At least four rows must remain available for fitting.")

    fit_mask = ~holdout_mask
    if args.frequency_mode == "plateau":
        slope, frequency = plateau_frequency(time[fit_mask], x[fit_mask])
        _, coefficients, fitted_training = fit_at_frequency(
            x[fit_mask], y[fit_mask], frequency
        )
    else:
        slope = None
        _, frequency, coefficients, fitted_training = find_best_fit(
            x[fit_mask], y[fit_mask], args.max_frequency, args.grid_size
        )
    fitted = predict(x, coefficients, frequency)

    a, b, sine_coefficient, cosine_coefficient = coefficients
    amplitude = np.hypot(sine_coefficient, cosine_coefficient)
    phase = np.arctan2(cosine_coefficient, sine_coefficient)
    training_correlation, training_rmse = metrics(y[fit_mask], fitted_training)

    if slope is None:
        print(f"Best fit over 0 <= D <= {args.max_frequency:g}")
    else:
        print(f"Plateau monthly slope = {slope:.12f}")
        print("D = 2*pi / abs(plateau monthly slope)")
    print(f"A = {a:.12f}")
    print(f"B = {b:.12f}")
    print(f"C = {amplitude:.12f}")
    print(f"D = {frequency:.12f}")
    print(f"E = {phase:.12f} radians")
    print(f"Training rows = {fit_mask.sum()}")
    print(f"Training correlation coefficient = {training_correlation:.12f}")
    print(f"Training RMSE = {training_rmse:.12f}")
    if has_holdout:
        holdout_correlation, holdout_rmse = metrics(
            y[holdout_mask], fitted[holdout_mask]
        )
        print(f"Holdout rows = {holdout_mask.sum()}")
        print(f"Holdout correlation coefficient = {holdout_correlation:.12f}")
        print(f"Holdout RMSE = {holdout_rmse:.12f}")
        window_start, window_stop = args.holdout_start_year, args.holdout_stop_year
        plot_title = (
            "Column 3 fit from column 4: "
            f"holdout r = {holdout_correlation:.6f}, RMSE = {holdout_rmse:.6f}"
        )
        window_label = f"{window_start:g}-{window_stop:g} holdout window"
    else:
        window_start, window_stop = args.start_year, args.stop_year
        plot_title = (
            "Column 3 fit from column 4: "
            f"r = {training_correlation:.6f}, RMSE = {training_rmse:.6f}"
        )
        window_label = f"{window_start:g}-{window_stop:g} reference window"

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(16, 6),
        gridspec_kw={"width_ratios": [2.5, 1]},
        facecolor="#e6f2ff",
    )
    axes[0].plot(time, y, color="blue", linewidth=1, alpha=0.75, label="Column 3")
    axes[0].plot(time, fitted, color="red", linewidth=1.25, label="Sinusoidal fit")
    axes[0].axvspan(
        window_start,
        window_stop,
        color="black",
        alpha=0.08,
        label=window_label,
    )
    axes[0].set_title(plot_title)
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Column 3")
    axes[0].legend(loc="lower right")
    axes[0].grid(alpha=0.2)

    axes[1].scatter(
        y[fit_mask], fitted[fit_mask], s=10, color="green", alpha=0.65, label="Training"
    )
    if has_holdout:
        axes[1].scatter(
            y[holdout_mask],
            fitted[holdout_mask],
            s=12,
            color="red",
            alpha=0.8,
            label="Holdout",
        )
    lower, upper = min(y.min(), fitted.min()), max(y.max(), fitted.max())
    axes[1].plot([lower, upper], [lower, upper], "k--", linewidth=1)
    axes[1].set_title("Observed vs. fitted")
    axes[1].set_xlabel("Column 3")
    axes[1].set_ylabel("Fit")
    axes[1].grid(alpha=0.2)
    if has_holdout:
        axes[1].legend(loc="upper left")

    figure.tight_layout()
    figure.savefig(args.output, dpi=160, bbox_inches="tight")


if __name__ == "__main__":
    main()
