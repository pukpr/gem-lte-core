#!/usr/bin/env python3
"""Fit Singapore daily 30 hPa zonal wind from interpolated LTE hidden-layer values.

For each daily observation, LTE column 4 is linearly interpolated at the
day's midpoint in LTE's decimal-month time coordinate. The fitted model is:

    u = A + B*x + C*sin(D*x + E)

When applying fixed LTE column-3 coefficients, the monthly linear baseline is
interpolated to daily timestamps while the sine is evaluated from the daily
interpolated hidden layer. This retains the unresolved phase oscillations and
matches the monthly model exactly at LTE month-start anchors.

Example:
    python3 fit_singapore_qvbo_by_lte_hidden.py

Apply coefficients derived from LTE column 3 instead of fitting Singapore:
    python3 fit_singapore_qvbo_by_lte_hidden.py --coefficient-source lte-column3
"""

import argparse
import csv
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


LTE_COLUMN3_COEFFICIENTS = np.array(
    (
        -0.237212668547,
        0.130830150283,
        0.342717824310 * np.cos(2.196508156418),
        0.342717824310 * np.sin(2.196508156418),
    )
)
LTE_COLUMN3_FREQUENCY = 5062.151726330685


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lte-input",
        type=Path,
        default=Path("lte_results.csv"),
        help="Headerless LTE results file (default: lte_results.csv).",
    )
    parser.add_argument(
        "--qvbo-input",
        type=Path,
        default=Path("singapore_qvbo_30hpa_daily.csv"),
        help="Daily Singapore QVBO CSV (default: singapore_qvbo_30hpa_daily.csv).",
    )
    parser.add_argument(
        "--interpolated-output",
        type=Path,
        default=Path("singapore_qvbo_30hpa_lte_interpolated.csv"),
        help="Daily target and interpolated LTE column-4 values.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("singapore_qvbo_30hpa_lte_hidden_fit.png"),
        help="Output plot path.",
    )
    parser.add_argument(
        "--frequency-mode",
        choices=("plateau", "search"),
        default="plateau",
        help="Derive D from LTE monthly plateau steps or search D (default: plateau).",
    )
    parser.add_argument(
        "--coefficient-source",
        choices=("daily", "lte-column3"),
        default="daily",
        help=(
            "Fit coefficients to daily Singapore data, or fit LTE column 3 from "
            "column 4 and apply those fixed coefficients to Singapore (default: daily)."
        ),
    )
    parser.add_argument(
        "--max-frequency",
        type=float,
        default=50.0,
        help="Largest D to consider with --frequency-mode search (default: 50).",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=5001,
        help="D samples for --frequency-mode search (default: 5001).",
    )
    parser.add_argument(
        "--holdout-start",
        type=date.fromisoformat,
        help="First held-out daily observation, YYYY-MM-DD.",
    )
    parser.add_argument(
        "--holdout-end",
        type=date.fromisoformat,
        help="Last held-out daily observation, YYYY-MM-DD.",
    )
    return parser.parse_args()


def month_fraction(value: date):
    """Return a day's position within its calendar month."""
    if value.month == 12:
        next_month = date(value.year + 1, 1, 1)
    else:
        next_month = date(value.year, value.month + 1, 1)
    month_start = date(value.year, value.month, 1)
    return (value - month_start).days / (next_month - month_start).days


def dates_to_lte_time(values, lte_time):
    """Map dates onto the stored LTE monthly timestamp grid.

    LTE timestamps are rounded decimal months. Interpolating between adjacent
    stored values, rather than reconstructing year + month/12, makes every
    first-of-month timestamp exactly equal to its source LTE timestamp.
    """
    first_year = int(np.floor(lte_time[0]))
    first_month = int(round((lte_time[0] - first_year) * 12.0)) + 1
    mapped = []
    for value in values:
        index = (value.year - first_year) * 12 + value.month - first_month
        if index < 0 or index >= len(lte_time):
            raise ValueError(f"{value.isoformat()} is outside the LTE monthly range.")
        fraction = month_fraction(value)
        if index == len(lte_time) - 1:
            if fraction:
                raise ValueError(f"{value.isoformat()} is beyond the final LTE month.")
            mapped.append(lte_time[index])
        else:
            mapped.append(lte_time[index] + fraction * (lte_time[index + 1] - lte_time[index]))
    return np.array(mapped, dtype=float)


def load_qvbo(path: Path):
    dates = []
    wind = []
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {"date", "u_30hpa_ms"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain {', '.join(sorted(required))}.")
        for row in reader:
            dates.append(date.fromisoformat(row["date"]))
            wind.append(float(row["u_30hpa_ms"]))
    if not dates:
        raise ValueError(f"{path} contains no daily observations.")
    return np.array(dates, dtype=object), np.array(wind, dtype=float)


def fit_at_frequency(x, y, frequency):
    design = np.column_stack(
        (np.ones_like(x), x, np.sin(frequency * x), np.cos(frequency * x))
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coefficients
    return np.mean((y - fitted) ** 2), coefficients, fitted


def predict(x, coefficients, frequency):
    design = np.column_stack(
        (np.ones_like(x), x, np.sin(frequency * x), np.cos(frequency * x))
    )
    return design @ coefficients


def plateau_frequency(time, hidden):
    """Set D to 2*pi divided by the median low-magnitude monthly step."""
    order = np.argsort(time)
    time = time[order]
    hidden = hidden[order]
    monthly_pairs = np.isclose(np.diff(time), 1.0 / 12.0, rtol=0.0, atol=2e-6)
    increments = np.abs(np.diff(hidden)[monthly_pairs])
    if increments.size == 0:
        raise ValueError("No adjacent monthly LTE hidden-layer values are available.")
    plateau_steps = increments[increments <= np.quantile(increments, 0.25)]
    slope = np.median(plateau_steps)
    if slope <= 0:
        raise ValueError("The derived LTE plateau slope must be non-zero.")
    return slope, 2.0 * np.pi / slope


def find_best_fit(x, y, maximum_frequency, grid_size):
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
            frequency = minimize_scalar(
                lambda value: fit_at_frequency(x, y, value)[0],
                bounds=(frequencies[index - 1], frequencies[index + 1]),
                method="bounded",
                options={"xatol": 1e-12},
            ).x
        mse, coefficients, fitted = fit_at_frequency(x, y, frequency)
        if best is None or mse < best[0]:
            best = mse, frequency, coefficients, fitted
    return best


def metrics(observed, fitted):
    return np.corrcoef(observed, fitted)[0, 1], np.sqrt(np.mean((observed - fitted) ** 2))


def main():
    args = parse_args()
    if args.max_frequency <= 0 or args.grid_size < 3:
        raise ValueError("--max-frequency must be positive and --grid-size must be at least 3.")
    has_holdout = args.holdout_start is not None
    if has_holdout != (args.holdout_end is not None):
        raise ValueError("--holdout-start and --holdout-end must be supplied together.")
    if has_holdout and args.holdout_start > args.holdout_end:
        raise ValueError("--holdout-start must not be later than --holdout-end.")

    lte = np.loadtxt(args.lte_input, delimiter=",", ndmin=2)
    if lte.shape[1] < 4:
        raise ValueError(f"{args.lte_input} must contain at least four columns.")
    lte = lte[np.argsort(lte[:, 0])]
    lte_time = lte[:, 0]
    lte_hidden = lte[:, 3]

    observation_dates, observed = load_qvbo(args.qvbo_input)
    observation_time = dates_to_lte_time(observation_dates, lte_time)
    hidden = np.interp(observation_time, lte_time, lte_hidden)

    holdout = np.zeros(observed.shape, dtype=bool)
    if has_holdout:
        holdout = np.array(
            [args.holdout_start <= value <= args.holdout_end for value in observation_dates]
        )
        if holdout.sum() < 2 or (~holdout).sum() < 4:
            raise ValueError("The holdout must contain two rows and leave four fitting rows.")
    training = ~holdout

    if args.coefficient_source == "lte-column3":
        coefficients = LTE_COLUMN3_COEFFICIENTS
        frequency = LTE_COLUMN3_FREQUENCY
        slope = None
        coefficient_label = "LTE column-3 coefficients applied to Singapore"
    else:
        coefficient_x = hidden[training]
        coefficient_y = observed[training]
        coefficient_time = lte_time
        coefficient_label = "Singapore-fitted coefficients"
        if args.frequency_mode == "plateau":
            holdout_start_time, holdout_end_time = (
                dates_to_lte_time((args.holdout_start, args.holdout_end), lte_time)
                if has_holdout
                else (None, None)
            )
            source_training = ~(
                (coefficient_time >= holdout_start_time)
                & (coefficient_time <= holdout_end_time)
            ) if has_holdout else np.ones(coefficient_time.shape, dtype=bool)
            slope, frequency = plateau_frequency(
                lte_time[source_training], lte_hidden[source_training]
            )
            _, coefficients, fitted_training = fit_at_frequency(
                coefficient_x, coefficient_y, frequency
            )
        else:
            slope = None
            _, frequency, coefficients, fitted_training = find_best_fit(
                coefficient_x, coefficient_y, args.max_frequency, args.grid_size
            )
    if args.coefficient_source == "lte-column3":
        baseline = np.interp(
            observation_time,
            lte_time,
            coefficients[0] + coefficients[1] * lte_hidden,
        )
        # D times a plateau step is approximately 2*pi. Evaluate this term at
        # daily hidden-layer values to retain the intervening phase windings.
        fitted = baseline + coefficients[2] * np.sin(
            frequency * hidden
        ) + coefficients[3] * np.cos(frequency * hidden)
    else:
        fitted = predict(hidden, coefficients, frequency)
    correlation, rmse = metrics(observed[training], fitted[training])
    a, b, sine_coefficient, cosine_coefficient = coefficients
    amplitude = np.hypot(sine_coefficient, cosine_coefficient)
    phase = np.arctan2(cosine_coefficient, sine_coefficient)

    if slope is not None:
        print(f"Plateau monthly slope = {slope:.12f}")
        print("D = 2*pi / abs(plateau monthly slope)")
    print(f"A = {a:.12f}")
    print(f"B = {b:.12f}")
    print(f"C = {amplitude:.12f}")
    print(f"D = {frequency:.12f}")
    print(f"E = {phase:.12f} radians")
    print(f"Coefficient source = {args.coefficient_source}")
    print(f"Evaluation rows = {training.sum()}")
    print(f"Evaluation correlation coefficient = {correlation:.12f}")
    print(f"Evaluation RMSE = {rmse:.12f}")
    if has_holdout:
        holdout_correlation, holdout_rmse = metrics(observed[holdout], fitted[holdout])
        print(f"Holdout rows = {holdout.sum()}")
        print(f"Holdout correlation coefficient = {holdout_correlation:.12f}")
        print(f"Holdout RMSE = {holdout_rmse:.12f}")

    with args.interpolated_output.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination)
        writer.writerow(
            (
                "date",
                "lte_decimal_time",
                "lte_hidden_column4",
                "u_30hpa_ms",
                "model_value",
            )
        )
        for day, decimal_time, hidden_value, wind, model_value in zip(
            observation_dates, observation_time, hidden, observed, fitted
        ):
            writer.writerow(
                (
                    day.isoformat(),
                    f"{decimal_time:.12f}",
                    f"{hidden_value:.12f}",
                    f"{wind:.6f}",
                    f"{model_value:.12f}",
                )
            )

    plot_dates = np.array(observation_dates, dtype="datetime64[D]")
    figure, axes = plt.subplots(
        1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [2.5, 1]}, facecolor="#e6f2ff"
    )
    data_axis = axes[0]
    model_axis = data_axis.twinx()
    data_line = data_axis.plot(
        plot_dates,
        observed,
        color="blue",
        linewidth=0.8,
        alpha=0.7,
        label="Singapore 30 hPa",
    )
    model_line = model_axis.plot(
        plot_dates,
        fitted,
        color="red",
        linewidth=1.1,
        label="LTE hidden-layer model",
    )
    if has_holdout:
        data_axis.axvspan(
            np.datetime64(args.holdout_start),
            np.datetime64(args.holdout_end),
            color="black",
            alpha=0.08,
            label="Holdout",
        )
    data_axis.set_title(coefficient_label)
    data_axis.set_xlabel("Date")
    data_axis.set_ylabel("Singapore zonal wind at 30 hPa (m/s)", color="blue")
    model_units = (
        "Model output (LTE column-3 units)"
        if args.coefficient_source == "lte-column3"
        else "Model zonal wind at 30 hPa (m/s)"
    )
    model_axis.set_ylabel(model_units, color="red")
    data_axis.tick_params(axis="y", colors="blue")
    model_axis.tick_params(axis="y", colors="red")
    data_axis.grid(alpha=0.2)
    data_axis.legend(data_line + model_line, ["Singapore 30 hPa", "LTE hidden-layer model"], loc="lower right")

    axes[1].scatter(observed[training], fitted[training], s=8, alpha=0.5, color="green", label="Training")
    if has_holdout:
        axes[1].scatter(observed[holdout], fitted[holdout], s=10, alpha=0.7, color="red", label="Holdout")
    lower, upper = min(observed.min(), fitted.min()), max(observed.max(), fitted.max())
    axes[1].plot([lower, upper], [lower, upper], "k--", linewidth=1)
    axes[1].set_title(f"Evaluation r = {correlation:.4f}, RMSE = {rmse:.4f}")
    axes[1].set_xlabel("Observed Singapore wind (m/s)")
    axes[1].set_ylabel("Fitted wind (m/s)")
    axes[1].grid(alpha=0.2)
    if has_holdout:
        axes[1].legend(loc="upper left")
    figure.tight_layout()
    figure.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"Wrote {args.interpolated_output}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
