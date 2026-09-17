#!/usr/bin/env python3
"""Plot hidden latent forcing for the standard LTE result directories."""

from __future__ import annotations

import argparse
from bisect import bisect_right
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


SERIES_DIRECTORIES = (
    "nino4",
    "pdo",
    "amo",
    "nao",
    "iode",
    "baltic",
    "tna"
#    "emi"
)

INSET_ALIGNMENTS = {
    "ul": (0.02, 0.74),
    "uc": (0.41, 0.74),
    "ur": (0.80, 0.74),
    "cl": (0.02, 0.39),
    "cc": (0.41, 0.39),
    "cr": (0.80, 0.39),
    "ll": (0.02, 0.04),
    "lc": (0.41, 0.04),
    "lr": (0.80, 0.04),
}
INSET_SIZE = (0.18, 0.22)


def read_time_and_forcing(path: Path) -> tuple[list[float], list[float]]:
    """Read decimal time (column 1) and hidden latent forcing (column 4)."""
    times: list[float] = []
    forcings: list[float] = []

    with path.open(newline="", encoding="utf-8") as csv_file:
        for line_number, row in enumerate(csv.reader(csv_file, skipinitialspace=True), start=1):
            if not row:
                continue
            if len(row) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least four columns")
            try:
                times.append(float(row[0]))
                forcings.append(float(row[3]))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid numeric value") from error

    if not times:
        raise ValueError(f"{path}: contains no data rows")
    return times, forcings


def normalize_excursions(values: list[float]) -> list[float]:
    """Remove the DC offset and scale zero-mean excursions to unit RMS."""
    mean = sum(values) / len(values)
    excursions = [value - mean for value in values]
    rms = math.sqrt(sum(value * value for value in excursions) / len(excursions))
    if rms == 0.0:
        raise ValueError("Cannot normalize a series with zero RMS excursions")
    return [value / rms for value in excursions]


def pearson_correlation(left: list[float], right: list[float]) -> float:
    """Return the Pearson correlation over the shared prefix of two series."""
    sample_count = min(len(left), len(right))
    if sample_count < 2:
        raise ValueError("At least two shared samples are required for correlation")

    left_values = left[:sample_count]
    right_values = right[:sample_count]
    left_mean = sum(left_values) / sample_count
    right_mean = sum(right_values) / sample_count
    covariance = sum(
        (left_value - left_mean) * (right_value - right_mean)
        for left_value, right_value in zip(left_values, right_values)
    )
    left_variance = sum((value - left_mean) ** 2 for value in left_values)
    right_variance = sum((value - right_mean) ** 2 for value in right_values)
    denominator = math.sqrt(left_variance * right_variance)
    if denominator == 0.0:
        raise ValueError("Cannot correlate a series with zero variance")
    return covariance / denominator


def add_correlation_inset(
    axis: plt.Axes,
    series: list[tuple[str, list[float], list[float]]],
    alignment: str,
) -> None:
    """Add a pairwise Pearson-correlation matrix inset to the comparison plot."""
    labels = [name for name, _, _ in series]
    forcings = [forcing for _, _, forcing in series]
    correlations = [
        [pearson_correlation(left, right) for right in forcings] for left in forcings
    ]

    inset = axis.inset_axes((*INSET_ALIGNMENTS[alignment], *INSET_SIZE))
    image = inset.imshow(correlations, vmin=-1, vmax=1, cmap="coolwarm")
    inset.set_xticks(range(len(labels)), labels, rotation=45, ha="right", fontsize=6)
    inset.set_yticks(range(len(labels)), labels, fontsize=6)
    inset.set_title("Pearson r", fontsize=6)
    for row, correlation_row in enumerate(correlations):
        for column, correlation in enumerate(correlation_row):
            inset.text(
                column,
                row,
                f"{correlation:.3f}",
                ha="center",
                va="center",
                fontsize=4,
                color="white" if abs(correlation) > 0.55 else "black",
            )
    inset.figure.colorbar(
        image, ax=inset, location="right", fraction=0.05, pad=0.01, shrink=0.7
    )


def average_normalized_waveform(
    series: list[tuple[str, list[float], list[float]]],
) -> tuple[list[float], list[float]]:
    """Return a monthly ensemble average across the time shared by all series."""
    start = max(time[0] for _, time, _ in series)
    stop = min(time[-1] for _, time, _ in series)
    if stop <= start:
        raise ValueError("Series do not share a common time range")

    sample_count = int(math.floor((stop - start) * 12)) + 1
    times = [start + index / 12 for index in range(sample_count)]
    normalized = [(time, normalize_excursions(forcing)) for _, time, forcing in series]
    average: list[float] = []

    for point in times:
        values: list[float] = []
        for time, forcing in normalized:
            index = bisect_right(time, point) - 1
            if index < 0 or index >= len(time) - 1:
                raise ValueError(f"Cannot interpolate average waveform at {point:.3f}")
            fraction = (point - time[index]) / (time[index + 1] - time[index])
            values.append(forcing[index] + fraction * (forcing[index + 1] - forcing[index]))
        average.append(sum(values) / len(values))
    return times, average


def moving_average(values: list[float], window: int) -> list[float]:
    radius = window // 2
    return [
        sum(values[max(0, index - radius) : min(len(values), index + radius + 1)])
        / (min(len(values), index + radius + 1) - max(0, index - radius))
        for index in range(len(values))
    ]


def cycle_peak_pair(
    times: list[float], values: list[float], period_years: float
) -> tuple[float, float]:
    """Find the strongest adjacent peak pair separated by the target period."""
    smoothed = moving_average(values, max(3, round(period_years * 3)))
    minimum_separation = round(period_years * 12 * 0.65)
    peak_indices: list[int] = []

    for index in range(1, len(smoothed) - 1):
        if smoothed[index] < smoothed[index - 1] or smoothed[index] <= smoothed[index + 1]:
            continue
        if not peak_indices or index - peak_indices[-1] >= minimum_separation:
            peak_indices.append(index)
        elif smoothed[index] > smoothed[peak_indices[-1]]:
            peak_indices[-1] = index

    tolerance = period_years * 0.25
    pairs = [
        (times[left], times[right], smoothed[left] + smoothed[right])
        for left, right in zip(peak_indices, peak_indices[1:])
        if abs((times[right] - times[left]) - period_years) <= tolerance
    ]
    if not pairs:
        raise ValueError(f"No peak intervals found near the {period_years:g}-year cycle")
    start, stop, _ = max(pairs, key=lambda pair: pair[2])
    return start, stop


def annotate_cycles(axis: plt.Axes, series: list[tuple[str, list[float], list[float]]]) -> None:
    """Mark peak-to-peak intervals for the target cycles of the average waveform."""
    times, average = average_normalized_waveform(series)
    cycle_definitions = (
        (3.8, 0.05, "3.8-year cycle", "tab:purple"),
        (18.6, 0.13, "18.6-year beat cycle", "tab:brown"),
    )

    for period, y_position, label, color in cycle_definitions:
        start, stop = cycle_peak_pair(times, average, period)
        transform = axis.get_xaxis_transform()
        axis.vlines(
            (start, stop),
            y_position - 0.012,
            y_position + 0.012,
            color=color,
            linewidth=1.2,
            transform=transform,
        )
        axis.annotate(
            "",
            xy=(stop, y_position),
            xytext=(start, y_position),
            xycoords=transform,
            arrowprops={"arrowstyle": "<->", "color": color, "lw": 1.2},
        )
        axis.text(
            (start + stop) / 2,
            y_position + 0.018,
            label,
            color=color,
            ha="center",
            va="bottom",
            transform=transform,
            fontsize=9,
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare hidden latent forcing across LTE result directories."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("hidden_latent_forcing_comparison.png"),
        help="output image path (default: %(default)s)",
    )
    parser.add_argument("--show", action="store_true", help="display the chart interactively")
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="remove each series DC offset and scale its excursions to unit RMS",
    )
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="mark peak-detected 3.8-year and 18.6-year cycles of the average waveform",
    )
    parser.add_argument(
        "--align",
        choices=INSET_ALIGNMENTS,
        default="uc",
        help="Pearson r inset alignment (default: %(default)s)",
    )
    args = parser.parse_args()

    figure, axis = plt.subplots(figsize=(14, 8), constrained_layout=True)
    pdo_axis = None if args.normalize else axis.twinx()
    colors = plt.get_cmap("tab20").colors
    series: list[tuple[str, list[float], list[float]]] = []

    for index, directory in enumerate(SERIES_DIRECTORIES):
        csv_path = Path(directory) / "lte_results.csv"
        if not csv_path.is_file():
            raise FileNotFoundError(f"Missing results file: {csv_path}")
        time, forcing = read_time_and_forcing(csv_path)
        series.append((directory, time, forcing))
        if args.normalize:
            forcing = normalize_excursions(forcing)
        target_axis = pdo_axis if directory == "pdo" and pdo_axis is not None else axis
        target_axis.plot(time, forcing, label=directory, color=colors[index], linewidth=1.1)

    axis.set_title("Hidden Latent Forcing Comparison")
    axis.set_xlabel("Decimal time")
    axis.set_ylabel(
        "Normalized hidden latent forcing (unit RMS)"
        if args.normalize
        else "Hidden latent forcing"
    )
    if pdo_axis is not None:
        pdo_axis.set_ylabel("PDO hidden latent forcing")
    axis.grid(True, alpha=0.3)
    add_correlation_inset(axis, series, args.align)
    if args.annotate:
        annotate_cycles(axis, series)
    lines = axis.get_lines() if pdo_axis is None else axis.get_lines() + pdo_axis.get_lines()
    axis.legend(lines, [line.get_label() for line in lines], title="Series", ncols=2, loc="upper right")
    figure.savefig(args.output, dpi=200)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
