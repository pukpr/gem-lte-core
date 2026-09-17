#!/usr/bin/env python3
"""Calculate Pearson correlation for a column in two lte_results.csv files.

Usage:
    ./compare_lte_results.py amo/lte_results.csv pdo/lte_results.csv model
    ./compare_lte_results.py amo/lte_results.csv pdo/lte_results.csv 4
    ./compare_lte_results.py amo/lte_results.csv pdo/lte_results.csv forcing --plot

Columns may be specified by their standard LTE result name or by a 1-based
column number: time, model, data, forcing, frequency, model_psd, data_psd.
Only finite values at corresponding timestamps in both files are compared.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


COLUMN_NAMES = (
    "time",
    "model",
    "data",
    "forcing",
    "frequency",
    "model_psd",
    "data_psd",
)


def parse_column(value: str) -> int:
    normalized = value.lower().replace("-", "_")
    if normalized in COLUMN_NAMES:
        return COLUMN_NAMES.index(normalized)

    try:
        column = int(value)
    except ValueError as exc:
        names = ", ".join(COLUMN_NAMES)
        raise argparse.ArgumentTypeError(
            f"column must be a number from 1 to {len(COLUMN_NAMES)} or one of: {names}"
        ) from exc

    if not 1 <= column <= len(COLUMN_NAMES):
        raise argparse.ArgumentTypeError(
            f"column number must be between 1 and {len(COLUMN_NAMES)}"
        )
    return column - 1


def read_results(path: Path, column: int) -> list[tuple[float, float]]:
    values: list[tuple[float, float]] = []
    with path.open(newline="", encoding="utf-8") as source:
        for line_number, row in enumerate(csv.reader(source), start=1):
            if not row or not any(field.strip() for field in row):
                continue
            if len(row) <= column:
                raise ValueError(
                    f"{path}:{line_number} has {len(row)} columns; "
                    f"column {column + 1} was requested"
                )
            try:
                timestamp = float(row[0])
                value = float(row[column])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number} contains a non-numeric value"
                ) from exc
            if values and timestamp <= values[-1][0]:
                raise ValueError(
                    f"{path}:{line_number} timestamp is not strictly increasing"
                )
            values.append((timestamp, value))
    return values


def matched_values(
    first: list[tuple[float, float]],
    second: list[tuple[float, float]],
    tolerance: float,
) -> tuple[list[float], list[float]]:
    left: list[float] = []
    right: list[float] = []
    first_index = 0
    second_index = 0

    while first_index < len(first) and second_index < len(second):
        first_time, first_value = first[first_index]
        second_time, second_value = second[second_index]
        if abs(first_time - second_time) <= tolerance:
            if math.isfinite(first_value) and math.isfinite(second_value):
                left.append(first_value)
                right.append(second_value)
            first_index += 1
            second_index += 1
        elif first_time < second_time:
            first_index += 1
        else:
            second_index += 1

    return left, right


def pearson_correlation(left: list[float], right: list[float]) -> float:
    if len(left) < 2:
        raise ValueError("at least two shared finite samples are required")

    left_mean = math.fsum(left) / len(left)
    right_mean = math.fsum(right) / len(right)
    covariance = math.fsum(
        (left_value - left_mean) * (right_value - right_mean)
        for left_value, right_value in zip(left, right)
    )
    left_sum_squares = math.fsum((value - left_mean) ** 2 for value in left)
    right_sum_squares = math.fsum((value - right_mean) ** 2 for value in right)
    denominator = math.sqrt(left_sum_squares * right_sum_squares)
    if denominator == 0:
        raise ValueError("correlation is undefined because one selected column is constant")
    return covariance / denominator


def plot_series(
    first: list[tuple[float, float]],
    second: list[tuple[float, float]],
    first_path: Path,
    second_path: Path,
    column: int,
    correlation: float,
) -> None:
    try:
        import matplotlib.pyplot as pyplot
    except ImportError as exc:
        raise RuntimeError("--plot requires matplotlib to be installed") from exc

    first_finite = [(time, value) for time, value in first if math.isfinite(value)]
    second_finite = [(time, value) for time, value in second if math.isfinite(value)]
    pyplot.plot(*zip(*first_finite), label=str(first_path))
    pyplot.plot(*zip(*second_finite), label=str(second_path))
    pyplot.xlabel("Time")
    pyplot.ylabel(COLUMN_NAMES[column])
    pyplot.title(f"{COLUMN_NAMES[column]} Pearson correlation: {correlation:.6g}")
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("first_file", type=Path)
    parser.add_argument("second_file", type=Path)
    parser.add_argument("column", type=parse_column)
    parser.add_argument("--plot", action="store_true", help="display both time series")
    parser.add_argument(
        "--time-tolerance",
        type=float,
        default=0.001,
        metavar="YEARS",
        help="maximum difference between corresponding timestamps (default: 0.001)",
    )
    args = parser.parse_args()
    if args.time_tolerance < 0:
        parser.error("--time-tolerance must be non-negative")

    try:
        first = read_results(args.first_file, args.column)
        second = read_results(args.second_file, args.column)
        left, right = matched_values(first, second, args.time_tolerance)
        correlation = pearson_correlation(left, right)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"column: {COLUMN_NAMES[args.column]} ({args.column + 1})")
    print(f"matched finite samples: {len(left)}")
    print(f"Pearson correlation: {correlation:.12g}")
    if args.plot:
        try:
            plot_series(
                first,
                second,
                args.first_file,
                args.second_file,
                args.column,
                correlation,
            )
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
