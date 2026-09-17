#!/usr/bin/env python3
"""Plot daily Singapore/Changi QVBO proxy data extracted at 30 hPa.

Example:
    python3 plot_singapore_qvbo_30hpa.py
"""

import argparse
import csv
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot daily Singapore QVBO zonal wind at 30 hPa."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("singapore_qvbo_30hpa_daily.csv"),
        help="Daily CSV from extract_singapore_qvbo_30hpa.sh.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("singapore_qvbo_30hpa_daily.png"),
        help="Output plot path (default: singapore_qvbo_30hpa_daily.png).",
    )
    parser.add_argument(
        "--start",
        type=date.fromisoformat,
        help="First date to plot, in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--end",
        type=date.fromisoformat,
        help="Last date to plot, in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--rolling-days",
        type=int,
        default=30,
        help="Calendar-day rolling-mean window (default: 30).",
    )
    return parser.parse_args()


def load_daily_wind(path: Path):
    days = []
    wind = []
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        expected_columns = {"date", "u_30hpa_ms", "soundings"}
        if reader.fieldnames is None or not expected_columns.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain: {', '.join(sorted(expected_columns))}")
        for row in reader:
            days.append(date.fromisoformat(row["date"]))
            wind.append(float(row["u_30hpa_ms"]))
    if not days:
        raise ValueError(f"{path} contains no daily observations.")
    return np.array(days, dtype="datetime64[D]"), np.array(wind, dtype=float)


def calendar_rolling_mean(days, values, window_days):
    """Return a trailing calendar-day mean that ignores days without observations."""
    start = days.min()
    end = days.max()
    all_days = np.arange(start, end + np.timedelta64(1, "D"), dtype="datetime64[D]")
    positions = (days - start).astype(int)
    complete = np.full(all_days.size, np.nan)
    complete[positions] = values

    valid = np.isfinite(complete)
    kernel = np.ones(window_days)
    totals = np.convolve(np.where(valid, complete, 0.0), kernel, mode="full")[
        : all_days.size
    ]
    counts = np.convolve(valid.astype(int), kernel, mode="full")[: all_days.size]
    return all_days, np.divide(totals, counts, out=np.full_like(totals, np.nan), where=counts > 0)


def main():
    args = parse_args()
    if args.rolling_days < 1:
        raise ValueError("--rolling-days must be at least one.")
    if args.start and args.end and args.start > args.end:
        raise ValueError("--start must not be later than --end.")

    days, wind = load_daily_wind(args.input)
    rolling_days, rolling_wind = calendar_rolling_mean(days, wind, args.rolling_days)
    plot_mask = np.ones(days.shape, dtype=bool)
    rolling_mask = np.ones(rolling_days.shape, dtype=bool)
    if args.start:
        start = np.datetime64(args.start)
        plot_mask &= days >= start
        rolling_mask &= rolling_days >= start
    if args.end:
        end = np.datetime64(args.end)
        plot_mask &= days <= end
        rolling_mask &= rolling_days <= end
    if not np.any(plot_mask):
        raise ValueError("The requested date range contains no observations.")

    figure, axis = plt.subplots(figsize=(16, 7), facecolor="#e6f2ff")
    axis.plot(
        days[plot_mask],
        wind[plot_mask],
        color="steelblue",
        linewidth=0.7,
        alpha=0.6,
        label="Daily mean",
    )
    axis.plot(
        rolling_days[rolling_mask],
        rolling_wind[rolling_mask],
        color="crimson",
        linewidth=1.4,
        label=f"{args.rolling_days}-day mean",
    )
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_title("Singapore/Changi QVBO proxy: zonal wind at 30 hPa")
    axis.set_xlabel("Date")
    axis.set_ylabel("Eastward zonal wind (m/s)")
    axis.grid(alpha=0.25)
    axis.legend(loc="upper right")
    figure.tight_layout()
    figure.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
