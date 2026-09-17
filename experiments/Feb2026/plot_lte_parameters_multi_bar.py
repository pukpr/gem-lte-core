#!/usr/bin/env python3
"""Plot fitted LTE parameters across the standard result directories."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SERIES_DIRECTORIES = (
    "nino4",
#    "emi",
#    "111",
    "pdo",
    "amo",
    "nao",
    "tna",
    "iode",
    "baltic"
#    "88",
#    "155",
#    "245",
)

# "offs" is the parameter-file name for the requested offset parameter.
PARAMETERS = (
    ("offset", ("offset", "offs")),
    ("impA", ("impA",)),
    ("impB", ("impB",)),
    ("delA", ("delA",)),
    ("delB", ("delB",)),
    ("asym", ("asym",)),
    ("ma", ("ma",)),
    ("mp", ("mp",)),
    ("shfT", ("shfT",)),
    ("init", ("init",)),
    ("dd", ("dd",)),
    ("fbS", ("fbS",)),
    ("fbC", ("fbC",)),
)


def load_parameters(directory: Path) -> dict[str, float]:
    path = directory / "lt.exe.p"
    if not path.is_file():
        raise FileNotFoundError(f"Missing parameter file: {path}")

    with path.open(encoding="utf-8") as parameter_file:
        data = json.load(parameter_file)
    if not isinstance(data, dict):
        raise TypeError(f"{path} must contain a JSON object")

    values: dict[str, float] = {}
    for label, keys in PARAMETERS:
        key = next((candidate for candidate in keys if candidate in data), None)
        if key is None:
            raise KeyError(f"{path} is missing parameter {label!r}")
        value = data[key]
        if not isinstance(value, (int, float)):
            raise TypeError(f"{path}: parameter {key!r} must be numeric")
        values[label] = float(value)
    return values


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot fitted LTE parameters from offset through fbC as grouped bars."
    )
    parser.add_argument(
        "subdirs",
        nargs="*",
        default=SERIES_DIRECTORIES,
        help="result directories to plot (default: the standard 12 directories)",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="use a linear x-axis instead of the default symmetric-log axis",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("lte_parameters_multi_bar.png"),
        help="output image path (default: %(default)s)",
    )
    parser.add_argument("--show", action="store_true", help="display the chart interactively")
    args = parser.parse_args()

    labels = [label for label, _ in PARAMETERS]
    values_by_directory = {
        directory: load_parameters(Path(directory)) for directory in args.subdirs
    }

    figure, axis = plt.subplots(figsize=(14, 8), constrained_layout=True)
    positions = np.arange(len(labels))
    bar_height = 0.8 / len(values_by_directory)
    colors = plt.get_cmap("tab20").colors

    for index, (directory, values) in enumerate(values_by_directory.items()):
        offset = (index - (len(values_by_directory) - 1) / 2) * bar_height
        axis.barh(
            positions + offset,
            [values[label] for label in labels],
            height=bar_height,
            color=colors[index],
            label=directory,
        )

    axis.axvline(0, color="black", linewidth=0.8)
    axis.set_yticks(positions, labels)
    axis.invert_yaxis()
    axis.set_xlabel("Fitted parameter value")
    axis.set_ylabel("Parameter")
    axis.set_title("LTE Parameters: offset through fbC")
    axis.grid(True, axis="x", linestyle=":", alpha=0.5)
    if not args.linear:
        axis.set_xscale("symlog", linthresh=1e-4)
        axis.set_xlabel("Fitted parameter value (symmetric log scale)")
    axis.legend(title="Subdirectory", ncols=2, loc="lower right")

    figure.savefig(args.output, dpi=200)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
