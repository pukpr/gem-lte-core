#!/usr/bin/env python3
"""brestexcl_massage_overlay_check.py — direct visual overlay of the
massage.py reconstruction (13pt boxcar + quadratic detrend, row-index
date convention) against real brestexcl.dat, to look for the shape of
the discrepancy (r=0.56) rather than just its magnitude."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare  # noqa: E402

ROOT = Path(__file__).resolve().parent


def boxcar_filter(values, window_size=13):
    n = len(values)
    filtered = np.zeros(n)
    half = window_size // 2
    for i in range(n):
        s, e = max(0, i - half), min(n, i + half + 1)
        filtered[i] = values[s:e].mean()
    return filtered


def main() -> int:
    times, values, full_idx = [], [], []
    with open(ROOT / "brestexcl" / "psmsl_brest_1.rlrdata") as f:
        for i, line in enumerate(f):
            t, v = float(line.split(";")[0]), float(line.split(";")[1])
            if v == -99999:
                continue
            times.append(t)
            values.append(v)
            full_idx.append(i)
    times, values, full_idx = np.array(times), np.array(values), np.array(full_idx)
    proj_dates = 1807.0 + full_idx / 12.0

    filtered = boxcar_filter(values, 13)
    p = np.polyfit(times, filtered, 2)
    massaged = filtered - np.polyval(p, times)

    prep = prepare("brestexcl")
    dates_orig, data_raw = prep["dates"], prep["data_raw"]
    nz = data_raw != 0.0

    recon_by_month = {round(d * 12): v for d, v in zip(proj_dates, massaged)}
    d_match, b_match, r_match = [], [], []
    for do, bv in zip(dates_orig[nz], data_raw[nz]):
        key = round(do * 12)
        if key in recon_by_month:
            d_match.append(do)
            b_match.append(bv)
            r_match.append(recon_by_month[key])
    d_match, b_match, r_match = np.array(d_match), np.array(b_match), np.array(r_match)

    fig, axes = plt.subplots(4, 1, figsize=(13, 13))

    axes[0].plot(d_match, b_match, color="0.25", linewidth=0.8, label="real brestexcl.dat")
    axes[0].plot(d_match, r_match, color="teal", linewidth=0.8, alpha=0.8,
                 label="massage.py reconstruction (13pt boxcar + quad detrend)")
    axes[0].set_title(f"full overlap, 1880-2020 (r={np.corrcoef(r_match,b_match)[0,1]:+.3f})",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    zmask = (d_match >= 1950) & (d_match < 1965)
    axes[1].plot(d_match[zmask], b_match[zmask], color="0.25", linewidth=1.2,
                 marker="o", markersize=3, label="real brestexcl.dat")
    axes[1].plot(d_match[zmask], r_match[zmask], color="teal", linewidth=1.2,
                 marker="o", markersize=3, alpha=0.8, label="reconstruction")
    axes[1].set_title("zoom: 1950-1965", fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")

    zmask2 = (d_match >= 2000) & (d_match < 2015)
    axes[2].plot(d_match[zmask2], b_match[zmask2], color="0.25", linewidth=1.2,
                 marker="o", markersize=3, label="real brestexcl.dat")
    axes[2].plot(d_match[zmask2], r_match[zmask2], color="teal", linewidth=1.2,
                 marker="o", markersize=3, alpha=0.8, label="reconstruction")
    axes[2].set_title("zoom: 2000-2015", fontsize=10)
    axes[2].legend(fontsize=8, loc="upper left")

    axes[3].scatter(r_match, b_match, s=4, alpha=0.3, color="0.3")
    lims = [min(r_match.min(), b_match.min()), max(r_match.max(), b_match.max())]
    axes[3].plot(lims, lims, color="crimson", linewidth=1, linestyle="--", label="y=x")
    axes[3].set_xlabel("reconstruction")
    axes[3].set_ylabel("real brestexcl.dat")
    axes[3].set_title("scatter: reconstruction vs real brestexcl.dat", fontsize=10)
    axes[3].legend(fontsize=8)

    fig.tight_layout()
    out = ROOT / "brestexcl" / "brestexcl_massage_overlay_check.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[saved] {out}")

    # a couple of extra numeric clues worth having alongside the plot
    print(f"skewness check: recon skew={float(((r_match-r_match.mean())**3).mean()/r_match.std()**3):+.3f}  "
          f"brestexcl skew={float(((b_match-b_match.mean())**3).mean()/b_match.std()**3):+.3f}")
    print(f"recon range=[{r_match.min():.1f},{r_match.max():.1f}]  "
          f"brestexcl range=[{b_match.min():.1f},{b_match.max():.1f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
