#!/usr/bin/env python3
"""bessel_order_diff.py -- isolate and plot the two Bessel sideband terms'
individual contributions to the forcing manifold (lte_results.csv column
4), for a given index.

The Ada Bessel step (gem-lte-primitives-solution.adb) is:

    M(I) := M(I) + eS*sin(2*pi*k*M(I)) + eC*cos(2*pi*k*M(I))     -- 1st order
                  + eS2*eS*sin(4*pi*k*M(I)) + eC2*eC*cos(4*pi*k*M(I))  -- 2nd order

applied to M = the raw IIR/comb output (BEFORE any Bessel step), with
eS=ImpA, eC=ImpB, k=the backbone winding rate, eS2=Offset, eC2=bg. The
"2nd order" term is literally the 1st order's own amplitude re-scaled by
eS2/eC2 and evaluated at DOUBLE the carrier's argument-rate (4*pi*k vs
2*pi*k) -- an FM sideband structure, not an independent free term (k2 is
passed into the Ada function but the body never reads it -- confirmed
dead by the compiler's own unused-parameter warning).

This script reconstructs three variants of the forcing pipeline up to
that point:
    R  = raw IIR output, no Bessel step at all
    F1 = R with ONLY the 1st-order Bessel term added (eS2=eC2=0)
    F2 = R with BOTH orders added (the actual production forcing --
         matches lte_results.csv column 4 / ws.load_index's own
         prep["forcing"] exactly)

and plots the two individual ADDITIONS:
    F1 - R   (pure 1st-order contribution)
    F2 - F1  (pure 2nd-order contribution, incremental)
    F2 - R   (their sum, i.e. the full Bessel step's total effect)

Usage:
    bessel_order_diff.py <index> [--root DIR] [--out FILE.png]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ws  # noqa: E402


def run(index: str, root: str | None, out: Path) -> None:
    prep = ws.load_index(ws.find_index_root(index, root), index)
    dates = prep["dates"]
    params = prep["params"]

    lpap = np.array(params["lpap"], dtype=float)
    periods = ws.doodson_periods(prep["year_startup"], prep["year_cand"])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp", "shfT", "init")}
    k1 = float(np.array(params["ltep"], dtype=float)[prep["nm"] - 1])

    tf = ws.tide_sum(dates, lpap[:, 1:3], periods, prep["yl"], 0.0, B["shfT"])
    comb = ws.impulse_delta(dates, B["delA"], B["delB"], B["asym"], 12)
    R = ws.iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
               start_date=prep["idate"], dates=dates)

    F1 = ws.bessel(R, B["impA"], B["impB"], k1, 0.0, 0.0)
    F2 = ws.bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])

    # sanity check: F2 must reproduce prep["forcing"] (lte_results.csv col 4) exactly
    max_check = float(np.max(np.abs(F2 - prep["forcing"])))
    print(f"[{index}] sanity check |F2 - prep['forcing']| max = {max_check:.3e} "
          f"(should be ~0)")
    print(f"[{index}] k1={k1:.4f}  impA={B['impA']:.4f}  impB={B['impB']:.4f}  "
          f"offs(eS2)={B['offs']:.4f}  bg(eC2)={B['bg']:.4f}")

    d1 = F1 - R      # pure 1st-order addition
    d2 = F2 - F1     # pure 2nd-order addition (incremental)
    dtot = F2 - R    # total Bessel addition

    rms_R = float(np.sqrt(np.mean(R ** 2)))
    for name, d in [("1st-order addition (F1-R)", d1),
                    ("2nd-order addition (F2-F1)", d2),
                    ("total Bessel addition (F2-R)", dtot)]:
        rms = float(np.sqrt(np.mean(d ** 2)))
        print(f"  {name:32s} RMS={rms:.5f}  ({100*rms/rms_R:.2f}% of raw-IIR RMS)")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    ax1.plot(dates, R, color="0.5", lw=0.8, label="R (raw IIR, no Bessel)")
    ax1.plot(dates, F1, color="tab:blue", lw=0.8, label="F1 (+ 1st order)")
    ax1.plot(dates, F2, color="tab:red", lw=0.8, alpha=0.8, label="F2 (+ 1st + 2nd order, = production Forcing)")
    ax1.set_ylabel("forcing manifold value")
    ax1.set_title(f"{index}: forcing manifold (lte_results.csv col 4) with/without Bessel sidebands")
    ax1.legend(fontsize=8)

    ax2.plot(dates, d1, color="tab:blue", lw=0.9, label="1st-order addition  (F1 - R)")
    ax2.plot(dates, d2, color="tab:green", lw=0.9, label="2nd-order addition  (F2 - F1)")
    ax2.plot(dates, dtot, color="tab:red", lw=0.7, alpha=0.6, label="total Bessel addition  (F2 - R)")
    ax2.axhline(0.0, color="0.7", lw=0.6)
    ax2.set_xlabel("year")
    ax2.set_ylabel("difference")
    ax2.set_title("Isolated Bessel contributions")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\nsaved {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index")
    ap.add_argument("--root", default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    out = args.out if args.out else Path(f"{args.index}_bessel_order_diff.png")
    run(args.index, args.root, out)


if __name__ == "__main__":
    main()
