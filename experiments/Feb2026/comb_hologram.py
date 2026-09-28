#!/usr/bin/env python3
"""comb_hologram.py -- visualize a region's k_amp_phase terms AS a
frequency comb (harmonics of the shared backbone) rather than as a
scalogram, and test whether the comb's phases carry the signature of a
periodic IMPULSE TRAIN.

Motivation: winding_scalogram.py asks "is there a persistent SINGLE
winding at rate M" -- one M at a time, localized in time. For a region
like kPersianGulf whose real power is spread across many comparable-
amplitude harmonics of one shared backbone (confirmed directly: 10 of
kPersianGulf's 11 fitted k's are EXACT integer multiples of its 0.2075
backbone -- k/backbone lands within 0.001 of an integer every time), no
single M dominates any given time window, so the scalogram shows
neighboring teeth trading off local dominance as their relative phase
precesses -- the "moire" effect. This script instead looks at the whole
comb at once, the way a spectrum analyzer -- not a single-frequency
lock-in amplifier -- would.

Physical motivation for the comb itself: GEM-LTE's Impulse_Delta gates
the tidal forcing with a once-a-year DELTA (a spike, not a sinusoid --
see gem-lte-primitives-solution.adb). A perfectly periodic impulse train
has a Fourier series with ALL harmonics at comparable amplitude and a
PHASE THAT RAMPS LINEARLY with harmonic number (a shifted delta's n-th
Fourier coefficient is exp(-i*n*shift), i.e. phase = n*shift). If a
region's comb is genuinely an image of that impulse-gating mechanism
(filtered/reshaped by the subsequent IIR+Bessel stages, which will blur
but not destroy a linear phase ramp), its fitted phases should show a
detectable tendency to increase linearly with harmonic order n =
round(k/backbone) -- tested here with a circular ("Rayleigh-style")
regression, not just eyeballed.

Three panels:
  1. Comb teeth: amplitude vs harmonic order n (a real frequency-comb
     display, not a time-frequency map).
  2. Phase ramp: phase vs n, with the best-fit circular-linear ramp
     overlaid and its concentration R (1.0 = perfectly linear, 0 =
     random) and Rayleigh-test p-value reported directly -- the
     impulse-train test.
  3. The "hologram": each term plotted as a phasor amp*exp(i*phase) in
     the complex plane, connected in order of increasing n. A genuine
     linear phase ramp draws a smooth spiral; unrelated/random phases
     draw a zigzag. This is the direct geometric picture of whatever the
     phase-ramp test in panel 2 finds numerically.

Usage:
    comb_hologram.py <windings.json> [--backbone K] [--out FILE.png]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BACKBONE_K = 0.2075


def load_terms(path: Path) -> list[tuple[float, float, float]]:
    data = json.loads(path.read_text())
    return [(float(k), float(a), float(p)) for k, a, p in data["k_amp_phase"]]


def detect_backbone(terms: list[tuple[float, float, float]]) -> float:
    ks = [k for k, _, _ in terms]
    best = min(ks, key=lambda k: abs(k - BACKBONE_K))
    return best if abs(best - BACKBONE_K) < 0.05 else float(np.median(np.abs(ks)))


def circular_linear_fit(n: np.ndarray, phase: np.ndarray, n_grid: int = 4001):
    """Best-fit phase = slope*n + intercept (mod 2*pi), by direct grid
    search over slope (exact for this purpose -- n is small-integer and
    the objective is not convex, a closed-form circular regression is
    overkill). Returns (slope, intercept, resultant_length R, Rayleigh p)."""
    best = None
    for slope in np.linspace(-np.pi, np.pi, n_grid):
        resid = phase - slope * n
        c = np.angle(np.mean(np.exp(1j * resid)))
        score = np.sum(np.cos(resid - c))
        if best is None or score > best[0]:
            best = (score, slope, c)
    _, slope, intercept = best
    resid = np.angle(np.exp(1j * (phase - slope * n - intercept)))
    R = np.hypot(np.mean(np.cos(resid)), np.mean(np.sin(resid)))
    n_pts = len(n)
    z = n_pts * R * R
    p = np.exp(-z) * (1.0 + (2.0 * z - z * z) / (4.0 * n_pts))  # Rayleigh test, small-n correction
    return slope, intercept, R, max(0.0, min(1.0, p))


def run(path: Path, backbone: float | None, out: Path) -> None:
    terms = sorted(load_terms(path), key=lambda t: t[0])
    k1 = backbone if backbone is not None else detect_backbone(terms)

    rows = []
    for k, amp, phase in terms:
        n_raw = k / k1
        n = round(n_raw)
        if n <= 0:
            continue  # the near-zero/DC-ish term isn't part of the harmonic comb
        rows.append((n, k, amp, phase, n_raw - n))
    rows.sort(key=lambda r: r[0])

    print(f"[{path.parent.name}] backbone k1={k1:.4f}")
    print(f"{'n':>4} {'k':>9} {'resid':>7} {'amp':>8} {'phase':>8}")
    for n, k, amp, phase, resid in rows:
        print(f"{n:4d} {k:9.4f} {resid:7.3f} {amp:8.4f} {phase:8.3f}")

    n_arr = np.array([r[0] for r in rows], dtype=float)
    amp_arr = np.array([r[2] for r in rows])
    phase_arr = np.array([r[3] for r in rows])

    slope, intercept, R, p = circular_linear_fit(n_arr, phase_arr)
    print(f"\ncircular-linear phase-vs-n fit: slope={slope:.3f} rad/harmonic  "
          f"R={R:.3f}  Rayleigh p={p:.4f}")
    if p < 0.05:
        print("  -> significant linear phase ramp: consistent with a periodic "
              "impulse-train origin (Impulse_Delta), not independent random phases")
    else:
        print("  -> no significant linear ramp detected at this data volume")

    fig = plt.figure(figsize=(15, 5))
    ax1, ax2, ax3 = fig.subplots(1, 3)

    ax1.stem(n_arr, amp_arr, basefmt=" ")
    for n, k, amp, phase, resid in rows:
        ax1.annotate(f"{k:.3f}", (n, amp), fontsize=7, ha="center", va="bottom",
                     rotation=90, xytext=(0, 3), textcoords="offset points")
    ax1.set_xlabel("harmonic order n = round(k / backbone)")
    ax1.set_ylabel("amplitude")
    ax1.set_title(f"Comb teeth (backbone k1={k1:.4f})")

    ax2.scatter(n_arr, phase_arr, color="tab:blue", zorder=3, label="fitted phase")
    n_line = np.linspace(n_arr.min(), n_arr.max(), 200)
    ramp = np.angle(np.exp(1j * (slope * n_line + intercept)))
    # split the wrapped line into segments so it doesn't draw spurious
    # vertical jumps across the +/-pi seam
    jumps = np.where(np.abs(np.diff(ramp)) > np.pi)[0]
    segs = np.split(np.arange(len(n_line)), jumps + 1)
    for i, seg in enumerate(segs):
        ax2.plot(n_line[seg], ramp[seg], color="tab:red", lw=1.2,
                  label="best-fit linear ramp" if i == 0 else None)
    ax2.set_xlabel("harmonic order n")
    ax2.set_ylabel("phase (rad)")
    ax2.set_ylim(-np.pi, np.pi)
    ax2.set_title(f"Phase vs harmonic: R={R:.2f}, Rayleigh p={p:.3f}")
    ax2.legend(fontsize=8)

    phasors = amp_arr * np.exp(1j * phase_arr)
    ax3.plot(phasors.real, phasors.imag, "-", color="0.6", lw=1.0, zorder=1)
    sc = ax3.scatter(phasors.real, phasors.imag, c=n_arr, cmap="viridis", zorder=3)
    for n, z in zip(n_arr, phasors):
        ax3.annotate(f"n={int(n)}", (z.real, z.imag), fontsize=7,
                     xytext=(4, 4), textcoords="offset points")
    ax3.axhline(0, color="0.8", lw=0.6)
    ax3.axvline(0, color="0.8", lw=0.6)
    ax3.set_aspect("equal")
    ax3.set_xlabel("Re(amp * e^{i*phase})")
    ax3.set_ylabel("Im(amp * e^{i*phase})")
    ax3.set_title("Comb hologram: phasors connected in harmonic order\n"
                   "(spiral = linear phase ramp, zigzag = no ramp)")
    fig.colorbar(sc, ax=ax3, label="harmonic order n", fraction=0.046, pad=0.04)

    fig.suptitle(f"{path.parent.name}: harmonic comb structure "
                 f"({len(rows)} teeth, backbone k1={k1:.4f})", fontweight="bold")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\nsaved {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windings", type=Path, help="lt.exe.windings.json")
    ap.add_argument("--backbone", type=float, default=None,
                     help="override the auto-detected backbone k (default: "
                          "nearest term to the project's universal ~0.2075)")
    ap.add_argument("--out", type=Path, default=None,
                     help="output plot path (default: <region>_comb_hologram.png)")
    args = ap.parse_args()
    out = args.out if args.out else Path(f"{args.windings.parent.name}_comb_hologram.png")
    run(args.windings, args.backbone, out)


if __name__ == "__main__":
    main()
