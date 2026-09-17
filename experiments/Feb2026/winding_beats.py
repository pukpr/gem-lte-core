#!/usr/bin/env python3
"""Ridge detection + beat-interference analysis on winding scalograms.

Algorithm (as proposed):
  1. Compute the winding scalogram (see winding_scalogram.py) of a signal
     against a Forcing manifold.
  2. Identify strong, *unbroken* horizontal bands — winding numbers M where
     power stays above the matched-noise floor across nearly the whole
     record, not just a transient blip.
  3. Where two such bands sit close enough together that they're hard to
     resolve as two separate ridges, treat them as a candidate
     near-degenerate pair and look for a beat signature: since both terms
     share the same clock (Forcing(t)), the exact instantaneous envelope of
     their sum is

         R(t) = sqrt(A1^2 + A2^2 + 2*A1*A2*cos(dtheta(t))),
         dtheta(t) = 2*pi*(M1-M2)*Forcing(t) + (phase1-phase2)

     — an exact identity, not an approximation. R(t) -> A1+A2 wherever the
     two components constructively interfere, R(t) -> |A1-A2| wherever they
     cancel. Plotting R(t) directly shows where the two nearby windings
     amplify or cancel each other over the record.

Two runs:
  (a) PDO — a worked example using PDO's own regression fit (ground truth
      M/amplitude/phase from compute_winding), reproducing and visualizing
      the PDO/NINO4 near-degenerate pair already described in
      WINDING_NARRATIVE.md.
  (b) nino34 — a genuinely out-of-sample test: nino34/nino34.dat's raw data
      (never separately fit for this purpose) run against nino4's own
      already-calibrated Forcing manifold (column 4 of nino4/lte_results.csv)
      — no per-index optimization at all, just reusing a manifold fit to a
      *different* dataset.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns, standardize
from winding_scalogram import winding_transform, noise_floor
from param_survey import compute_winding

ROOT = Path(__file__).resolve().parent
M_MAX, DM, SIGMA, T0_STEP = 5.0, 0.01, 15.0, 5.0


def compute_snr_grid(year, x, forcing, m_grid, t0_grid, sigma):
    G, edge_mask, _ = winding_transform(year, x, forcing, m_grid, t0_grid, sigma)
    floor = noise_floor(year, forcing, m_grid, t0_grid, sigma)
    snr = np.abs(G) ** 2 / floor[:, None]
    return np.log2(np.maximum(snr, 1e-6)), edge_mask


def detect_ridges(m_grid, log_snr, edge_mask, min_prominence=1.0,
                   min_persistence=0.6, min_distance_pts=3):
    """Strong = high mean SNR; unbroken = above the noise floor most of the
    (non-edge-affected) record. Returns indices into m_grid."""
    non_edge = ~edge_mask
    strength = log_snr[:, non_edge].mean(axis=1)
    persistence = (log_snr[:, non_edge] > 0).mean(axis=1)
    peaks, _ = find_peaks(strength, prominence=min_prominence,
                           distance=min_distance_pts)
    keep = np.array([p for p in peaks if persistence[p] >= min_persistence],
                     dtype=int)
    return keep, strength, persistence


def beat_envelope(forcing, m1, amp1, phase1, m2, amp2, phase2):
    dtheta = 2.0 * np.pi * (m1 - m2) * forcing + (phase1 - phase2)
    R = np.sqrt(amp1 ** 2 + amp2 ** 2 + 2.0 * amp1 * amp2 * np.cos(dtheta))
    return R, dtheta


def closest_pair(m, amp, phase):
    absm = np.abs(m)
    order = np.argsort(absm)
    gaps = np.diff(absm[order])
    j = np.argmin(gaps)
    i1, i2 = order[j], order[j + 1]
    return (m[i1], amp[i1], phase[i1]), (m[i2], amp[i2], phase[i2])


def pair_near(m, amp, phase, target1, target2):
    """Pick the two fitted components nearest two given M targets — used to
    reproduce a *specific* known pair rather than whichever gap happens to
    be globally smallest."""
    absm = np.abs(m)
    i1 = int(np.argmin(np.abs(absm - target1)))
    i2 = int(np.argmin(np.abs(absm - target2)))
    return (m[i1], amp[i1], phase[i1]), (m[i2], amp[i2], phase[i2])


def plot_ridge_profile(ax, m_grid, strength, persistence, ridges,
                        highlight=None, title=""):
    ax2 = ax.twinx()
    ax.plot(m_grid, strength, color="tab:blue", linewidth=1.0)
    ax.scatter(m_grid[ridges], strength[ridges], color="red", s=18,
               zorder=5, label="detected ridge")
    ax2.plot(m_grid, persistence, color="tab:gray", linewidth=0.8,
              alpha=0.6, linestyle=":")
    if highlight is not None:
        for m in highlight:
            ax.axvline(abs(m), color="green", linestyle="--", linewidth=1.0,
                       alpha=0.7)
    ax.set_ylabel("mean log2 SNR (strength)", color="tab:blue")
    ax2.set_ylabel("time fraction above noise floor", color="tab:gray")
    ax.set_xlabel("winding number M")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)


def plot_beat(ax_top, ax_bot, year, forcing, m1, a1, p1, m2, a2, p2,
              actual=None, actual_label=""):
    s1 = a1 * np.sin(2.0 * np.pi * m1 * forcing + p1)
    s2 = a2 * np.sin(2.0 * np.pi * m2 * forcing + p2)
    ssum = s1 + s2
    R, dtheta = beat_envelope(forcing, m1, a1, p1, m2, a2, p2)
    c = np.cos(dtheta)

    ax_top.plot(year, ssum, color="black", linewidth=0.7,
                label=f"reconstructed sum: M1={m1:.4f}+M2={m2:.4f}")
    ax_top.plot(year, R, color="tab:red", linewidth=1.2,
                label="envelope R(t) (exact)")
    ax_top.plot(year, -R, color="tab:red", linewidth=1.2)
    ax_top.fill_between(year, -R, R, color="tab:red", alpha=0.08)
    ax_top.legend(loc="upper right", fontsize=7)
    ax_top.set_ylabel("reconstructed amplitude")

    # constructive (c near +1) vs destructive (c near -1) shading
    ax_bot.fill_between(year, 0, 1, where=c > 0.7, color="tab:green",
                        alpha=0.4, transform=ax_bot.get_xaxis_transform(),
                        label="constructive (in phase)")
    ax_bot.fill_between(year, 0, 1, where=c < -0.7, color="tab:red",
                        alpha=0.4, transform=ax_bot.get_xaxis_transform(),
                        label="destructive (cancel)")
    if actual is not None:
        ax2 = ax_bot.twinx()
        ax2.plot(year, actual, color="tab:blue", linewidth=0.5, alpha=0.8,
                 label=actual_label)
        ax2.set_ylabel(actual_label, color="tab:blue")
        ax2.legend(loc="upper right", fontsize=7)
    ax_bot.set_ylim(0, 1)
    ax_bot.set_yticks([])
    ax_bot.set_xlabel("year")
    ax_bot.legend(loc="upper left", fontsize=7)


def demo_pdo():
    idx = "pdo"
    year, dt, model, obs, forcing = load_columns(idx)
    w = compute_winding(idx)
    m, amp, phase = (np.asarray(w["m"]), np.asarray(w["amp"]),
                      np.asarray(w["phase"]))

    m_grid = np.arange(0.0, M_MAX + DM / 2, DM)
    t0_grid = np.arange(year[0] + SIGMA / 2, year[-1] - SIGMA / 2 + 1e-9,
                         T0_STEP)
    log_snr, edge_mask = compute_snr_grid(year, standardize(model), forcing,
                                           m_grid, t0_grid, SIGMA)
    ridges, strength, persistence = detect_ridges(m_grid, log_snr, edge_mask)

    print(f"[{idx}] detected {len(ridges)} unbroken ridges:")
    for r in ridges:
        print(f"   M={m_grid[r]:.3f}  strength={strength[r]:+.2f}  "
              f"persistence={persistence[r]:.2f}")

    (gm1, ga1, gp1), (gm2, ga2, gp2) = closest_pair(m, amp, phase)
    print(f"[{idx}] globally closest fitted pair (for reference): "
          f"M1={gm1:.4f} (amp={ga1:.4f}), M2={gm2:.4f} (amp={ga2:.4f}), "
          f"dM={abs(gm1 - gm2):.4f}")

    # The specific pair already discussed in WINDING_NARRATIVE.md's "Deeper
    # dive" section (PDO's dominant mode vs. the shared M(NM) harmonic it
    # nearly overlaps) — used here as the worked example, rather than
    # whichever gap is globally smallest.
    (m1, a1, p1), (m2, a2, p2) = pair_near(m, amp, phase, 0.449, 0.415)
    dm = abs(m1 - m2)
    print(f"[{idx}] worked-example pair (from WINDING_NARRATIVE.md): "
          f"M1={m1:.4f} (amp={a1:.4f}), M2={m2:.4f} (amp={a2:.4f}), "
          f"dM={dm:.4f}")

    nearest_to_m1 = ridges[np.argmin(np.abs(m_grid[ridges] - abs(m1)))]
    nearest_to_m2 = ridges[np.argmin(np.abs(m_grid[ridges] - abs(m2)))]
    if nearest_to_m1 == nearest_to_m2:
        print(f"[{idx}] the ridge detector assigns BOTH M1={m1:.4f} and "
              f"M2={m2:.4f} to the same single ridge at "
              f"M={m_grid[nearest_to_m1]:.3f} — true dM={dm:.4f} is below "
              f"its resolving power. {NOT_RESOLVED_MSG}")
    else:
        print(f"[{idx}] the ridge detector resolves M1={m1:.4f} and "
              f"M2={m2:.4f} as two separate ridges at "
              f"M={m_grid[nearest_to_m1]:.3f} and "
              f"M={m_grid[nearest_to_m2]:.3f} (dM={dm:.4f}).")

    fig, axes = plt.subplots(3, 1, figsize=(11, 11),
                              gridspec_kw={"height_ratios": [1, 1.3, 0.5]})
    plot_ridge_profile(axes[0], m_grid, strength, persistence, ridges,
                        highlight=[m1, m2],
                        title=f"{idx}: ridge strength vs winding number "
                              f"(green dashed = the close pair M1={m1:.3f}, "
                              f"M2={m2:.3f})")
    plot_beat(axes[1], axes[2], year, forcing, m1, a1, p1, m2, a2, p2,
              actual=standardize(model), actual_label="standardized Model")
    axes[1].set_title(f"{idx}: exact beat envelope of the near-degenerate "
                       f"pair (dM={dm:.4f})", fontweight="bold")
    fig.suptitle(f"Winding ridge + beat analysis — {idx} (worked example)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = ROOT / idx / "winding_beat_demo.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}\n")


NOT_RESOLVED_MSG = "NOT resolved (merges into one apparent peak)"


def cross_index_nino34():
    src_idx, ref_idx = "nino34", "nino4"
    year_ref, dt_ref, model_ref, obs_ref, forcing_ref = load_columns(ref_idx)

    raw = np.loadtxt(ROOT / src_idx / f"{src_idx}.dat")
    year_src, val_src = raw[:, 0], raw[:, 1]
    val_on_ref_grid = np.interp(year_ref, year_src, val_src)
    x = standardize(val_on_ref_grid)

    m_grid = np.arange(0.0, M_MAX + DM / 2, DM)
    t0_grid = np.arange(year_ref[0] + SIGMA / 2,
                         year_ref[-1] - SIGMA / 2 + 1e-9, T0_STEP)
    log_snr, edge_mask = compute_snr_grid(year_ref, x, forcing_ref, m_grid,
                                           t0_grid, SIGMA)
    ridges, strength, persistence = detect_ridges(m_grid, log_snr, edge_mask)

    w_ref = compute_winding(ref_idx)
    ref_m = np.abs(np.asarray(w_ref["m"]))

    print(f"[{src_idx} raw data, tested against {ref_idx}'s own Forcing "
          f"manifold] detected {len(ridges)} unbroken ridges:")
    matched = []
    for r in ridges:
        m_here = m_grid[r]
        nearest_ref = ref_m[np.argmin(np.abs(ref_m - m_here))]
        d = abs(m_here - nearest_ref)
        tag = f"  <-- matches {ref_idx} M={nearest_ref:.3f} (dM={d:.3f})" \
            if d < 0.05 else ""
        if tag:
            matched.append(m_here)
        print(f"   M={m_here:.3f}  strength={strength[r]:+.2f}  "
              f"persistence={persistence[r]:.2f}{tag}")
    print(f"[{src_idx}] {len(matched)}/{len(ridges)} detected ridges land "
          f"within 0.05 of one of {ref_idx}'s own fitted winding numbers, "
          f"using only {src_idx}'s raw data and {ref_idx}'s manifold — "
          f"no separate fit of {src_idx} was used for this test.")

    # closest fitted pair from the REFERENCE index, checked against nino34's
    # own raw data via reconstruction using the reference's amp/phase.
    m_ref, amp_ref, phase_ref = (np.asarray(w_ref["m"]),
                                  np.asarray(w_ref["amp"]),
                                  np.asarray(w_ref["phase"]))
    (m1, a1, p1), (m2, a2, p2) = closest_pair(m_ref, amp_ref, phase_ref)
    dm = abs(m1 - m2)
    print(f"[{ref_idx}] closest fitted pair reused on {src_idx}: "
          f"M1={m1:.4f}, M2={m2:.4f}, dM={dm:.4f}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 11),
                              gridspec_kw={"height_ratios": [1, 1.3, 0.5]})
    plot_ridge_profile(axes[0], m_grid, strength, persistence, ridges,
                        highlight=ref_m,
                        title=f"{src_idx} raw data winding against "
                              f"{ref_idx}'s Forcing (green dashed = "
                              f"{ref_idx}'s own fitted windings)")
    plot_beat(axes[1], axes[2], year_ref, forcing_ref, m1, a1, p1, m2, a2, p2,
              actual=x, actual_label=f"standardized {src_idx} raw data")
    axes[1].set_title(f"{ref_idx}'s near-degenerate pair (dM={dm:.4f}), "
                       f"reconstructed and compared against {src_idx}'s own "
                       f"raw data", fontweight="bold")
    fig.suptitle(f"Cross-index winding test — {src_idx} data against "
                 f"{ref_idx}'s manifold (no separate fit of {src_idx})",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = ROOT / src_idx / f"winding_beat_vs_{ref_idx}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}\n")


if __name__ == "__main__":
    demo_pdo()
    cross_index_nino34()
