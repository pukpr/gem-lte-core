#!/usr/bin/env python3
"""30-Year Rolling-Window Non-Stationarity Diagnostic on baltic_3pt_na.dat.

Fits TIDES=POLY (Z10) on 8 overlapping 30-year windows (step 15 yr, 1885-2020)
plus 1 full-record global fit, then diagnoses:
  1. Which Z10 phasors (Amp, Phase) drift / rotate across 30-yr windows (beat -> missing sideband!).
  2. How each local 30-yr Pre-Bessel manifold F_w(t) deviates from F_global(t):
     Delta_F(t) = F_local(t) - F_global(t) stitched across 1885-2020.
  3. Exact spectral peaks of Delta_F(t) and top missing Doodson/orbital constituents
     needed to stitch the 30-year local solutions into one stationary global model.
"""
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/home/paul/gemini/gem-lte-core/experiments")
ENSO_OPT = ROOT / "Feb2026/enso_opt"
DLOD_DAT = ROOT / "Feb2026/dlod3.dat"
R5_DIR = ROOT / "Oct2026_generator/legs/baltic_3pt_overnight/repro_0693/R5_atl4_cold_H10_l03"
BALTIC_DAT = (R5_DIR / "baltic_3pt.dat").resolve()
OUT_ROOT = ROOT / "Oct2026_generator/legs/compact_tidal_poly/windowed_30yr_baltic"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

TIER = "Z10"
TERM_LABELS = [
    "Mf(2L)", "Mf_nod", "2F", "Mm(Mp)", "MSm(ev)",
    "Mt", "Mt_nod", "MSf(2D)", "Ssa(2Ls)", "2D+Mp(9.61d)"
]

WINDOWS = [
    (1885.0, 1915.0),
    (1900.0, 1930.0),
    (1915.0, 1945.0),
    (1930.0, 1960.0),
    (1945.0, 1975.0),
    (1960.0, 1990.0),
    (1975.0, 2005.0),
    (1990.0, 2020.0),
]


def setup_single_leg_dir(wdir: Path, tier: str):
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True)
    resp = """-FSTART                        0.1
FLIP                           0.0015
SPREAD_MIN                     0.000001
SPREAD_MAX                     0.45
MAXH                           45
NH                             "1 1 1 1 1 1 1 1 1 1"
NM                             3
YEAR                           0.00405305
DECAY                          48.0
DTW                            3
IDATE                          1880.0
INIT_DATE                      1880.0
STRICT_IDATE                   TRUE
THRESHOLD_ACTION               CATCHUP
PART                           FALSE
FILTER                         0.0
JERK                           0.0
CLIMATE_INDEX                  baltic_3pt.dat
TIDES                          POLY
POLY_TIER                      """ + tier + "\n"
    (wdir / "lt.exe.resp").write_text(resp)

    # Build single-leg seed lt.exe.p from R5's pre_refine Baltic leg (leg 2) + shared tidal seed
    p_base = json.loads((R5_DIR / "lt.exe.p.pre_refine").read_text())
    p_leg2 = json.loads((R5_DIR / "lt.exe.p.2.pre_refine").read_text())
    for k, v in p_leg2.items():
        if k != "lpap":
            p_base[k] = v
    (wdir / "lt.exe.p").write_text(json.dumps(p_base, indent=2))
    (wdir / "baltic_3pt.dat").symlink_to(BALTIC_DAT)
    (wdir / "lt.exe").symlink_to(ENSO_OPT)


def run_window(tag: str, t_start: float, t_end: float, iters: int = 60, lam: float = 0.30):
    wdir = OUT_ROOT / f"win_{tag}"
    setup_single_leg_dir(wdir, TIER)
    lt_bin = wdir / "lt.exe"

    env = os.environ.copy()
    env.update({
        "DLOD_DAT": str(DLOD_DAT),
        "DLOD_REF": "FALSE",
        "EXCLUDE": "FALSE",
        "TRAIN_START": str(t_start),
        "TRAIN_END": str(t_end),
        "REFINE": "TRUE",
        "REFINE_LAMBDA": str(lam),
        "REFINE_HOMOTOPY": "TRUE",
        "REFINE_ROUNDS": "4",
        "REFINE_ITERS": str(iters),
        "REFINE_APPLY": "TRUE",
        "NUMBER_OF_PROCESSORS": "4",
    })

    t0 = time.time()
    p_ref = subprocess.run(
        [str(lt_bin)], cwd=wdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=240
    )
    elapsed = time.time() - t0
    (wdir / "refine.log").write_text(p_ref.stdout)
    if p_ref.returncode != 0:
        raise RuntimeError(f"win_{tag} REFINE failed:\n{p_ref.stdout[-1000:]}")

    m_sum = re.search(r"REFINE summary: mean ([0-9.]+) -> ([0-9.]+); dLOD ([0-9.]+) -> ([0-9.]+)", p_ref.stdout)
    dlod1 = float(m_sum.group(4)) if m_sum else float("nan")

    # TEST_ONLY (Post-Bessel) on full record
    env_t = env.copy()
    env_t["REFINE"] = "FALSE"
    env_t["DLOD_REF"] = "TRUE"
    env_t["TEST_ONLY"] = "TRUE"
    subprocess.run([str(lt_bin)], cwd=wdir, env=env_t,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)

    # Pre-Bessel dump
    pre_dir = wdir / "pre_bessel"
    pre_dir.mkdir()
    for fn in ("lt.exe.resp", "lt.exe.p"):
        shutil.copy2(wdir / fn, pre_dir / fn)
    (pre_dir / "baltic_3pt.dat").symlink_to(BALTIC_DAT)
    (pre_dir / "lt.exe").symlink_to(ENSO_OPT)
    env_pre = env_t.copy()
    env_pre["MANIFOLD_PRE"] = "TRUE"
    subprocess.run([str(pre_dir / "lt.exe")], cwd=pre_dir, env=env_pre,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)

    arr_post = np.loadtxt(wdir / "lte_results.csv", delimiter=",")
    arr_pre = np.loadtxt(pre_dir / "lte_results.csv", delimiter=",")
    t, y, m, f_post = arr_post[:, 0], arr_post[:, 1], arr_post[:, 2], arr_post[:, 3]
    f_pre = arr_pre[:, 3]
    valid = y != 0.0
    in_w = valid & (t >= t_start) & (t <= t_end)
    out_w = valid & (~in_w)
    cc_in = float(np.corrcoef(y[in_w], m[in_w])[0, 1])
    cc_out = float(np.corrcoef(y[out_w], m[out_w])[0, 1]) if np.sum(out_w) > 10 else float("nan")

    # Adjacent +/-15 yr OOS CC
    adj = valid & (((t >= t_start - 15.0) & (t < t_start)) | ((t > t_end) & (t <= t_end + 15.0)))
    cc_adj = float(np.corrcoef(y[adj], m[adj])[0, 1]) if np.sum(adj) > 10 else float("nan")

    # Extract fitted lpap phasors (first 10 triplets) and cw1/cw2/Bessel params
    p_fit = json.loads((wdir / "lt.exe.p").read_text())
    lpap = p_fit["lpap"]
    phasors = []
    for k in range(10):
        per, amp, ph = lpap[k]
        phasors.append((float(per), float(amp), float(ph), float(amp) * math.cos(float(ph)), float(amp) * math.sin(float(ph))))

    return {
        "tag": tag,
        "t_start": t_start,
        "t_end": t_end,
        "t_mid": 0.5 * (t_start + t_end),
        "elapsed_s": round(elapsed, 1),
        "dlod1": dlod1,
        "cc_in": cc_in,
        "cc_adj": cc_adj,
        "cc_out": cc_out,
        "cw1": float(p_fit.get("cw1", 0.0)),
        "cw2": float(p_fit.get("cw2", 0.0)),
        "impA": float(p_fit.get("impA", 0.0)),
        "impB": float(p_fit.get("impB", 0.0)),
        "besK": float(p_fit.get("besK", 0.0)),
        "phasors": phasors,
        "t": t, "y": y, "m": m, "f_post": f_post, "f_pre": f_pre,
    }


def main():
    print("=" * 105)
    print(f"30-YEAR ROLLING-WINDOW NON-STATIONARITY DIAGNOSTIC: baltic_3pt_na.dat (TIDES=POLY {TIER})")
    print("=" * 105)
    print(f"{'Window':<14} {'Time':>6} {'dLOD1':>7} | {'In-Win CC':>10} {'Adj15y CC':>10} {'Rest CC':>9} | {'cw1':>7} {'cw2':>7} {'|cw|':>6}")
    print("-" * 105)

    glob = run_window("global_1885_2020", 1885.0, 2020.0)
    cw_mag = math.hypot(glob["cw1"], glob["cw2"])
    print(f"{'GLOBAL 135y':<14} {glob['elapsed_s']:5.1f}s {glob['dlod1']:7.4f} | "
          f"{glob['cc_in']:+10.4f} {'   n/a':>10} {'   n/a':>9} | "
          f"{glob['cw1']:+7.2f} {glob['cw2']:+7.2f} {cw_mag:6.2f}", flush=True)

    wins = []
    for a, b in WINDOWS:
        r = run_window(f"{int(a)}_{int(b)}", a, b)
        wins.append(r)
        cm = math.hypot(r["cw1"], r["cw2"])
        print(f"{int(a)}-{int(b):<9} {r['elapsed_s']:5.1f}s {r['dlod1']:7.4f} | "
              f"{r['cc_in']:+10.4f} {r['cc_adj']:+10.4f} {r['cc_out']:+9.4f} | "
              f"{r['cw1']:+7.2f} {r['cw2']:+7.2f} {cm:6.2f}", flush=True)

    # ── 1. Non-Stationarity Index (NSI) of each Z10 phasor across the 8 windows ──
    print("\n" + "=" * 105)
    print("DIAGNOSTIC 1: PHASOR NON-STATIONARITY ACROSS 30-YEAR WINDOWS (high NSI => missing beat/sideband!)")
    print("=" * 105)
    print(f"{'Term':<14} {'Period(d)':>10} | {'Glob Amp':>9} {'Mean Amp':>9} {'Std|dz|':>9} {'NSI':>7} | {'Phase Drift (deg)':>18}")
    print("-" * 105)
    nsi_list = []
    t_mids = np.array([w["t_mid"] for w in wins])
    for k in range(10):
        z_glob = complex(glob["phasors"][k][3], glob["phasors"][k][4])
        zs = np.array([complex(w["phasors"][k][3], w["phasors"][k][4]) for w in wins])
        z_mean = np.mean(zs)
        std_dz = float(np.sqrt(np.mean(np.abs(zs - z_mean) ** 2)))
        mean_amp = float(np.mean(np.abs(zs)))
        nsi = std_dz / max(mean_amp, 1e-9)
        ph_deg = np.rad2deg(np.unwrap(np.angle(zs)))
        ph_span = float(np.max(ph_deg) - np.min(ph_deg))
        nsi_list.append({
            "term": TERM_LABELS[k],
            "period_d": glob["phasors"][k][0],
            "glob_amp": abs(z_glob),
            "mean_amp": mean_amp,
            "std_dz": std_dz,
            "nsi": nsi,
            "ph_span_deg": ph_span,
            "amps": [float(abs(z)) for z in zs],
            "phs_deg": [float(x) for x in ph_deg],
        })
        print(f"{TERM_LABELS[k]:<14} {glob['phasors'][k][0]:10.4f} | "
              f"{abs(z_glob):9.5f} {mean_amp:9.5f} {std_dz:9.5f} {nsi:7.3f} | {ph_span:14.1f} deg")

    # ── 2. Stitch local Pre-Bessel F_w(t) using raised-cosine weights & compute Delta_F(t) ──
    t = glob["t"]
    f_glob = glob["f_pre"]
    num = np.zeros_like(t)
    den = np.zeros_like(t)
    m_stitched = np.zeros_like(t)
    for w in wins:
        a, b = w["t_start"], w["t_end"]
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        u = (t - mid) / half
        mask = np.abs(u) <= 1.0
        wt = np.zeros_like(t)
        wt[mask] = 0.5 * (1.0 + np.cos(np.pi * u[mask]))
        # Align constant offset of local f_pre to global f_pre inside window
        off = float(np.mean(w["f_pre"][mask] - f_glob[mask]))
        num += wt * (w["f_pre"] - off)
        m_stitched += wt * w["m"]
        den += wt
    ok = (den > 1e-6) & (glob["y"] != 0.0)
    f_stitched = np.where(ok, num / np.maximum(den, 1e-12), f_glob)
    m_stitched = np.where(ok, m_stitched / np.maximum(den, 1e-12), glob["m"])
    delta_f = f_stitched[ok] - f_glob[ok]
    t_ok = t[ok]
    cc_stitched_model = float(np.corrcoef(glob["y"][ok], m_stitched[ok])[0, 1])

    # FFT of Delta_F(t) (monthly series -> frequencies in cycles/year!)
    dt = float(np.median(np.diff(t_ok)))
    df_zero = delta_f - np.mean(delta_f)
    win_hann = np.hanning(len(df_zero))
    spec = np.abs(np.fft.rfft(df_zero * win_hann))
    freqs = np.fft.rfftfreq(len(df_zero), d=dt)  # cycles/year
    order = np.argsort(spec)[::-1]
    top_peaks = []
    for idx in order:
        f_cpy = float(freqs[idx])
        if f_cpy < 0.015:
            continue
        # Avoid duplicates within 0.03 cpy
        if any(abs(f_cpy - p["f_cpy"]) < 0.03 for p in top_peaks):
            continue
        top_peaks.append({
            "f_cpy": f_cpy,
            "period_yr": 1.0 / f_cpy,
            "amp": float(spec[idx] / (0.5 * np.sum(win_hann))),
        })
        if len(top_peaks) >= 10:
            break

    print("\n" + "=" * 105)
    print(f"DIAGNOSTIC 2: STITCHED 30-YR LOCAL MODEL CC = {cc_stitched_model:+.4f} vs GLOBAL = {glob['cc_in']:+.4f}")
    print("Top Spectral Peaks of Stitched Manifold Residual Delta_F(t) = F_local(t) - F_global(t):")
    print("=" * 105)
    for i, pk in enumerate(top_peaks, start=1):
        print(f"  #{i:2d}: f = {pk['f_cpy']:.4f} cpy  (Period = {pk['period_yr']:6.2f} yr)  rel_amp = {pk['amp']:.5f}")

    # Save JSON summary
    out_json = {
        "global": {k: v for k, v in glob.items() if k not in ("t", "y", "m", "f_post", "f_pre")},
        "windows": [{k: v for k, v in w.items() if k not in ("t", "y", "m", "f_post", "f_pre")} for w in wins],
        "phasor_nsi": nsi_list,
        "stitched_model_cc": cc_stitched_model,
        "delta_f_top_peaks": top_peaks,
    }
    (OUT_ROOT / "windowed_30yr_results.json").write_text(json.dumps(out_json, indent=2))

    # ── 3. 4-Panel Diagnostic Figure ─────────────────────────────────────────
    fig, axes = plt.subplots(4, 1, figsize=(15, 13.5), gridspec_kw={"height_ratios": [1.1, 1.1, 1.0, 1.1]})

    # Panel 1: In-window CC vs Adjacent 15-yr OOS CC vs Global CC
    ax = axes[0]
    ax.plot(t_mids, [w["cc_in"] for w in wins], "o-", color="#1f77b4", lw=2, label="Local 30-yr In-Window CC")
    ax.plot(t_mids, [w["cc_adj"] for w in wins], "s--", color="#ff7f0e", lw=1.8, label="Adjacent ±15-yr OOS CC")
    ax.axhline(glob["cc_in"], color="#d62728", ls=":", lw=1.6, label=f"Global 135-yr POLY Z10 CC = {glob['cc_in']:+.3f}")
    ax.axhline(cc_stitched_model, color="#2ca02c", ls="-.", lw=1.4, label=f"Stitched 30-yr Local Model CC = {cc_stitched_model:+.3f}")
    ax.set_ylabel("Correlation (CC)")
    ax.set_title("30-Year Rolling-Window Non-Stationarity Diagnostic on baltic_3pt_na.dat (TIDES=POLY Z10)", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", ncol=2, fontsize=8.5)
    ax.grid(alpha=0.25)

    # Panel 2: Non-stationarity of top-drifting Z10 phasors & cw1/cw2 magnitude
    ax = axes[1]
    top_nsi = sorted(range(10), key=lambda k: nsi_list[k]["nsi"], reverse=True)[:5]
    pal = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
    for c_i, k in enumerate(top_nsi):
        amps = np.array(nsi_list[k]["amps"])
        ax.plot(t_mids, amps / max(np.mean(amps), 1e-9), "o-", color=pal[c_i], lw=1.7,
                label=f"{TERM_LABELS[k]} (NSI={nsi_list[k]['nsi']:.2f}, Δφ={nsi_list[k]['ph_span_deg']:.0f}°)")
    cw_mags = np.array([math.hypot(w["cw1"], w["cw2"]) for w in wins])
    ax.plot(t_mids, cw_mags / max(np.mean(cw_mags), 1e-9), "k^--", lw=1.5, label=f"Chandler |cw| (mean={np.mean(cw_mags):.1f})")
    ax.set_ylabel("Normalized Amp")
    ax.set_xlabel("Window Midpoint Year")
    ax.legend(loc="upper right", ncol=3, fontsize=8.2)
    ax.grid(alpha=0.25)

    # Panel 3: Stitched Delta_F(t) = F_local(t) - F_global(t)
    ax = axes[2]
    ax.plot(t_ok, delta_f, color="#1f77b4", lw=1.2, label="Stitched ΔF(t) = F_local(t) − F_global(t)  [missing manifold drift]")
    ax.axhline(0.0, color="#666666", lw=0.8)
    ax.set_ylabel("ΔF (Pre-Bessel)")
    ax.set_xlabel("Year")
    ax.legend(loc="upper left", fontsize=8.5)
    ax.grid(alpha=0.25)

    # Panel 4: Spectrum of Delta_F(t) showing exact missing frequencies!
    ax = axes[3]
    mask_f = (freqs >= 0.015) & (freqs <= 2.5)
    ax.plot(freqs[mask_f], spec[mask_f] / (0.5 * np.sum(win_hann)), color="#d62728", lw=1.4)
    # Annotate known physical frequencies: 18.61y, 8.85y, 6.0y, 4.43y, Chandler (0.844 cpy), 1.0 cpy (Sa), 2.0 cpy (Ssa)
    ref_lines = [
        (1.0 / 18.613, "18.6y (N)", "#1f77b4"),
        (1.0 / 8.850, "8.85y (p)", "#2ca02c"),
        (1.0 / 6.000, "6.0y (p+N)", "#9467bd"),
        (0.8442, "CW/2F (1.18y)", "#ff7f0e"),
        (1.0, "1.0y (Sa)", "#8c564b"),
        (2.0, "0.5y (Ssa)", "#e377c2"),
    ]
    ymax = float(np.max(spec[mask_f] / (0.5 * np.sum(win_hann))))
    for f_ref, lbl, col in ref_lines:
        ax.axvline(f_ref, color=col, ls="--", lw=1.1, alpha=0.8)
        ax.text(f_ref + 0.015, ymax * 0.88, lbl, color=col, fontsize=8, rotation=90, va="top")
    ax.set_xlabel("Frequency (cycles/year)")
    ax.set_ylabel("ΔF Spectral Amp")
    ax.set_title("Spectrum of Missing Manifold Residual ΔF(t) = F_local(t) − F_global(t)", fontsize=10.5)
    ax.grid(alpha=0.25)

    plt.tight_layout()
    fig_path = OUT_ROOT / "baltic_30yr_nonstationarity_diag.png"
    plt.savefig(fig_path, dpi=150)
    print(f"\nSaved diagnostic figure -> {fig_path}")


if __name__ == "__main__":
    main()
