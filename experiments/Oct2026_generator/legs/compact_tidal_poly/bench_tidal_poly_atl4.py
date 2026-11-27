#!/usr/bin/env python3
"""Benchmark GEM.Tidal_Poly (TIDES=POLY) tiers on the exact 4-leg ATL4 repro_0693
benchmark (NAO + baltic_3pt_na + AMO + TNA, 1960-1970 held out, NH=10, MAXH=45,
lam=0.30, cold dLOD + homotopy) against R5_atl4_cold_H10_l03 (52-line TABLE).
"""
import json
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
OUT_ROOT = ROOT / "Oct2026_generator/legs/compact_tidal_poly/atl4_poly_bench"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

TIERS = ["Z9", "Z10", "Z15", "Z10_D6", "Z15_D6"]
LEG_NAMES = ["NAO", "Baltic", "AMO", "TNA"]
LEG_FILES = ["nao.dat", "baltic_3pt.dat", "amo.dat", "tna.dat"]


def score_csv(csv_path: Path, ho_lo=1960.0, ho_hi=1970.0):
    arr = np.loadtxt(csv_path, delimiter=",")
    t, y, m, f = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    valid = y != 0.0
    ho = valid & (t >= ho_lo) & (t <= ho_hi)
    tr = valid & (~ho)
    cc_tr = float(np.corrcoef(y[tr], m[tr])[0, 1])
    cc_ho = float(np.corrcoef(y[ho], m[ho])[0, 1])
    return t, y, m, f, cc_tr, cc_ho


def max_winding_amp(w_json: Path) -> float:
    if not w_json.exists():
        return float("nan")
    d = json.loads(w_json.read_text())
    amps = [abs(w.get("amp", 0.0)) for w in d.get("windings", [])]
    return max(amps) if amps else 0.0


def run_tier(tier: str):
    wdir = OUT_ROOT / f"poly_{tier}"
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True)

    # Copy resp files, adding TIDES POLY and POLY_TIER <tier>
    resp_txt = (R5_DIR / "lt.exe.resp").read_text()
    resp_txt += f"\nTIDES                          POLY\nPOLY_TIER                      {tier}\n"
    (wdir / "lt.exe.resp").write_text(resp_txt)
    for k in range(1, 5):
        shutil.copy2(R5_DIR / f"lt.exe.resp.{k}", wdir / f"lt.exe.resp.{k}")
        shutil.copy2(R5_DIR / f"lt.exe.p.{k}.pre_refine", wdir / f"lt.exe.p.{k}")
    shutil.copy2(R5_DIR / "lt.exe.p.pre_refine", wdir / "lt.exe.p")
    for fn in LEG_FILES:
        src = (R5_DIR / fn).resolve()
        (wdir / fn).symlink_to(src)
    lt_bin = wdir / "lt.exe"
    lt_bin.symlink_to(ENSO_OPT)

    env = os.environ.copy()
    env.update({
        "DLOD_DAT": str(DLOD_DAT),
        "DLOD_REF": "FALSE",
        "EXCLUDE": "TRUE",
        "TRAIN_START": "1960.0",
        "TRAIN_END": "1970.0",
        "REFINE": "TRUE",
        "REFINE_LAMBDA": "0.30",
        "REFINE_HOMOTOPY": "TRUE",
        "REFINE_ROUNDS": "4",
        "REFINE_ITERS": "80",
        "REFINE_APPLY": "TRUE",
        "NUMBER_OF_PROCESSORS": "4",
    })

    t0 = time.time()
    p_ref = subprocess.run(
        [str(lt_bin)], cwd=wdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=360
    )
    elapsed = time.time() - t0
    (wdir / "refine.log").write_text(p_ref.stdout)
    if p_ref.returncode != 0:
        raise RuntimeError(f"{tier} REFINE failed (rc={p_ref.returncode}):\n{p_ref.stdout[-1500:]}")

    # Parse dLOD calibration CC, analytic J spot-check, and final refine dLOD CC
    m_cal = re.search(r"TIDES=POLY tier=\S+.*dLOD CC=([0-9.]+)", p_ref.stdout)
    dlod_cal = float(m_cal.group(1)) if m_cal else float("nan")
    m_j = re.search(r"max diff = ([0-9.]+)", p_ref.stdout)
    j_diff = float(m_j.group(1)) if m_j else float("nan")
    m_sum = re.search(r"REFINE summary: mean ([0-9.]+) -> ([0-9.]+); dLOD ([0-9.]+) -> ([0-9.]+)", p_ref.stdout)
    dlod_fit = float(m_sum.group(4)) if m_sum else float("nan")

    # Run TEST_ONLY=TRUE (Post-Bessel) to verify .p round-trip and dump lte_results.*.csv
    env_test = env.copy()
    env_test["REFINE"] = "FALSE"
    env_test["DLOD_REF"] = "TRUE"
    env_test["TEST_ONLY"] = "TRUE"
    p_test = subprocess.run(
        [str(lt_bin)], cwd=wdir, env=env_test,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120
    )
    (wdir / "score.log").write_text(p_test.stdout)
    if p_test.returncode != 0:
        raise RuntimeError(f"{tier} TEST_ONLY failed:\n{p_test.stdout[-1000:]}")

    # Also dump Pre-Bessel manifold for Baltic (leg 2)
    env_pre = env_test.copy()
    env_pre["MANIFOLD_PRE"] = "TRUE"
    pre_dir = wdir / "pre_bessel"
    pre_dir.mkdir()
    for fn in ["lt.exe.resp", "lt.exe.p"] + [f"lt.exe.resp.{k}" for k in range(1, 5)] + [f"lt.exe.p.{k}" for k in range(1, 5)]:
        shutil.copy2(wdir / fn, pre_dir / fn)
    for fn in LEG_FILES:
        (pre_dir / fn).symlink_to((R5_DIR / fn).resolve())
    pre_bin = pre_dir / "lt.exe"
    pre_bin.symlink_to(ENSO_OPT)
    subprocess.run([str(pre_bin)], cwd=pre_dir, env=env_pre,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)

    legs = {}
    for k, name in enumerate(LEG_NAMES, start=1):
        _, _, _, _, tr, ho = score_csv(wdir / f"lte_results.{k}.csv")
        legs[name] = {
            "tr": tr,
            "ho": ho,
            "max_amp": max_winding_amp(wdir / f"lt.exe.windings.{k}.json"),
        }

    # Pre-Bessel and Post-Bessel Baltic manifold correlation vs R5
    _, _, _, f_post, _, _ = score_csv(wdir / "lte_results.2.csv")
    _, _, _, f_pre, _, _ = score_csv(pre_dir / "lte_results.2.csv")
    return {
        "tier": tier,
        "elapsed_s": round(elapsed, 1),
        "dlod_cal": dlod_cal,
        "dlod_fit": dlod_fit,
        "j_diff": j_diff,
        "legs": legs,
        "mean_tr": float(np.mean([legs[n]["tr"] for n in LEG_NAMES])),
        "mean_ho": float(np.mean([legs[n]["ho"] for n in LEG_NAMES])),
        "n_pos_ho": sum(1 for n in LEG_NAMES if legs[n]["ho"] > 0),
        "f_post": f_post,
        "f_pre": f_pre,
    }


def main():
    # Load R5 (52-line TABLE) reference scores & Baltic Pre/Post manifolds
    r5_legs = {}
    for k, name in enumerate(LEG_NAMES, start=1):
        _, _, _, _, tr, ho = score_csv(R5_DIR / f"lte_results.{k}.csv")
        r5_legs[name] = {"tr": tr, "ho": ho, "max_amp": max_winding_amp(R5_DIR / f"lt.exe.windings.{k}.json")}
    t_b, y_b, m_r5, f_r5_post, tr_r5, ho_r5 = score_csv(R5_DIR / "lte_results.2.csv")

    # Get R5 Pre-Bessel manifold
    r5_pre_dir = OUT_ROOT / "r5_pre_bessel"
    if not (r5_pre_dir / "lte_results.2.csv").exists():
        r5_pre_dir.mkdir(parents=True, exist_ok=True)
        for fn in ["lt.exe.resp", "lt.exe.p"] + [f"lt.exe.resp.{k}" for k in range(1, 5)] + [f"lt.exe.p.{k}" for k in range(1, 5)]:
            shutil.copy2(R5_DIR / fn, r5_pre_dir / fn)
        for fn in LEG_FILES:
            p = r5_pre_dir / fn
            if not p.exists():
                p.symlink_to((R5_DIR / fn).resolve())
        r5_bin = r5_pre_dir / "lt.exe"
        if not r5_bin.exists():
            r5_bin.symlink_to(ENSO_OPT)
        env_r5 = os.environ.copy()
        env_r5.update({
            "DLOD_DAT": str(DLOD_DAT), "DLOD_REF": "TRUE", "EXCLUDE": "TRUE",
            "TRAIN_START": "1960.0", "TRAIN_END": "1970.0", "TEST_ONLY": "TRUE",
            "MANIFOLD_PRE": "TRUE",
        })
        subprocess.run([str(r5_bin)], cwd=r5_pre_dir, env=env_r5,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    _, _, _, f_r5_pre, _, _ = score_csv(r5_pre_dir / "lte_results.2.csv")

    results = []
    print(f"{'Tier':<8} {'Time':>6} {'dLOD0':>7} {'dLOD1':>7} | {'NAO tr/ho':>13} {'Baltic tr/ho':>14} {'AMO tr/ho':>13} {'TNA tr/ho':>13} | {'MeanHO':>7} {'PreCC':>7}")
    print("-" * 116)
    print(f"{'R5(ALL)':<8} {'127.0s':>6} {'0.9717':>7} {'0.9544':>7} | "
          f"{r5_legs['NAO']['tr']:+.3f}/{r5_legs['NAO']['ho']:+.3f} "
          f"{r5_legs['Baltic']['tr']:+.3f}/{r5_legs['Baltic']['ho']:+.3f}  "
          f"{r5_legs['AMO']['tr']:+.3f}/{r5_legs['AMO']['ho']:+.3f} "
          f"{r5_legs['TNA']['tr']:+.3f}/{r5_legs['TNA']['ho']:+.3f} | "
          f"{np.mean([r5_legs[n]['ho'] for n in LEG_NAMES]):+.3f} {'1.0000':>7}")

    for tier in TIERS:
        r = run_tier(tier)
        r["pre_cc_vs_r5"] = float(np.corrcoef(r["f_pre"], f_r5_pre)[0, 1])
        r["post_cc_vs_r5"] = float(np.corrcoef(r["f_post"], f_r5_post)[0, 1])
        results.append(r)
        L = r["legs"]
        print(f"{tier:<8} {r['elapsed_s']:5.1f}s {r['dlod_cal']:7.4f} {r['dlod_fit']:7.4f} | "
              f"{L['NAO']['tr']:+.3f}/{L['NAO']['ho']:+.3f} "
              f"{L['Baltic']['tr']:+.3f}/{L['Baltic']['ho']:+.3f}  "
              f"{L['AMO']['tr']:+.3f}/{L['AMO']['ho']:+.3f} "
              f"{L['TNA']['tr']:+.3f}/{L['TNA']['ho']:+.3f} | "
              f"{r['mean_ho']:+.3f} {r['pre_cc_vs_r5']:7.4f}", flush=True)

    # Save JSON summary (excluding large arrays)
    out_json = {
        "r5_all": {
            "legs": r5_legs,
            "mean_tr": float(np.mean([r5_legs[n]["tr"] for n in LEG_NAMES])),
            "mean_ho": float(np.mean([r5_legs[n]["ho"] for n in LEG_NAMES])),
            "elapsed_s": 127.0,
            "dlod_cal": 0.9717,
            "dlod_fit": 0.9544,
        },
        "poly_tiers": [{k: v for k, v in r.items() if k not in ("f_post", "f_pre")} for r in results],
    }
    (OUT_ROOT / "atl4_poly_results.json").write_text(json.dumps(out_json, indent=2))

    # Plot 6-panel visual comparison matching baltic_0710_atl4_cold_H10.png
    fig, axes = plt.subplots(6, 1, figsize=(15, 14.5), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1.3, 1.3, 1.3, 1.0, 1.0]})
    show_tiers = ["Z9", "Z10", "Z15_D6"]
    res_map = {r["tier"]: r for r in results}

    # Panel 0: R5 (52-line ALL)
    ax = axes[0]
    ax.axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    ax.plot(t_b, y_b, color="#444444", lw=0.8, alpha=0.8, label="Baltic 3pt obs")
    ax.plot(t_b, m_r5, color="#d62728", lw=1.1, label=f"R5 (52-line TABLE, 104p): tr={tr_r5:+.3f}, ho={ho_r5:+.3f}, dLOD=0.954")
    ax.set_ylabel("Baltic (sd)")
    ax.set_title("ATL4 (NAO + Baltic + AMO + TNA, 1960-1970 held out): GEM.Tidal_Poly (Meeus-free) vs R5 (52-line TABLE)", fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    ax.grid(alpha=0.25)

    colors = {"Z9": "#1f77b4", "Z10": "#2ca02c", "Z15": "#9467bd", "Z10_D6": "#ff7f0e", "Z15_D6": "#8c564b"}
    for idx, tier in enumerate(show_tiers, start=1):
        r = res_map[tier]
        _, _, m_t, _, tr_t, ho_t = score_csv(OUT_ROOT / f"poly_{tier}" / "lte_results.2.csv")
        ax = axes[idx]
        ax.axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
        ax.plot(t_b, y_b, color="#444444", lw=0.8, alpha=0.8, label="Baltic 3pt obs")
        ax.plot(t_b, m_t, color=colors[tier], lw=1.1,
                label=f"TIDES=POLY {tier} ({r['elapsed_s']:.0f}s): Baltic tr={tr_t:+.3f}, ho={ho_t:+.3f} | ATL4 mean_ho={r['mean_ho']:+.3f}, dLOD={r['dlod_fit']:.3f}")
        ax.set_ylabel("Baltic (sd)")
        ax.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
        ax.grid(alpha=0.25)

    # Panel 4: Post-Bessel manifold overlay (R5 vs POLY tiers)
    ax = axes[4]
    ax.axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    ax.plot(t_b, f_r5_post, color="#d62728", lw=1.2, alpha=0.85, label="R5 (52-line TABLE) Post-Bessel F")
    for tier in show_tiers:
        r = res_map[tier]
        ax.plot(t_b, r["f_post"], color=colors[tier], lw=0.9, ls="--", alpha=0.85,
                label=f"POLY {tier} (CC vs R5={r['post_cc_vs_r5']:.4f})")
    ax.set_ylabel("Post-Bessel F")
    ax.legend(loc="upper left", ncol=2, fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25)

    # Panel 5: Pre-Bessel manifold overlay (R5 vs POLY tiers)
    ax = axes[5]
    ax.axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    ax.plot(t_b, f_r5_pre, color="#d62728", lw=1.2, alpha=0.85, label="R5 (52-line TABLE) Pre-Bessel F")
    for tier in show_tiers:
        r = res_map[tier]
        ax.plot(t_b, r["f_pre"], color=colors[tier], lw=0.9, ls="--", alpha=0.85,
                label=f"POLY {tier} (CC vs R5={r['pre_cc_vs_r5']:.4f})")
    ax.set_ylabel("Pre-Bessel F")
    ax.set_xlabel("Year")
    ax.legend(loc="upper left", ncol=2, fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25)

    plt.xlim(1880, 2025)
    plt.tight_layout()
    fig_path = OUT_ROOT / "baltic_0710_atl4_poly_vs_r5.png"
    plt.savefig(fig_path, dpi=150)
    print(f"\nSaved figure -> {fig_path}")


if __name__ == "__main__":
    main()
