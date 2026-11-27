#!/usr/bin/env python3
"""Benchmark GEM.Tidal_Poly multi-year free-drift tiers (Z10_N5, Z10_D6_N5) on:
  1. ATL4 (NAO + Baltic + AMO + TNA, 1960-1970 held out, exact R5 repro_0693 setup)
  2. IP3  (NINO4 + PDO + IODE, 1960-1970 held out, exact indopacific_family setup)
comparing directly against 52-line TABLE (ALL) and sub-annual-only POLY tiers.
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
IP3_DIR = ROOT / "Oct2026_generator/legs/indopacific_family"
OUT_ROOT = ROOT / "Oct2026_generator/legs/compact_tidal_poly/n5_bench"
OUT_ROOT.mkdir(parents=True, exist_ok=True)


def score_csv(csv_path: Path, ho_lo=1960.0, ho_hi=1970.0):
    arr = np.loadtxt(csv_path, delimiter=",")
    t, y, m, f = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    valid = y != 0.0
    ho = valid & (t >= ho_lo) & (t <= ho_hi)
    tr = valid & (~ho)
    cc_tr = float(np.corrcoef(y[tr], m[tr])[0, 1])
    cc_ho = float(np.corrcoef(y[ho], m[ho])[0, 1])
    return t, y, m, f, cc_tr, cc_ho


def run_suite_case(suite_name: str, src_dir: Path, leg_names: list, leg_files: list,
                   tag: str, tides_mode: str, poly_tier: str = None,
                   iters: int = 80, lam: float = 0.30):
    wdir = OUT_ROOT / f"{suite_name}_{tag}"
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True)

    resp_txt = (src_dir / "lt.exe.resp").read_text()
    if tides_mode == "POLY":
        resp_txt += f"\nTIDES                          POLY\nPOLY_TIER                      {poly_tier}\n"
    (wdir / "lt.exe.resp").write_text(resp_txt)
    n_legs = len(leg_names)
    for k in range(1, n_legs + 1):
        shutil.copy2(src_dir / f"lt.exe.resp.{k}", wdir / f"lt.exe.resp.{k}")
        shutil.copy2(src_dir / f"lt.exe.p.{k}.pre_refine", wdir / f"lt.exe.p.{k}")
    shutil.copy2(src_dir / "lt.exe.p.pre_refine", wdir / "lt.exe.p")
    for fn in leg_files:
        (wdir / fn).symlink_to((src_dir / fn).resolve())
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
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=480
    )
    elapsed = time.time() - t0
    (wdir / "refine.log").write_text(p_ref.stdout)
    if p_ref.returncode != 0:
        raise RuntimeError(f"{suite_name}_{tag} REFINE failed:\n{p_ref.stdout[-1200:]}")

    m_cal = re.search(r"dLOD CC=([0-9.]+)|([0-9.]+)=CC", p_ref.stdout)
    dlod_cal = float(m_cal.group(1) or m_cal.group(2)) if m_cal else float("nan")
    m_sum = re.search(r"REFINE summary: mean ([0-9.]+) -> ([0-9.]+); dLOD ([0-9.]+) -> ([0-9.]+)", p_ref.stdout)
    dlod_fit = float(m_sum.group(4)) if m_sum else float("nan")

    # TEST_ONLY (Post-Bessel)
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
        raise RuntimeError(f"{suite_name}_{tag} TEST_ONLY failed:\n{p_test.stdout[-1000:]}")

    # Pre-Bessel dump
    pre_dir = wdir / "pre_bessel"
    pre_dir.mkdir()
    for fn in ["lt.exe.resp", "lt.exe.p"] + [f"lt.exe.resp.{k}" for k in range(1, n_legs + 1)] + [f"lt.exe.p.{k}" for k in range(1, n_legs + 1)]:
        shutil.copy2(wdir / fn, pre_dir / fn)
    for fn in leg_files:
        (pre_dir / fn).symlink_to((src_dir / fn).resolve())
    pre_bin = pre_dir / "lt.exe"
    pre_bin.symlink_to(ENSO_OPT)
    env_pre = env_test.copy()
    env_pre["MANIFOLD_PRE"] = "TRUE"
    subprocess.run([str(pre_bin)], cwd=pre_dir, env=env_pre,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)

    legs = {}
    for k, name in enumerate(leg_names, start=1):
        _, _, _, _, tr, ho = score_csv(wdir / f"lte_results.{k}.csv")
        legs[name] = {"tr": tr, "ho": ho}

    # Leg 2 is Baltic in ATL4 and PDO in IP3!
    t2, y2, m2, f2_post, _, _ = score_csv(wdir / "lte_results.2.csv")
    _, _, _, f2_pre, _, _ = score_csv(pre_dir / "lte_results.2.csv")
    return {
        "tag": tag,
        "elapsed_s": round(elapsed, 1),
        "dlod_cal": dlod_cal,
        "dlod_fit": dlod_fit,
        "legs": legs,
        "mean_tr": float(np.mean([legs[n]["tr"] for n in leg_names])),
        "mean_ho": float(np.mean([legs[n]["ho"] for n in leg_names])),
        "t2": t2, "y2": y2, "m2": m2, "f2_post": f2_post, "f2_pre": f2_pre,
    }


def main():
    # ── 1. ATL4: Z10_N5 and Z10_D6_N5 vs R5 (ALL) ──────────────────────────
    atl_names = ["NAO", "Baltic", "AMO", "TNA"]
    atl_files = ["nao.dat", "baltic_3pt.dat", "amo.dat", "tna.dat"]
    print("=" * 115)
    print("SUITE 1: ATL4 (NAO + Baltic + AMO + TNA, 1960-1970 held out) — Multi-year N5 free of sub-annual dLOD prior")
    print("=" * 115)
    print(f"{'Tier':<11} {'Time':>6} {'dLOD0':>7} {'dLOD1':>7} | {'NAO tr/ho':>13} {'Baltic tr/ho':>14} {'AMO tr/ho':>13} {'TNA tr/ho':>13} | {'MeanHO':>7} {'PreCC':>7} {'PostCC':>7}")
    print("-" * 115)

    t_b, y_b, m_r5, f_r5_post, _, _ = score_csv(R5_DIR / "lte_results.2.csv")
    r5_pre_csv = ROOT / "Oct2026_generator/legs/compact_tidal_poly/atl4_poly_bench/r5_pre_bessel/lte_results.2.csv"
    _, _, _, f_r5_pre, _, _ = score_csv(r5_pre_csv)
    r5_legs = {}
    for k, n in enumerate(atl_names, start=1):
        _, _, _, _, tr, ho = score_csv(R5_DIR / f"lte_results.{k}.csv")
        r5_legs[n] = {"tr": tr, "ho": ho}
    print(f"{'R5(ALL)':<11} {'127.0s':>6} {'0.9717':>7} {'0.9544':>7} | "
          f"{r5_legs['NAO']['tr']:+.3f}/{r5_legs['NAO']['ho']:+.3f} "
          f"{r5_legs['Baltic']['tr']:+.3f}/{r5_legs['Baltic']['ho']:+.3f}  "
          f"{r5_legs['AMO']['tr']:+.3f}/{r5_legs['AMO']['ho']:+.3f} "
          f"{r5_legs['TNA']['tr']:+.3f}/{r5_legs['TNA']['ho']:+.3f} | "
          f"{np.mean([r5_legs[n]['ho'] for n in atl_names]):+.3f} {'1.0000':>7} {'1.0000':>7}", flush=True)

    atl_res = []
    for tier in ["Z10_N5", "Z10_D6_N5"]:
        r = run_suite_case("atl4", R5_DIR, atl_names, atl_files, tier, "POLY", tier, iters=80, lam=0.30)
        r["pre_cc"] = float(np.corrcoef(r["f2_pre"], f_r5_pre)[0, 1])
        r["post_cc"] = float(np.corrcoef(r["f2_post"], f_r5_post)[0, 1])
        atl_res.append(r)
        L = r["legs"]
        print(f"{tier:<11} {r['elapsed_s']:5.1f}s {r['dlod_cal']:7.4f} {r['dlod_fit']:7.4f} | "
              f"{L['NAO']['tr']:+.3f}/{L['NAO']['ho']:+.3f} "
              f"{L['Baltic']['tr']:+.3f}/{L['Baltic']['ho']:+.3f}  "
              f"{L['AMO']['tr']:+.3f}/{L['AMO']['ho']:+.3f} "
              f"{L['TNA']['tr']:+.3f}/{L['TNA']['ho']:+.3f} | "
              f"{r['mean_ho']:+.3f} {r['pre_cc']:7.4f} {r['post_cc']:7.4f}", flush=True)

    # ── 2. IP3 (NINO4 + PDO + IODE, 1960-1970 held out) ────────────────────
    ip_names = ["NINO4", "PDO", "IODE"]
    ip_files = ["nino4.dat", "pdo.dat", "iode.dat"]
    print("\n" + "=" * 115)
    print("SUITE 2: IP3 (NINO4 + PDO + IODE, 1960-1970 held out) — PDO evaluation across TABLE vs POLY tiers")
    print("=" * 115)
    print(f"{'Tier':<11} {'Time':>6} {'dLOD0':>7} {'dLOD1':>7} | {'NINO4 tr/ho':>13} {'PDO tr/ho':>14} {'IODE tr/ho':>13} | {'MeanHO':>7} {'PreCC':>7} {'PostCC':>7}")
    print("-" * 115)

    ip_all = run_suite_case("ip3", IP3_DIR, ip_names, ip_files, "ALL", "TABLE", None, iters=80, lam=0.30)
    ip_all["pre_cc"] = 1.0
    ip_all["post_cc"] = 1.0
    L = ip_all["legs"]
    print(f"{'IP3(ALL)':<11} {ip_all['elapsed_s']:5.1f}s {ip_all['dlod_cal']:7.4f} {ip_all['dlod_fit']:7.4f} | "
          f"{L['NINO4']['tr']:+.3f}/{L['NINO4']['ho']:+.3f} "
          f"{L['PDO']['tr']:+.3f}/{L['PDO']['ho']:+.3f}  "
          f"{L['IODE']['tr']:+.3f}/{L['IODE']['ho']:+.3f} | "
          f"{ip_all['mean_ho']:+.3f} {'1.0000':>7} {'1.0000':>7}", flush=True)

    ip_res = [ip_all]
    for tier in ["Z10", "Z10_N5", "Z10_D6", "Z10_D6_N5"]:
        r = run_suite_case("ip3", IP3_DIR, ip_names, ip_files, tier, "POLY", tier, iters=80, lam=0.30)
        r["pre_cc"] = float(np.corrcoef(r["f2_pre"], ip_all["f2_pre"])[0, 1])
        r["post_cc"] = float(np.corrcoef(r["f2_post"], ip_all["f2_post"])[0, 1])
        ip_res.append(r)
        L = r["legs"]
        print(f"{tier:<11} {r['elapsed_s']:5.1f}s {r['dlod_cal']:7.4f} {r['dlod_fit']:7.4f} | "
              f"{L['NINO4']['tr']:+.3f}/{L['NINO4']['ho']:+.3f} "
              f"{L['PDO']['tr']:+.3f}/{L['PDO']['ho']:+.3f}  "
              f"{L['IODE']['tr']:+.3f}/{L['IODE']['ho']:+.3f} | "
              f"{r['mean_ho']:+.3f} {r['pre_cc']:7.4f} {r['post_cc']:7.4f}", flush=True)

    # Save JSON
    def strip_arrs(r):
        return {k: v for k, v in r.items() if k not in ("t2", "y2", "m2", "f2_post", "f2_pre")}
    (OUT_ROOT / "n5_atl4_pdo_results.json").write_text(json.dumps({
        "atl4_r5": r5_legs,
        "atl4_poly": [strip_arrs(r) for r in atl_res],
        "ip3": [strip_arrs(r) for r in ip_res],
    }, indent=2))

    # Plot 1: ATL4 Baltic comparison (R5 vs Z10_N5 vs Z10_D6_N5)
    fig, axes = plt.subplots(5, 1, figsize=(15, 12.5), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1.3, 1.3, 1.0, 1.0]})
    axes[0].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    axes[0].plot(t_b, y_b, color="#444444", lw=0.8, alpha=0.8, label="Baltic 3pt obs")
    axes[0].plot(t_b, m_r5, color="#d62728", lw=1.1,
                 label=f"R5 (52-line TABLE): Baltic tr={r5_legs['Baltic']['tr']:+.3f}, ho={r5_legs['Baltic']['ho']:+.3f} | ATL4 mean_ho=+0.226")
    axes[0].set_ylabel("Baltic (sd)")
    axes[0].set_title("ATL4 (1960-1970 held out): Free Multi-Year Nodal/Perigee (N5) + Sub-Annual dLOD Prior vs R5 (ALL)", fontsize=12, fontweight="bold")
    axes[0].legend(loc="upper left", fontsize=8.5)
    axes[0].grid(alpha=0.25)

    cols = {"Z10_N5": "#1f77b4", "Z10_D6_N5": "#2ca02c"}
    for i, r in enumerate(atl_res, start=1):
        axes[i].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
        axes[i].plot(t_b, y_b, color="#444444", lw=0.8, alpha=0.8, label="Baltic 3pt obs")
        axes[i].plot(t_b, r["m2"], color=cols[r["tag"]], lw=1.1,
                     label=f"POLY {r['tag']} ({r['elapsed_s']:.0f}s): Baltic tr={r['legs']['Baltic']['tr']:+.3f}, ho={r['legs']['Baltic']['ho']:+.3f} | ATL4 mean_ho={r['mean_ho']:+.3f}, dLOD={r['dlod_fit']:.3f}")
        axes[i].set_ylabel("Baltic (sd)")
        axes[i].legend(loc="upper left", fontsize=8.5)
        axes[i].grid(alpha=0.25)

    axes[3].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    axes[3].plot(t_b, f_r5_post, color="#d62728", lw=1.2, label="R5 Post-Bessel F")
    for r in atl_res:
        axes[3].plot(t_b, r["f2_post"], color=cols[r["tag"]], lw=0.9, ls="--",
                     label=f"POLY {r['tag']} (CC vs R5={r['post_cc']:.4f})")
    axes[3].set_ylabel("Post-Bessel F")
    axes[3].legend(loc="upper left", ncol=3, fontsize=8.5)
    axes[3].grid(alpha=0.25)

    axes[4].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    axes[4].plot(t_b, f_r5_pre, color="#d62728", lw=1.2, label="R5 Pre-Bessel F")
    for r in atl_res:
        axes[4].plot(t_b, r["f2_pre"], color=cols[r["tag"]], lw=0.9, ls="--",
                     label=f"POLY {r['tag']} (CC vs R5={r['pre_cc']:.4f})")
    axes[4].set_ylabel("Pre-Bessel F")
    axes[4].set_xlabel("Year")
    axes[4].legend(loc="upper left", ncol=3, fontsize=8.5)
    axes[4].grid(alpha=0.25)
    plt.xlim(1880, 2025)
    plt.tight_layout()
    plt.savefig(OUT_ROOT / "atl4_n5_comparison.png", dpi=150)

    # Plot 2: IP3 PDO comparison (ALL vs Z10 vs Z10_N5 vs Z10_D6_N5)
    fig, axes = plt.subplots(6, 1, figsize=(15, 14.5), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1.3, 1.3, 1.3, 1.0, 1.0]})
    ip_map = {r["tag"]: r for r in ip_res}
    t_p, y_p = ip_all["t2"], ip_all["y2"]
    ip_show = [("ALL", "#d62728"), ("Z10", "#ff7f0e"), ("Z10_N5", "#1f77b4"), ("Z10_D6_N5", "#2ca02c")]
    for idx, (tag, col) in enumerate(ip_show):
        r = ip_map[tag]
        ax = axes[idx]
        ax.axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
        ax.plot(t_p, y_p, color="#444444", lw=0.8, alpha=0.8, label="PDO obs")
        ax.plot(t_p, r["m2"], color=col, lw=1.1,
                label=f"{'IP3 (52-line TABLE)' if tag == 'ALL' else 'POLY ' + tag} ({r['elapsed_s']:.0f}s): PDO tr={r['legs']['PDO']['tr']:+.3f}, ho={r['legs']['PDO']['ho']:+.3f} | IP3 mean_ho={r['mean_ho']:+.3f}, dLOD={r['dlod_fit']:.3f}")
        ax.set_ylabel("PDO")
        if idx == 0:
            ax.set_title("IP3 (NINO4 + PDO + IODE, 1960-1970 held out): GEM.Tidal_Poly (Free Multi-Year N5) vs 52-Line TABLE", fontsize=12, fontweight="bold")
        ax.legend(loc="upper left", fontsize=8.5)
        ax.grid(alpha=0.25)

    axes[4].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    axes[4].plot(t_p, ip_all["f2_post"], color="#d62728", lw=1.2, label="IP3(ALL) Post-Bessel F")
    for tag, col in ip_show[1:]:
        r = ip_map[tag]
        axes[4].plot(t_p, r["f2_post"], color=col, lw=0.9, ls="--",
                     label=f"POLY {tag} (CC vs ALL={r['post_cc']:.4f})")
    axes[4].set_ylabel("Post-Bessel F")
    axes[4].legend(loc="upper left", ncol=2, fontsize=8.5)
    axes[4].grid(alpha=0.25)

    axes[5].axvspan(1960, 1970, color="#ffe8cc", alpha=0.7)
    axes[5].plot(t_p, ip_all["f2_pre"], color="#d62728", lw=1.2, label="IP3(ALL) Pre-Bessel F")
    for tag, col in ip_show[1:]:
        r = ip_map[tag]
        axes[5].plot(t_p, r["f2_pre"], color=col, lw=0.9, ls="--",
                     label=f"POLY {tag} (CC vs ALL={r['pre_cc']:.4f})")
    axes[5].set_ylabel("Pre-Bessel F")
    axes[5].set_xlabel("Year")
    axes[5].legend(loc="upper left", ncol=2, fontsize=8.5)
    axes[5].grid(alpha=0.25)
    plt.xlim(1880, 2025)
    plt.tight_layout()
    plt.savefig(OUT_ROOT / "pdo_n5_comparison.png", dpi=150)
    print(f"\nSaved figures -> {OUT_ROOT}/atl4_n5_comparison.png and pdo_n5_comparison.png")


if __name__ == "__main__":
    main()
