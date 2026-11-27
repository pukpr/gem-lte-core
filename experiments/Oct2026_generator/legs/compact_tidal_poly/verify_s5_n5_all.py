#!/usr/bin/env python3
"""Verify Z10_S5_N5 & Z10_S5_D6_N5 (with light N5 ridge = 0.05) on:
  A. Exact 30-year rolling-window Baltic non-stationarity diagnostic (before vs after Z10)
  B. ATL4 (NAO + Baltic + AMO + TNA, 1960-1970 held out)
  C. IP3  (NINO4 + PDO + IODE, 1960-1970 held out)
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
IP3_DIR = ROOT / "Oct2026_generator/legs/indopacific_family"
BALTIC_DAT = (R5_DIR / "baltic_3pt.dat").resolve()
Z10_WIN_JSON = ROOT / "Oct2026_generator/legs/compact_tidal_poly/windowed_30yr_baltic/windowed_30yr_results.json"
OUT_ROOT = ROOT / "Oct2026_generator/legs/compact_tidal_poly/s5_n5_verify"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

TERM_LABELS_10 = [
    "Mf(2L)", "Mf_nod", "2F", "Mm(Mp)", "MSm(ev)",
    "Mt", "Mt_nod", "MSf(2D)", "Ssa(2Ls)", "2D+Mp(9.61d)"
]
WINDOWS = [
    (1885.0, 1915.0), (1900.0, 1930.0), (1915.0, 1945.0), (1930.0, 1960.0),
    (1945.0, 1975.0), (1960.0, 1990.0), (1975.0, 2005.0), (1990.0, 2020.0),
]


def score_csv(csv_path: Path, ho_lo=1960.0, ho_hi=1970.0):
    arr = np.loadtxt(csv_path, delimiter=",")
    t, y, m, f = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    valid = y != 0.0
    ho = valid & (t >= ho_lo) & (t <= ho_hi)
    tr = valid & (~ho)
    return t, y, m, f, float(np.corrcoef(y[tr], m[tr])[0, 1]), float(np.corrcoef(y[ho], m[ho])[0, 1])


def run_baltic_window(tag: str, t_start: float, t_end: float, tier: str = "Z10_S5_N5"):
    wdir = OUT_ROOT / f"bwin_{tier}_{tag}"
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True)
    resp = f"""-FSTART                        0.1
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
POLY_TIER                      {tier}
"""
    (wdir / "lt.exe.resp").write_text(resp)
    p_base = json.loads((R5_DIR / "lt.exe.p.pre_refine").read_text())
    p_leg2 = json.loads((R5_DIR / "lt.exe.p.2.pre_refine").read_text())
    for k, v in p_leg2.items():
        if k != "lpap":
            p_base[k] = v
    (wdir / "lt.exe.p").write_text(json.dumps(p_base, indent=2))
    (wdir / "baltic_3pt.dat").symlink_to(BALTIC_DAT)
    lt_bin = wdir / "lt.exe"
    lt_bin.symlink_to(ENSO_OPT)

    env = os.environ.copy()
    env.update({
        "DLOD_DAT": str(DLOD_DAT), "DLOD_REF": "FALSE", "EXCLUDE": "FALSE",
        "TRAIN_START": str(t_start), "TRAIN_END": str(t_end),
        "REFINE": "TRUE", "REFINE_LAMBDA": "0.30", "REFINE_HOMOTOPY": "TRUE",
        "REFINE_ROUNDS": "4", "REFINE_ITERS": "60", "REFINE_APPLY": "TRUE",
        "NUMBER_OF_PROCESSORS": "4",
    })
    t0 = time.time()
    p_ref = subprocess.run([str(lt_bin)], cwd=wdir, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=240)
    elapsed = time.time() - t0
    if p_ref.returncode != 0:
        raise RuntimeError(f"bwin_{tag} failed:\n{p_ref.stdout[-800:]}")
    m_sum = re.search(r"REFINE summary: mean ([0-9.]+) -> ([0-9.]+); dLOD ([0-9.]+) -> ([0-9.]+)", p_ref.stdout)
    dlod1 = float(m_sum.group(4)) if m_sum else float("nan")

    env_t = env.copy()
    env_t.update({"REFINE": "FALSE", "DLOD_REF": "TRUE", "TEST_ONLY": "TRUE"})
    subprocess.run([str(lt_bin)], cwd=wdir, env=env_t,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    arr = np.loadtxt(wdir / "lte_results.csv", delimiter=",")
    t, y, m = arr[:, 0], arr[:, 1], arr[:, 2]
    valid = y != 0.0
    in_w = valid & (t >= t_start) & (t <= t_end)
    adj = valid & (((t >= t_start - 15.0) & (t < t_start)) | ((t > t_end) & (t <= t_end + 15.0)))
    cc_in = float(np.corrcoef(y[in_w], m[in_w])[0, 1])
    cc_adj = float(np.corrcoef(y[adj], m[adj])[0, 1]) if np.sum(adj) > 10 else float("nan")
    p_fit = json.loads((wdir / "lt.exe.p").read_text())
    lpap = p_fit["lpap"]
    phasors = []
    for k in range(10):
        per, amp, ph = lpap[k]
        phasors.append((float(per), float(amp), float(ph), float(amp) * math.cos(float(ph)), float(amp) * math.sin(float(ph))))
    return {"tag": tag, "t_mid": 0.5 * (t_start + t_end), "elapsed_s": round(elapsed, 1),
            "dlod1": dlod1, "cc_in": cc_in, "cc_adj": cc_adj,
            "cw_mag": math.hypot(float(p_fit.get("cw1", 0.0)), float(p_fit.get("cw2", 0.0))),
            "phasors": phasors}


def run_suite(suite_name: str, src_dir: Path, leg_names: list, leg_files: list, tier: str):
    wdir = OUT_ROOT / f"{suite_name}_{tier}"
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True)
    resp_txt = (src_dir / "lt.exe.resp").read_text()
    resp_txt += f"\nTIDES                          POLY\nPOLY_TIER                      {tier}\n"
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
        "DLOD_DAT": str(DLOD_DAT), "DLOD_REF": "FALSE", "EXCLUDE": "TRUE",
        "TRAIN_START": "1960.0", "TRAIN_END": "1970.0",
        "REFINE": "TRUE", "REFINE_LAMBDA": "0.30", "REFINE_HOMOTOPY": "TRUE",
        "REFINE_ROUNDS": "4", "REFINE_ITERS": "80", "REFINE_APPLY": "TRUE",
        "NUMBER_OF_PROCESSORS": "4",
    })
    t0 = time.time()
    p_ref = subprocess.run([str(lt_bin)], cwd=wdir, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=480)
    elapsed = time.time() - t0
    (wdir / "refine.log").write_text(p_ref.stdout)
    if p_ref.returncode != 0:
        raise RuntimeError(f"{suite_name}_{tier} failed:\n{p_ref.stdout[-1000:]}")
    m_sum = re.search(r"REFINE summary: mean ([0-9.]+) -> ([0-9.]+); dLOD ([0-9.]+) -> ([0-9.]+)", p_ref.stdout)
    dlod1 = float(m_sum.group(4)) if m_sum else float("nan")

    env_t = env.copy()
    env_t.update({"REFINE": "FALSE", "DLOD_REF": "TRUE", "TEST_ONLY": "TRUE"})
    subprocess.run([str(lt_bin)], cwd=wdir, env=env_t,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    legs = {}
    for k, name in enumerate(leg_names, start=1):
        _, _, _, _, tr, ho = score_csv(wdir / f"lte_results.{k}.csv")
        legs[name] = {"tr": tr, "ho": ho}
    return {
        "tier": tier, "elapsed_s": round(elapsed, 1), "dlod1": dlod1,
        "legs": legs, "mean_ho": float(np.mean([legs[n]["ho"] for n in leg_names])),
        "n_pos": sum(1 for n in leg_names if legs[n]["ho"] > 0),
    }


def main():
    z10_old = json.loads(Z10_WIN_JSON.read_text())

    # ── A. 30-Year Rolling-Window Baltic Diagnostic: Z10_S5_N5 vs Z10 ────────
    print("=" * 105)
    print("PART A: 30-YEAR WINDOWED BALTIC DIAGNOSTIC — POLY Z10_S5_N5 (20t, light N5 ridge=0.05) vs Z10 (10t)")
    print("=" * 105)
    glob_new = run_baltic_window("global_1885_2020", 1885.0, 2020.0, "Z10_S5_N5")
    wins_new = [run_baltic_window(f"{int(a)}_{int(b)}", a, b, "Z10_S5_N5") for a, b in WINDOWS]

    old_adj = [w["cc_adj"] for w in z10_old["windows"]]
    new_adj = [w["cc_adj"] for w in wins_new]
    print(f"Global 135-yr CC : Z10 = {z10_old['global']['cc_in']:+.4f}  ->  Z10_S5_N5 = {glob_new['cc_in']:+.4f}  (dLOD1={glob_new['dlod1']:.4f})")
    print(f"Mean 30-yr In-Win: Z10 = {np.mean([w['cc_in'] for w in z10_old['windows']]):+.4f}  ->  Z10_S5_N5 = {np.mean([w['cc_in'] for w in wins_new]):+.4f}")
    print(f"Mean Adj ±15y OOS: Z10 = {np.mean(old_adj):+.4f}  ->  Z10_S5_N5 = {np.mean(new_adj):+.4f}")
    print("\nPhasor Non-Stationarity Index (NSI) & Phase Drift (Before Z10 -> After Z10_S5_N5):")
    print(f"{'Term':<14} | {'Z10 NSI':>8} {'New NSI':>8} | {'Z10 Δφ':>8} {'New Δφ':>8}")
    print("-" * 58)
    new_nsi = []
    for k in range(10):
        zs = np.array([complex(w["phasors"][k][3], w["phasors"][k][4]) for w in wins_new])
        std_dz = float(np.sqrt(np.mean(np.abs(zs - np.mean(zs)) ** 2)))
        nsi = std_dz / max(float(np.mean(np.abs(zs))), 1e-9)
        ph_deg = np.rad2deg(np.unwrap(np.angle(zs)))
        ph_span = float(np.max(ph_deg) - np.min(ph_deg))
        new_nsi.append({"term": TERM_LABELS_10[k], "nsi": nsi, "ph_span": ph_span})
        o = z10_old["phasor_nsi"][k]
        print(f"{TERM_LABELS_10[k]:<14} | {o['nsi']:8.3f} {nsi:8.3f} | {o['ph_span_deg']:7.1f}° {ph_span:7.1f}°", flush=True)

    # ── B. ATL4 & IP3 (PDO) Verification ─────────────────────────────────────
    print("\n" + "=" * 105)
    print("PART B: ATL4 (NAO + Baltic + AMO + TNA, 1960-1970 held out) — Z10_S5_N5 & Z10_S5_D6_N5")
    print("=" * 105)
    atl_names = ["NAO", "Baltic", "AMO", "TNA"]
    atl_files = ["nao.dat", "baltic_3pt.dat", "amo.dat", "tna.dat"]
    atl_out = []
    for tier in ["Z10_S5_N5", "Z10_S5_D6_N5"]:
        r = run_suite("atl4", R5_DIR, atl_names, atl_files, tier)
        atl_out.append(r)
        L = r["legs"]
        print(f"{tier:<14} {r['elapsed_s']:5.1f}s dLOD1={r['dlod1']:.4f} | "
              f"NAO {L['NAO']['tr']:+.3f}/{L['NAO']['ho']:+.3f}  "
              f"Baltic {L['Baltic']['tr']:+.3f}/{L['Baltic']['ho']:+.3f}  "
              f"AMO {L['AMO']['tr']:+.3f}/{L['AMO']['ho']:+.3f}  "
              f"TNA {L['TNA']['tr']:+.3f}/{L['TNA']['ho']:+.3f} | MeanHO={r['mean_ho']:+.3f} ({r['n_pos']}/4)", flush=True)

    print("\n" + "=" * 105)
    print("PART C: IP3 (NINO4 + PDO + IODE, 1960-1970 held out) — Z10_S5_N5 & Z10_S5_D6_N5")
    print("=" * 105)
    ip_names = ["NINO4", "PDO", "IODE"]
    ip_files = ["nino4.dat", "pdo.dat", "iode.dat"]
    ip_out = []
    for tier in ["Z10_S5_N5", "Z10_S5_D6_N5"]:
        r = run_suite("ip3", IP3_DIR, ip_names, ip_files, tier)
        ip_out.append(r)
        L = r["legs"]
        print(f"{tier:<14} {r['elapsed_s']:5.1f}s dLOD1={r['dlod1']:.4f} | "
              f"NINO4 {L['NINO4']['tr']:+.3f}/{L['NINO4']['ho']:+.3f}  "
              f"PDO {L['PDO']['tr']:+.3f}/{L['PDO']['ho']:+.3f}  "
              f"IODE {L['IODE']['tr']:+.3f}/{L['IODE']['ho']:+.3f} | MeanHO={r['mean_ho']:+.3f} ({r['n_pos']}/3)", flush=True)

    # Save JSON
    (OUT_ROOT / "s5_n5_verify_results.json").write_text(json.dumps({
        "baltic_30yr_global_old": z10_old["global"]["cc_in"],
        "baltic_30yr_global_new": glob_new["cc_in"],
        "baltic_30yr_adj_old": old_adj,
        "baltic_30yr_adj_new": new_adj,
        "nsi_old": z10_old["phasor_nsi"],
        "nsi_new": new_nsi,
        "atl4": atl_out,
        "ip3": ip_out,
    }, indent=2))

    # Plot Before vs After NSI + Adj15y OOS
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.2))
    x = np.arange(10)
    w = 0.36
    ax1.bar(x - w/2, [o["nsi"] for o in z10_old["phasor_nsi"]], w, color="#d62728", alpha=0.85, label="Before: POLY Z10 (10t)")
    ax1.bar(x + w/2, [n["nsi"] for n in new_nsi], w, color="#2ca02c", alpha=0.85, label="After: POLY Z10_S5_N5 (20t)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(TERM_LABELS_10, rotation=30, ha="right", fontsize=8.5)
    ax1.set_ylabel("Non-Stationarity Index (NSI)")
    ax1.set_title("Phasor Non-Stationarity Across 30-Yr Baltic Windows", fontweight="bold")
    ax1.legend(fontsize=8.5)
    ax1.grid(axis="y", alpha=0.25)

    t_mids = [w["t_mid"] for w in wins_new]
    ax2.plot(t_mids, old_adj, "s--", color="#d62728", lw=1.8, label=f"Z10 Adj ±15y OOS (mean={np.mean(old_adj):+.3f})")
    ax2.plot(t_mids, new_adj, "o-", color="#2ca02c", lw=2.0, label=f"Z10_S5_N5 Adj ±15y OOS (mean={np.mean(new_adj):+.3f})")
    ax2.axhline(0.0, color="#666666", lw=0.8)
    ax2.set_xlabel("30-Yr Window Midpoint Year")
    ax2.set_ylabel("Adjacent ±15-Yr OOS CC")
    ax2.set_title(f"Global 135y CC: {z10_old['global']['cc_in']:+.3f} → {glob_new['cc_in']:+.3f}", fontweight="bold")
    ax2.legend(fontsize=8.5)
    ax2.grid(alpha=0.25)

    plt.tight_layout()
    fig_path = OUT_ROOT / "baltic_30yr_z10_vs_z10_s5_n5.png"
    plt.savefig(fig_path, dpi=150)
    print(f"\nSaved comparison figure -> {fig_path}")


if __name__ == "__main__":
    main()
