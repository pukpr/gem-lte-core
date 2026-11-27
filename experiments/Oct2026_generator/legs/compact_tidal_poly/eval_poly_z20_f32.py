#!/usr/bin/env python3
"""Evaluate Z20 (Z18 + idx 34 9.6137d) and F32 (F30 + idx 34 9.6137d) on ATL3 & IP2,
merge into abc_results.json, and plot complete 5-mask comparison (ALL, F32, Z20, F30, Z18).
"""
import concurrent.futures
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

LEGS_ROOT = Path("/home/paul/gemini/gem-lte-core/experiments/Oct2026_generator/legs")
ATL = LEGS_ROOT / "atlantic_family"
IND7 = LEGS_ROOT / "7indices"
HERE = LEGS_ROOT / "compact_tidal_poly"
WORK = HERE / "abc_eval"

DATA_MAP = {
    "nao":        ATL / "nao.dat",
    "baltic_3pt": ATL / "baltic_3pt_na.dat",
    "amo":        ATL / "amo.dat",
    "nino4":      IND7 / "nino4.dat",
    "pdo":        IND7 / "pdo.dat",
}

ATL_P = json.loads((ATL / "lt.exe.p").read_text())
H7 = ATL_P["harm"][:7]
H10 = [2, 3, 4, 6, 9, 13, 18, 24, 31, 41]
FOLDS_5 = [("1970_80", 1970, 1980), ("1980_90", 1980, 1990), ("1990_00", 1990, 2000),
           ("2000_10", 2000, 2010), ("2010_20", 2010, 2020)]

def corr(a, b):
    if len(a) < 4 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])

def run_job(job):
    family, legs, harm, maxh, mask, wname, h0, h1 = job
    slug = f"{family}_{mask}_{wname}"
    d = WORK / slug
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for i, nm in enumerate(legs, 1):
        os.symlink(DATA_MAP[nm], d / f"{nm}.dat")
        (d / f"lt.exe.resp.{i}").write_text(f"CLIMATE_INDEX                  {nm}.dat\n")
    nh_str = " ".join(["1"] * len(harm))
    (d / "lt.exe.resp").write_text(
        f"-FSTART                        0.1\nFLIP                           0.0015\n"
        f"SPREAD_MIN                     0.000001\nSPREAD_MAX                     0.45\n"
        f"MAXH                           {maxh}\n"
        f'NH                             "{nh_str}"\nNM                             3\n'
        f"YEAR                           0.00405305\nDECAY                          48.0\n"
        f"DTW                            3\nIDATE                          1880.0\n"
        f"INIT_DATE                      1880.0\nSTRICT_IDATE                   TRUE\n"
        f"THRESHOLD_ACTION               CATCHUP\nPART                           FALSE\n"
        f"FILTER                         0.0\nJERK                           0.0\n"
        f"LEGS                           {len(legs)}\nLEG_WINDINGS                   TRUE\n"
        f"LEG_P                          PARAMS\n"
    )
    p = json.loads(json.dumps(ATL_P))
    p["harm"] = [float(h) for h in harm]
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2))
    for i in range(1, len(legs) + 1):
        lp = json.loads((ATL / f"lt.exe.p.{min(i, 4)}").read_text())
        lp["harm"] = [float(h) for h in harm]
        (d / f"lt.exe.p.{i}").write_text(json.dumps(lp, indent=2))

    env = os.environ.copy()
    env.update({
        "DLOD_REF": "FALSE", "EXCLUDE": "true", "TRAIN_START": str(h0), "TRAIN_END": str(h1),
        "F9": "1", "IDATE": "1880.0", "LOCKA": "FALSE", "LOCKT": "FALSE", "METRIC": "CC",
        "TREND": "true", "TIDE_LEGACY": "FALSE", "NUMBER_OF_PROCESSORS": "4",
        "TIMEOUT": "0", "TEST_ONLY": "FALSE", "REFINE": "TRUE", "REFINE_ANCHOR_SEED": "FALSE",
        "REFINE_HOMOTOPY": "TRUE", "REFINE_LAM": "0.30", "REFINE_ROUNDS": "4",
        "REFINE_ITERS": "40", "REFINE_STAGE1": "0", "REFINE_TIDAL_MASK": mask,
    })
    t0 = time.time()
    r1 = subprocess.run("ulimit -s 65536 && ../../../lt.exe", shell=True, cwd=d, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=480)
    dt = time.time() - t0
    (d / "refine.log").write_text(r1.stdout)
    env["DLOD_REF"] = "TRUE"; env["REFINE"] = "FALSE"; env["TEST_ONLY"] = "TRUE"
    r2 = subprocess.run("ulimit -s 65536 && ../../../lt.exe", shell=True, cwd=d, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=240)
    (d / "score.log").write_text(r2.stdout)

    p_fin = json.loads((d / "lt.exe.p").read_text())
    n_active = sum(1 for row in p_fin["lpap"] if abs(row[1]) > 1e-12)
    m_dl = re.search(r"([+-]?\d+\.\d+):dLOD:", r2.stdout)
    dlod_cc = float(m_dl.group(1)) if m_dl else float("nan")
    scores, trs, hos = {}, [], []
    for i, nm in enumerate(legs, 1):
        a = np.loadtxt(d / f"lte_results.{i}.csv", delimiter=",")
        t, m, y = a[:, 0], a[:, 1], a[:, 2]
        ok = y != 0; ho = ok & (t >= h0) & (t < h1); tr = ok & ~ho
        r_tr, r_ho = corr(m[tr], y[tr]), corr(m[ho], y[ho])
        scores[nm] = {"tr": r_tr, "ho": r_ho}
        trs.append(r_tr); hos.append(r_ho)
    return {
        "slug": slug, "family": family, "mask": mask, "wname": wname, "window": [h0, h1],
        "time_s": dt, "n_active": n_active, "dlod_cc": dlod_cc,
        "mean_tr": float(np.mean(trs)), "mean_ho": float(np.mean(hos)), "scores": scores,
    }

res_path = WORK / "abc_results.json"
results = json.loads(res_path.read_text())

jobs = []
for mask in ["Z20", "F32"]:
    for wname, h0, h1 in [("canon_60_70", 1960, 1970)] + FOLDS_5:
        jobs.append(("ATL3", ["nao", "baltic_3pt", "amo"], H10, 45, mask, wname, h0, h1))
    for wname, h0, h1 in [("canon_00_05", 2000, 2005)] + FOLDS_5:
        jobs.append(("IP2", ["nino4", "pdo"], H7, 30, mask, wname, h0, h1))

print("=" * 140, flush=True)
print(f"ADDING idx=34 (9.6137 d, 127.3 yr secular-arch term): Z20 (10 lines, 20p) & F32 (16 lines, 32p) — {len(jobs)} runs", flush=True)
print("=" * 140, flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    futs = [ex.submit(run_job, j) for j in jobs]
    for fut in concurrent.futures.as_completed(futs):
        r = fut.result()
        results[r["slug"]] = r
        leg_s = " ".join(f"{k}:{v['tr']:+.3f}/{v['ho']:+.2f}" for k, v in r["scores"].items())
        print(
            f"{r['slug']:24s} | {r['time_s']:5.1f}s | L={r['n_active']:2d} dLOD={r['dlod_cc']:.4f} | "
            f"mean tr={r['mean_tr']:+.3f} ho={r['mean_ho']:+.3f} | {leg_s}",
            flush=True
        )

res_path.write_text(json.dumps(results, indent=2))

MASKS = ["ALL", "F32", "Z20", "F30", "Z18"]
print("\n" + "=" * 140, flush=True)
print("COMPLETE 5-MASK SUMMARY (ALL=104p, F32=32p, Z20=20p, F30=30p, Z18=18p): 5-FOLD ROLLING OOS CV (1970-2020)", flush=True)
print("=" * 140, flush=True)

summary = {"ATL3": {}, "IP2": {}}
for fam, canon_w, legs in [("ATL3", "canon_60_70", ["nao","baltic_3pt","amo"]),
                            ("IP2",  "canon_00_05", ["nino4","pdo"])]:
    print(f"\n--- {fam} ({' + '.join(legs)}) ---", flush=True)
    for mask in MASKS:
        c_row = results[f"{fam}_{mask}_{canon_w}"]
        cv_rows = [results[f"{fam}_{mask}_{wn}"] for wn, _, _ in FOLDS_5]
        t_mean = float(np.mean([r["time_s"] for r in cv_rows]))
        dl_mean = float(np.mean([r["dlod_cc"] for r in cv_rows]))
        tr_mean = float(np.mean([r["mean_tr"] for r in cv_rows]))
        ho_mean = float(np.mean([r["mean_ho"] for r in cv_rows]))
        summary[fam][mask] = {"L": c_row["n_active"], "t": t_mean, "dlod": dl_mean,
                              "tr": tr_mean, "ho": ho_mean, "c_tr": c_row["mean_tr"], "c_ho": c_row["mean_ho"]}
        per_leg = " | ".join(
            f"{nm}: tr={np.mean([r['scores'][nm]['tr'] for r in cv_rows]):+.3f} ho={np.mean([r['scores'][nm]['ho'] for r in cv_rows]):+.3f}"
            for nm in legs
        )
        print(
            f"  {mask:3s} (L={c_row['n_active']:2d}, {t_mean:5.1f}s/fold, dLOD={dl_mean:.4f}) | "
            f"Canon tr={c_row['mean_tr']:+.3f} ho={c_row['mean_ho']:+.3f} | "
            f"5-Fold tr={tr_mean:+.3f} ho={ho_mean:+.3f} || {per_leg}",
            flush=True
        )

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
x = np.arange(len(MASKS))
bw = 0.36
mlbls = ["ALL\n(52L, 104p)", "F32 (+9.61d)\n(16L, 32p)", "Z20 (+9.61d)\n(10L, 20p)", "F30 (no 9.61d)\n(15L, 30p)", "Z18 (no 9.61d)\n(9L, 18p)"]

for ax, fam, title in [(ax1, "ATL3", "Atlantic Family: NAO + Baltic_3pt + AMO"),
                        (ax2, "IP2", "Indo-Pacific Family: NINO4 + PDO")]:
    trs = [summary[fam][m]["tr"] for m in MASKS]
    hos = [summary[fam][m]["ho"] for m in MASKS]
    ts  = [summary[fam][m]["t"]  for m in MASKS]
    ax.bar(x - bw/2, trs, bw, color="#1f77b4", alpha=0.88, label="5-Fold Mean Train CC")
    ax.bar(x + bw/2, hos, bw, color="#d95f02", alpha=0.90, label="5-Fold Mean OOS Held CC")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x); ax.set_xticklabels(mlbls, fontsize=8.3)
    ax.set_title(title, fontsize=10.5)
    ax.grid(axis="y", alpha=0.3)
    for i in range(len(MASKS)):
        ax.text(i - bw/2, trs[i] + 0.015, f"{trs[i]:+.2f}\n({ts[i]:.0f}s)", ha="center", va="bottom", fontsize=7.3)
        yt = hos[i] + 0.015 if hos[i] >= 0 else hos[i] - 0.065
        ax.text(i + bw/2, yt, f"{hos[i]:+.2f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold")
ax1.set_ylabel("5-Fold Mean Correlation CC (1970–2020)")
ax1.set_ylim(-0.20, 0.98)
ax2.legend(loc="upper right", fontsize=8.5)
fig.suptitle("Phase 1 Compact Gravitational Tidal Basis Evaluation (Cold dLOD + Homotopy, 5-Fold OOS CV)", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out_png = HERE / "abc_compact_tidal_comparison.png"
fig.savefig(out_png, dpi=120)
plt.close(fig)
print("\nSaved", out_png, flush=True)
