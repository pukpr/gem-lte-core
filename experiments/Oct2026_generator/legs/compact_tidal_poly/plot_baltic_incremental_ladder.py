#!/usr/bin/env python3
"""Direct visual incremental-ladder comparison on single-leg baltic_3pt_na.dat (1880-2025, 1960-1970 held out):
  Z18 (9 zonal lines) -> Z20 (+9.6137d, 10 lines) -> Z30 (+draconic/T/Mtm, 15 lines) -> ALL (52 lines)
Plots 6 stacked panels matching baltic_0710_atl4_cold_H10.png:
  1-4: baltic_3pt_na standardized fit for Z18, Z20, Z30, ALL
  5:   Post-Bessel manifolds overlaid
  6:   Pre-Bessel manifolds overlaid (MANIFOLD_PRE=TRUE) with dLOD ref CC
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
HERE = LEGS_ROOT / "compact_tidal_poly"
WORK = HERE / "baltic_ladder"
WORK.mkdir(parents=True, exist_ok=True)
if not (WORK / "dlod3.dat").exists():
    os.symlink(LEGS_ROOT / "dlod3.dat", WORK / "dlod3.dat")

ATL_P = json.loads((ATL / "lt.exe.p").read_text())
H10 = [2, 3, 4, 6, 9, 13, 18, 24, 31, 41]

def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])

def run_tier(mask):
    d = WORK / f"b3_{mask}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    os.symlink(ATL / "baltic_3pt_na.dat", d / "baltic_3pt.dat")
    nh_str = " ".join(["1"] * len(H10))
    (d / "lt.exe.resp").write_text(
        f"-FSTART                        0.1\nFLIP                           0.0015\n"
        f"SPREAD_MIN                     0.000001\nSPREAD_MAX                     0.45\n"
        f"MAXH                           45\n"
        f'NH                             "{nh_str}"\nNM                             3\n'
        f"YEAR                           0.00405305\nDECAY                          48.0\n"
        f"DTW                            3\nIDATE                          1880.0\n"
        f"INIT_DATE                      1880.0\nSTRICT_IDATE                   TRUE\n"
        f"THRESHOLD_ACTION               CATCHUP\nPART                           FALSE\n"
        f"FILTER                         0.0\nJERK                           0.0\n"
        f"CLIMATE_INDEX                  baltic_3pt.dat\n"
    )
    p = json.loads(json.dumps(ATL_P))
    p["harm"] = [float(h) for h in H10]
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2))

    env = os.environ.copy()
    env.update({
        "DLOD_REF": "FALSE", "EXCLUDE": "true", "TRAIN_START": "1960", "TRAIN_END": "1970",
        "F9": "1", "IDATE": "1880.0", "LOCKA": "FALSE", "LOCKT": "FALSE", "METRIC": "CC",
        "TREND": "true", "TIDE_LEGACY": "FALSE", "NUMBER_OF_PROCESSORS": "4",
        "TIMEOUT": "0", "TEST_ONLY": "FALSE", "REFINE": "TRUE", "REFINE_ANCHOR_SEED": "FALSE",
        "REFINE_HOMOTOPY": "TRUE", "REFINE_LAM": "0.30", "REFINE_ROUNDS": "4",
        "REFINE_ITERS": "60", "REFINE_STAGE1": "0", "REFINE_TIDAL_MASK": mask,
    })
    t0 = time.time()
    r1 = subprocess.run("ulimit -s 65536 && ../../../lt.exe", shell=True, cwd=d, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=240)
    dt = time.time() - t0
    (d / "refine.log").write_text(r1.stdout)

    # Post-Bessel score
    env["DLOD_REF"] = "TRUE"; env["REFINE"] = "FALSE"; env["TEST_ONLY"] = "TRUE"; env["MANIFOLD_PRE"] = "FALSE"
    r2 = subprocess.run("ulimit -s 65536 && ../../../lt.exe", shell=True, cwd=d, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
    (d / "score_post.log").write_text(r2.stdout)
    a_post = np.loadtxt(d / "lte_results.csv", delimiter=",")

    # Pre-Bessel score
    env["MANIFOLD_PRE"] = "TRUE"
    r3 = subprocess.run("ulimit -s 65536 && ../../../lt.exe", shell=True, cwd=d, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
    a_pre = np.loadtxt(d / "lte_results.csv", delimiter=",")

    m_dl = re.search(r"([+-]?\d+\.\d+):dLOD:", r2.stdout)
    dlod_cc = float(m_dl.group(1)) if m_dl else float("nan")
    p_fin = json.loads((d / "lt.exe.p").read_text())
    n_act = sum(1 for row in p_fin["lpap"] if abs(row[1]) > 1e-12)
    t, m, y, F_post = a_post[:, 0], a_post[:, 1], a_post[:, 2], a_post[:, 3]
    F_pre = a_pre[:, 3]
    ok = y != 0; ho = ok & (t >= 1960) & (t < 1970); tr = ok & ~ho
    return {
        "mask": mask, "n_act": n_act, "time_s": dt, "dlod_cc": dlod_cc,
        "tr": corr(m[tr], y[tr]), "ho": corr(m[ho], y[ho]),
        "t": t, "m": m, "y": y, "F_post": F_post, "F_pre": F_pre,
    }

TIERS = [
    ("Z18", "Z18 (9 zonal lines, 18p: Mf, Mm, Mt, MSm, MSf, 2D, Ssa)", "#1baf7a"),
    ("Z20", "Z20 (10 zonal lines, 20p: Z18 + 9.6137d secular arch)",     "#2a78d6"),
    ("Z30", "Z30 (15 zonal lines, 30p: Z20 + T, D, 27.44d, 26.99d, Mtm)", "#9467bd"),
    ("ALL", "ALL (52 empirical lines, 104p: full table baseline)",       "#eb6834"),
]

res_map = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    for r in ex.map(run_tier, [x[0] for x in TIERS]):
        res_map[r["mask"]] = r
        print(f"{r['mask']:4s} (L={r['n_act']:2d}, {r['time_s']:4.1f}s) | Baltic tr={r['tr']:+.4f} ho(60-70)={r['ho']:+.4f} | dLOD={r['dlod_cc']:.4f}", flush=True)

fig, axes = plt.subplots(6, 1, figsize=(14.5, 18.2), sharex=True)
F_post_all = res_map["ALL"]["F_post"]
F_pre_all  = res_map["ALL"]["F_pre"]

for i, (mask, desc, col) in enumerate(TIERS):
    r = res_map[mask]
    t, m, y = r["t"], r["m"], r["y"]
    ok = y != 0
    y_std = (y[ok] - np.mean(y[ok])) / np.std(y[ok])
    m_std = (m[ok] - np.mean(m[ok])) / np.std(m[ok])
    ax = axes[i]
    ax.plot(t[ok], y_std, color="#2a78d6", lw=0.8, alpha=0.85, label="baltic_3pt_na (standardized)")
    ax.plot(t[ok], m_std, color="#eb6834", lw=0.9, alpha=0.92, label=f"{mask} model (standardized)")
    ax.axvspan(1960, 1970, color="#e5e5e5", alpha=0.8)
    ax.set_title(
        f"{desc}  |  train {r['tr']:+.3f}, held 1960–1970 {r['ho']:+.3f}  |  dLOD {r['dlod_cc']:.4f}  ({r['time_s']:.0f}s)",
        loc="left", fontsize=9.8,
    )
    ax.legend(loc="upper left", fontsize=8.2, ncol=2)
    ax.set_xlim(1880, 2026)

    cc_post = corr(r["F_post"], F_post_all)
    cc_pre  = corr(r["F_pre"],  F_pre_all)
    axes[4].plot(t, r["F_post"], color=col, lw=1.0, alpha=0.88,
                 label=f"{mask} (L={r['n_act']}, CC vs ALL {cc_post:.4f})")
    axes[5].plot(t, r["F_pre"], color=col, lw=1.0, alpha=0.88,
                 label=f"{mask} (L={r['n_act']}, CC vs ALL {cc_pre:.4f}, dLOD {r['dlod_cc']:.3f})")

axes[4].axvspan(1960, 1970, color="#e5e5e5", alpha=0.8)
axes[4].set_title("Post-Bessel manifolds across incremental tidal tiers (Z18 → Z20 → Z30 → ALL)", loc="left", fontsize=10.2)
axes[4].legend(loc="upper left", fontsize=8.2, ncol=2)

axes[5].axvspan(1960, 1970, color="#e5e5e5", alpha=0.8)
axes[5].set_title("Pre-Bessel manifolds across incremental tidal tiers (MANIFOLD_PRE=TRUE)", loc="left", fontsize=10.2)
axes[5].legend(loc="upper left", fontsize=8.2, ncol=2)
axes[5].set_xlabel("Year (1880–2025; shaded: held-out 1960–1970)")

fig.tight_layout()
out_png = HERE / "baltic_3pt_incremental_ladder.png"
fig.savefig(out_png, dpi=115)
plt.close(fig)
print("Saved", out_png)
