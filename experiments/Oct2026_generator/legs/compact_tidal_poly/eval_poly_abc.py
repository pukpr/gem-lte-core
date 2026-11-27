#!/usr/bin/env python3
"""Phase 1 Systematic A/B/C Evaluation of Compact Gravitational Tidal Representations:
  - ALL (52 lines, 104 params)
  - F30 (15 lines: 9 zonal m=0 + 6 diurnal/semidiurnal m=1,2, 30 params)
  - Z18 ( 9 lines: zonal m=0 only, 18 params)

Families (strictly Cold dLOD + Homotopy, REFINE_LAM=0.30, zero seed leakage):
  1. Atlantic Family: NAO + baltic_3pt_na + AMO (3 legs, LEG_P=PARAMS, NM=3, NH=10, MAXH=45)
     Canonical 1960-1970 + 5-Fold Rolling OOS CV (1970-2020, five 10-yr folds)
  2. Indo-Pacific Family: NINO4 + PDO (2 legs, LEG_P=PARAMS, NM=3, NH=7, MAXH=30)
     Canonical 2000-2005 + 5-Fold Rolling OOS CV (1970-2020, five 10-yr folds)
"""
import concurrent.futures
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
import numpy as np

LEGS_ROOT = Path("/home/paul/gemini/gem-lte-core/experiments/Oct2026_generator/legs")
ATL = LEGS_ROOT / "atlantic_family"
IND7 = LEGS_ROOT / "7indices"
WORK = LEGS_ROOT / "compact_tidal_poly" / "abc_eval"
WORK.mkdir(parents=True, exist_ok=True)
if not (WORK / "dlod3.dat").exists():
    os.symlink(LEGS_ROOT / "dlod3.dat", WORK / "dlod3.dat")

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
        f"-FSTART                        0.1\n"
        f"FLIP                           0.0015\n"
        f"SPREAD_MIN                     0.000001\n"
        f"SPREAD_MAX                     0.45\n"
        f"MAXH                           {maxh}\n"
        f'NH                             "{nh_str}"\n'
        f"NM                             3\n"
        f"YEAR                           0.00405305\n"
        f"DECAY                          48.0\n"
        f"DTW                            3\n"
        f"IDATE                          1880.0\n"
        f"INIT_DATE                      1880.0\n"
        f"STRICT_IDATE                   TRUE\n"
        f"THRESHOLD_ACTION               CATCHUP\n"
        f"PART                           FALSE\n"
        f"FILTER                         0.0\n"
        f"JERK                           0.0\n"
        f"LEGS                           {len(legs)}\n"
        f"LEG_WINDINGS                   TRUE\n"
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
        "DLOD_REF": "FALSE",
        "EXCLUDE": "true",
        "TRAIN_START": str(h0),
        "TRAIN_END": str(h1),
        "F9": "1",
        "IDATE": "1880.0",
        "LOCKA": "FALSE",
        "LOCKT": "FALSE",
        "METRIC": "CC",
        "TREND": "true",
        "TIDE_LEGACY": "FALSE",
        "NUMBER_OF_PROCESSORS": "4",
        "TIMEOUT": "0",
        "TEST_ONLY": "FALSE",
        "REFINE": "TRUE",
        "REFINE_ANCHOR_SEED": "FALSE",
        "REFINE_HOMOTOPY": "TRUE",
        "REFINE_LAM": "0.30",
        "REFINE_ROUNDS": "4",
        "REFINE_ITERS": "40",
        "REFINE_STAGE1": "0",
        "REFINE_TIDAL_MASK": mask,
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

    # Verify active non-zero lines in saved lt.exe.p
    p_fin = json.loads((d / "lt.exe.p").read_text())
    n_active = sum(1 for row in p_fin["lpap"] if abs(row[1]) > 1e-12)
    m_dl = re.search(r"([+-]?\d+\.\d+):dLOD:", r2.stdout)
    dlod_cc = float(m_dl.group(1)) if m_dl else float("nan")

    scores = {}
    trs, hos = [], []
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

FOLDS_5 = [("1970_80", 1970, 1980), ("1980_90", 1980, 1990), ("1990_00", 1990, 2000),
           ("2000_10", 2000, 2010), ("2010_20", 2010, 2020)]

jobs = []
for mask in ["ALL", "F30", "Z18"]:
    # Atlantic Family: NAO + baltic_3pt + AMO
    for wname, h0, h1 in [("canon_60_70", 1960, 1970)] + FOLDS_5:
        jobs.append(("ATL3", ["nao", "baltic_3pt", "amo"], H10, 45, mask, wname, h0, h1))
    # Indo-Pacific Family: NINO4 + PDO
    for wname, h0, h1 in [("canon_00_05", 2000, 2005)] + FOLDS_5:
        jobs.append(("IP2", ["nino4", "pdo"], H7, 30, mask, wname, h0, h1))

print("=" * 140, flush=True)
print(f"PHASE 1 A/B/C EVALUATION: ALL (104p) vs F30 (30p) vs Z18 (18p) on ATL3 (NAO+Baltic+AMO) & IP2 (NINO4+PDO) — {len(jobs)} runs (4x parallel)", flush=True)
print("=" * 140, flush=True)

results = {}
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

print("\n" + "=" * 140, flush=True)
print("SUMMARY: 5-FOLD ROLLING OOS CV (1970-2020, five 10-yr folds) + CANONICAL HOLD-OUT", flush=True)
print("=" * 140, flush=True)
for fam, canon_w, legs in [("ATL3", "canon_60_70", ["nao","baltic_3pt","amo"]),
                            ("IP2",  "canon_00_05", ["nino4","pdo"])]:
    print(f"\n--- {fam} ({' + '.join(legs)}) ---", flush=True)
    for mask in ["ALL", "F30", "Z18"]:
        c_row = results[f"{fam}_{mask}_{canon_w}"]
        cv_rows = [results[f"{fam}_{mask}_{wn}"] for wn, _, _ in FOLDS_5]
        t_mean = np.mean([r["time_s"] for r in cv_rows])
        dl_mean = np.mean([r["dlod_cc"] for r in cv_rows])
        tr_mean = np.mean([r["mean_tr"] for r in cv_rows])
        ho_mean = np.mean([r["mean_ho"] for r in cv_rows])
        per_leg = " | ".join(
            f"{nm}: tr={np.mean([r['scores'][nm]['tr'] for r in cv_rows]):+.3f} ho={np.mean([r['scores'][nm]['ho'] for r in cv_rows]):+.3f}"
            for nm in legs
        )
        print(
            f"  {mask:3s} (L={c_row['n_active']:2d}, {t_mean:4.1f}s/fold, dLOD={dl_mean:.4f}) | "
            f"Canon tr={c_row['mean_tr']:+.3f} ho={c_row['mean_ho']:+.3f} | "
            f"5-Fold tr={tr_mean:+.3f} ho={ho_mean:+.3f} || {per_leg}",
            flush=True
        )

(WORK / "abc_results.json").write_text(json.dumps(results, indent=2))
print("\nSaved abc_results.json", flush=True)
