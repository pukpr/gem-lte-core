#!/usr/bin/env python3
"""Paired pass over the 89 LAYER=1 quads (2026-10-03): for each quad, from its Oct2026_layer fit,
  ctrl  300 s CC search, IR = 0 (symmetric, i.e. off)
  warm  300 s CC search, IR_SIDE = 1, IR seeded at -0.2 (warm-side-only delay differential)
then a TEST_ONLY re-score of each. Results installed separately in ctrl/<quad> and warm/<quad>
so both sets can be compared globally (held-out 2000-05, 1880 back-cast) at equal effort.
Pacific quads first. Frozen binary ./lt.exe (Feb2026/enso_opt 11:47 build). Resumable.
Usage: nohup setsid python3 run_ir_pass.py > pass.log 2>&1 &"""
import json, os, re, shutil, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; EXP = HERE.parent; SRC = EXP / "Oct2026_layer"; FEB = EXP / "Feb2026"
sys.path.insert(0, str(EXP / "Sep2026")); import sweep, sweep_cc_refine as r
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026"); import global_sst_from_quads as g
TIMEOUT = int(os.environ.get("IRP_TIMEOUT", 300)); WORKERS = 4; THREADS = "4"
KEEP = ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv", "lt.exe.windings.json", "dlod_compare.csv", "dlod_ref.dat")
LOCK = threading.Lock(); RESULTS = HERE / "results.jsonl"; RUNS = HERE / "_runs"
def log(m):
    with LOCK: print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)
def is_pacific(q):
    lat, lon = g.quad_centre(q); return abs(lat) < 60 and (lon >= 130 or lon <= -75)
def run_variant(q, tag, ir, side):
    d = RUNS / f"{q}__{tag}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(SRC / q / f, d / f)
    shutil.copy(SRC / q / f"{q}.dat", d / f"{d.name}.dat")
    p = json.loads((SRC / q / "lt.exe.p").read_text()); p["IR"] = ir; (d / "lt.exe.p").write_text(json.dumps(p, indent=2))
    sweep.set_resp_key(d / "lt.exe.resp", "IR_SIDE", str(side))
    env = {**dict(re.findall(r'^export (\w+)=(\S+)', (SRC / q / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{d.name}.dat"}
    r.run_lt(d, {**env, "TEST_ONLY": "false", "TIMEOUT": str(TIMEOUT), "NUMBER_OF_PROCESSORS": THREADS}, False, TIMEOUT)
    out = r.parse(r.run_lt(d, {**env, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x = a[:, 0], a[:, 1], a[:, 2]
    k = (x != 0) & (t >= 1950); tt = t[k] - t[k].mean()
    D = x[k] - np.polyval(np.polyfit(tt, x[k], 1), tt); M = m[k] - np.polyval(np.polyfit(tt, m[k], 1), tt)
    hi, lo = D >= np.quantile(D, .9), D <= np.quantile(D, .1)
    dst = HERE / tag / q; dst.mkdir(parents=True, exist_ok=True); shutil.copy(SRC / q / f"{q}.dat", dst / f"{q}.dat")
    for f in KEEP:
        if (d / f).exists(): shutil.copy2(d / f, dst / f)
    (dst / "lte_run.sh").write_text((dst / "lte_run.sh").read_text().replace(f"{d.name}.dat", f"{q}.dat").replace("TEST_ONLY=true", "TEST_ONLY=false"))
    t3 = out["triplet"]
    return dict(train=t3["train"], validate=t3["validate"], test=t3["test"], dlod=out["dlod"], IR=json.loads((d / "lt.exe.p").read_text())["IR"],
                warm=float(M[hi].mean() / D[hi].mean()), cold=float(M[lo].mean() / D[lo].mean()))
def solve(q):
    t0 = time.time(); rec = dict(quad=q, pacific=is_pacific(q))
    for tag, ir, side in (("ctrl", 0.0, 0), ("warm", -0.2, 1)):
        try: rec[tag] = run_variant(q, tag, ir, side)
        except Exception as exc: rec[tag] = dict(error=repr(exc)[:200])
    rec["elapsed_s"] = round(time.time() - t0)
    with LOCK, RESULTS.open("a") as fh: fh.write(json.dumps(rec) + "\n")
    c, w = rec["ctrl"], rec["warm"]
    if "train" in c and "train" in w:
        log(f"{q}{' (Pac)' if rec['pacific'] else ''}: ctrl {c['train']:.3f}/{c['validate']:.3f}/{c['test']:+.3f} | warm {w['train']:.3f}/{w['validate']:.3f}/{w['test']:+.3f} IR {w['IR']:+.3f} dLOD {w['dlod']:.4f} warmcap {c['warm']:.2f}->{w['warm']:.2f} ({rec['elapsed_s']}s)")
    else: log(f"{q}: error {c.get('error') or w.get('error')}")
def main():
    RUNS.mkdir(exist_ok=True)
    if not (HERE / "lt.exe").exists(): shutil.copy2(FEB / "enso_opt", HERE / "lt.exe")
    for link, tgt in (("lt.exe", HERE / "lt.exe"), ("dlod3.dat", FEB / "dlod3.dat")):
        if not (RUNS / link).exists(): (RUNS / link).symlink_to(tgt)
    quads = sorted(p.name for p in SRC.iterdir() if p.is_dir() and sweep.parse_grid_name(p.name))
    quads.sort(key=lambda q: (not is_pacific(q), q))
    done = {json.loads(l)["quad"] for l in RESULTS.read_text().splitlines() if l} if RESULTS.exists() else set()
    todo = [q for q in quads if q not in done]
    log(f"{len(quads)} quads ({sum(map(is_pacific, quads))} Pacific first), {len(todo)} to run; {TIMEOUT}s per variant, {WORKERS} workers x {THREADS} threads")
    with ThreadPoolExecutor(WORKERS) as ex: list(ex.map(solve, todo))
    (HERE / "PASS_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"); log("DONE")
if __name__ == "__main__": main()
