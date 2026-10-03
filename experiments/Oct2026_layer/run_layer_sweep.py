#!/usr/bin/env python3
"""LAYER=1 sweep of all 89 SST quads (2026-10-03), seeded from the user's Feb2026 LAYER=1
fits of amo, pdo and nino4 (each seed = that index's lt.exe.p AND lt.exe.resp).

Per quad:
  stage 1  for each seed: CC search (TIMEOUT s), then a TEST_ONLY re-score. A result is
           usable if dLOD >= 0.99. Choose the usable one with the best mean(train, validate);
           the 2000-05 test gap stays held out of the choice.
  stage 2  "higher harmonics as needed": from the stage-1 winner, add ONE more integer
           harmonic of the last regressed winding (NH gets another "1", harm gets the
           smallest multiplier not already used, starting at 2) and search again. Kept only
           if train rises >= 0.02, validate drops <= 0.03 and dLOD >= 0.99.
Settings: the seed resp (LAYER 1, NM, NH, IDATE/INIT_DATE 1880) plus EXCLUDE TRUE with the
2000-05 test gap, VALIDATE TRUE, ACCEL FALSE; env as in the seeds' lte_run.sh.
Data and dlod_ref from Feb2026_sweep/<quad>. Frozen binary ./lt.exe (copy of Feb2026/enso_opt).
Resumable (results.jsonl). Usage: nohup setsid python3 run_layer_sweep.py > sweep.log 2>&1 &
Env: LS_TIMEOUT (200), LS_WORKERS (4), LS_THREADS (4), LS_ONLY (comma list of quads)."""
import json, os, re, shutil, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; EXP = HERE.parent; FEB = EXP / "Feb2026"; DATA = EXP / "Feb2026_sweep"
sys.path.insert(0, str(EXP / "Sep2026")); import sweep, sweep_cc_refine as r
TIMEOUT = int(os.environ.get("LS_TIMEOUT", 200)); WORKERS = int(os.environ.get("LS_WORKERS", 4)); THREADS = os.environ.get("LS_THREADS", "4")
SEEDS = ("amo", "pdo", "nino4"); DLOD_MIN = 0.99
STD_RESP = {"EXCLUDE": "TRUE", "TRAIN_START": "2000.0", "TRAIN_END": "2005.0", "ACCEL": "FALSE", "VALIDATE": "TRUE"}
KEEP = ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv", "lt.exe.windings.json", "dlod_compare.csv", "dlod_ref.dat", "run.out")
LOCK = threading.Lock(); RESULTS = HERE / "results.jsonl"; RUNS = HERE / "_runs"
ENV = dict(re.findall(r'^export (\w+)=(\S+)', (FEB / "amo" / "lte_run.sh").read_text(), re.M))
def log(m):
    with LOCK: print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)
def make(q, tag, p_src, resp_src):
    d = RUNS / f"{q}__{tag}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    shutil.copy(DATA / q / f"{q}.dat", d / f"{d.name}.dat"); shutil.copy(DATA / q / "dlod_ref.dat", d / "dlod_ref.dat")
    shutil.copy(p_src, d / "lt.exe.p"); shutil.copy(resp_src, d / "lt.exe.resp")
    for k, v in STD_RESP.items(): sweep.set_resp_key(d / "lt.exe.resp", k, v)
    return d
def ov(d, test_only):
    return {**ENV, "CLIMATE_INDEX": f"{d.name}.dat", "ACCEL": "FALSE", "TEST_ONLY": "true" if test_only else "false",
            "TIMEOUT": "30" if test_only else str(TIMEOUT), "NUMBER_OF_PROCESSORS": "2" if test_only else THREADS}
def search_and_score(d):
    r.run_lt(d, ov(d, False), False, TIMEOUT)
    out = r.parse(r.run_lt(d, ov(d, True), True, 30)); t = out["triplet"]
    j = json.loads((d / "lt.exe.windings.json").read_text()); p = json.loads((d / "lt.exe.p").read_text())
    return dict(dir=d.name, train=t["train"], validate=t["validate"], test=t["test"], dlod=out["dlod"],
                usable=bool(out["dlod"] is not None and out["dlod"] >= DLOD_MIN), ltep=p["ltep"], harm=p.get("harm"),
                offs=p["offs"], bg=p["bg"], impA=p["impA"], impB=p["impB"], kap=j["k_amp_phase"])
def nh_of(resp):
    m = [l for l in resp.read_text().splitlines() if re.match(r'^NH\s', l)]
    return re.search(r'"([^"]*)"', m[-1]).group(1).split() if m else []
def add_harmonic(d):
    """one more integer harmonic of the last regressed winding"""
    rp = d / "lt.exe.resp"; nh = nh_of(rp) + ["1"]
    txt = re.sub(r'^NH\s.*\n', '', rp.read_text(), flags=re.M); rp.write_text(txt + f'NH                             "{" ".join(nh)}"\n')
    p = json.loads((d / "lt.exe.p").read_text()); h = [int(x) for x in (p.get("harm") or [])]
    h = h[:len(nh) - 1]; h.append(next(m for m in range(2, 50) if m not in h)); p["harm"] = [float(x) for x in h]
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return h
def install(q, best):
    src = RUNS / best["dir"]; dst = HERE / q; dst.mkdir(exist_ok=True)
    shutil.copy(DATA / q / f"{q}.dat", dst / f"{q}.dat")
    for f in KEEP:
        if (src / f).exists(): shutil.copy2(src / f, dst / f)
    sh = (dst / "lte_run.sh").read_text().replace(f"{src.name}.dat", f"{q}.dat").replace("TEST_ONLY=true", "TEST_ONLY=false")
    (dst / "lte_run.sh").write_text(sh)
def solve(q):
    t0 = time.time(); rs = []
    for s in SEEDS:
        try: rs.append(dict(seed=s, **search_and_score(make(q, s, FEB / s / "lt.exe.p", FEB / s / "lt.exe.resp"))))
        except Exception as exc: rs.append(dict(seed=s, error=repr(exc)[:200], usable=False))
    pool = [x for x in rs if x.get("usable")] or [x for x in rs if "train" in x]
    rec = dict(quad=q, seeds=rs)
    if not pool:
        rec.update(chosen=None, flag="no usable result")
    else:
        best = max(pool, key=lambda x: (x["train"] + x["validate"]) / 2)
        rec.update(chosen=best["seed"], flag=None if best.get("usable") else f"dLOD {best['dlod']} < {DLOD_MIN} (best kept)")
        try:   # stage 2: one more harmonic
            src = RUNS / best["dir"]; d = make(q, best["seed"] + "_h", src / "lt.exe.p", src / "lt.exe.resp")
            h = add_harmonic(d); h2 = search_and_score(d); rec["harmonic_try"] = dict(added=h, **h2)
            if h2["usable"] and h2["train"] >= best["train"] + 0.02 and h2["validate"] >= best["validate"] - 0.03:
                best = dict(seed=best["seed"] + "+h", **h2); rec["chosen"] = best["seed"]
        except Exception as exc: rec["harmonic_try"] = dict(error=repr(exc)[:200])
        install(q, best); rec["final"] = best
    rec["elapsed_s"] = round(time.time() - t0)
    with LOCK, RESULTS.open("a") as fh: fh.write(json.dumps(rec) + "\n")
    f = rec.get("final") or {}
    log(f"{q}: {rec.get('chosen')} train/val/test {f.get('train', 0):.3f}/{f.get('validate', 0):.3f}/{f.get('test', 0):.3f} dLOD {f.get('dlod')} "
        f"ltep {[round(x, 4) for x in f.get('ltep', [])]} {rec.get('flag') or ''} ({rec['elapsed_s']}s)")
    return rec
def main():
    RUNS.mkdir(exist_ok=True)
    if not (HERE / "lt.exe").exists(): shutil.copy2(FEB / "enso_opt", HERE / "lt.exe")
    for link, tgt in (("lt.exe", HERE / "lt.exe"), ("dlod3.dat", FEB / "dlod3.dat")):
        if not (RUNS / link).exists(): (RUNS / link).symlink_to(tgt)
    quads = sorted(p.name for p in DATA.iterdir() if p.is_dir() and sweep.parse_grid_name(p.name))
    if os.environ.get("LS_ONLY"): quads = [q for q in quads if q in os.environ["LS_ONLY"].split(",")]
    done = {json.loads(l)["quad"] for l in RESULTS.read_text().splitlines() if l} if RESULTS.exists() else set()
    todo = [q for q in quads if q not in done]
    log(f"{len(quads)} quads, {len(done)} done, {len(todo)} to run; seeds {SEEDS}; TIMEOUT {TIMEOUT}s, {WORKERS} workers x {THREADS} threads; dLOD >= {DLOD_MIN}")
    with ThreadPoolExecutor(WORKERS) as ex: list(ex.map(solve, todo))
    (HERE / "SWEEP_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"); log("DONE")
if __name__ == "__main__": main()
