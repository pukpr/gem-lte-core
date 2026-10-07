#!/usr/bin/env python3
"""Arm B (DECLARATION.md): full 300 s lt.exe CC search for the real series (draw 0) and 19 IAAFT
surrogates of each index, all starting from the same cross-seeded parameters (AMO <- PDO fit,
PDO <- AMO fit, NINO4 <- AMO fit; the index's own resp with NM/NH from the seed), then TEST_ONLY
re-score. Frozen binary ../Oct2026_layer_ir/lt.exe. Resumable; results in arm_b.jsonl."""
import json, re, shutil, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent; E = H.parent; S = H / "surr"; RUNS = H / "_runs"; OUT = H / "arm_b.jsonl"
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
SEEDOF = {"amo": "pdo", "pdo": "amo", "nino4": "amo"}; LOCK = threading.Lock()
def nh_line(resp): return [l for l in resp.read_text().splitlines() if re.match(r'^NH\s', l)][-1]
def run(job):
    i, d = job; name = f"{i}_s{d:02d}"; w = RUNS / name
    if w.exists(): shutil.rmtree(w)
    w.mkdir(parents=True); src = H / "seeds" / i; seed = H / "seeds" / SEEDOF[i]
    shutil.copy(seed / "lt.exe.p", w / "lt.exe.p"); shutil.copy(src / "lt.exe.resp", w / "lt.exe.resp"); shutil.copy(src / "dlod_ref.dat", w / "dlod_ref.dat")
    shutil.copy(S / f"{name}.dat", w / f"{name}.dat")
    rp = w / "lt.exe.resp"; nm = re.search(r'^NM\s+(\S+)', (seed / "lt.exe.resp").read_text(), re.M).group(1)
    txt = re.sub(r'^NH\s.*\n', '', rp.read_text(), flags=re.M); rp.write_text(txt + nh_line(seed / "lt.exe.resp") + "\n"); sweep.set_resp_key(rp, "NM", nm)
    env = {**dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{name}.dat"}
    t0 = time.time()
    try:
        r.run_lt(w, {**env, "TEST_ONLY": "false", "TIMEOUT": "300", "NUMBER_OF_PROCESSORS": "4"}, False, 300)
        out = r.parse(r.run_lt(w, {**env, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30)); t3 = out["triplet"]
        rec = dict(index=i, draw=d, train=t3["train"], validate=t3["validate"], test=t3["test"], dlod=out["dlod"], ltep=json.loads((w / "lt.exe.p").read_text())["ltep"])
    except Exception as exc: rec = dict(index=i, draw=d, error=repr(exc)[:200])
    rec["elapsed_s"] = round(time.time() - t0)
    with LOCK, OUT.open("a") as fh: fh.write(json.dumps(rec) + "\n")
    print(time.strftime("%H:%M:%S"), json.dumps(rec), flush=True)
if __name__ == "__main__":
    RUNS.mkdir(exist_ok=True)
    for f, tgt in (("lt.exe", E / "Oct2026_layer_ir" / "lt.exe"), ("dlod3.dat", E / "Feb2026" / "dlod3.dat")):
        if not (RUNS / f).exists(): (RUNS / f).symlink_to(tgt)
    done = {(json.loads(l)["index"], json.loads(l)["draw"]) for l in OUT.read_text().splitlines() if l} if OUT.exists() else set()
    jobs = [(i, d) for d in range(20) for i in ("amo", "pdo", "nino4") if (i, d) not in done]
    with ThreadPoolExecutor(4) as ex: list(ex.map(run, jobs))
    (H / "ARM_B_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"); print("DONE", flush=True)
