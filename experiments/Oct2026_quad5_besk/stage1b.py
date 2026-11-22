#!/usr/bin/env python3
"""Stage 1b (2026-11-21): stage 1 redone with the manifold integrated from 1880 (template.resp now has STRICT_IDATE TRUE
and INIT_DATE 1880, as Oct2026_quad5; stage 1's template lacked them, so its quads started the integrator in 1950).
  --score   every sample quad scored (TEST_ONLY, no search) with both the Atlantic seed (AMO leg, besk 0.2007) and the
            Indo-Pacific seed (NINO4 leg, besk 0.2306): s1b/<quad>_<atl|pac>/
  --search  each quad searched TIMEOUT s from the seed with the better TRAINING composite (s1b/<quad>_best/)
Runs go to s1b/; logs to s1b/stage1b.log.   python3 stage1b.py --score ;  python3 stage1b.py --search"""
import json, os, re, statistics as st, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE = Path(__file__).resolve().parent; S1B = HERE / "s1b"; E = HERE.parent
DLOD = E / "Oct2026_dlod"; MAKE = E / "Oct2026_generator/legs/make_leg_run.py"
SEEDS = {"atl": HERE / "seed_atlantic.p", "pac": HERE / "seed_pacific.p"}
TIMEOUT, THREADS, JOBS = 600, 5, 3
sample = json.loads((HERE / "sample.json").read_text()); nb = json.loads((E / "Oct2026_quad5/neighbours.json").read_text())
pat = re.compile(r"^(\S+): \S+ \| train ([+-][\d.]+).*?held out \S+ ([+-][\d.]+)", re.M)
def log(m):
    with open(S1B / "stage1b.log", "a") as f: f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + m + "\n")
def setup(q, tag, seed):
    args = [sys.executable, str(MAKE), f"{q}_{tag}", "--root", str(S1B), "--seed", str(seed), "--template", str(HERE / "template.resp"),
            "--extra-env", "ACCEL=TRUE MANIFOLD_PRE=TRUE", "--force", "--leg", f"{DLOD / q / (q + '.dat')}:C_{q[1:]}"]
    for n in nb[q]: args += ["--leg", f"{DLOD / n['quad'] / (n['quad'] + '.dat')}:{n['side']}_{n['quad'][1:]}"]
    r = subprocess.run(args, capture_output=True, text=True); return r.returncode == 0
def score(d, out="scores.out", threads=2):
    with open(d / out, "w") as f:
        subprocess.run(["bash", "lte_run.sh", "60", str(threads)], cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=dict(os.environ, TEST_ONLY="TRUE"))
    s = pat.findall((d / out).read_text()); return [(float(a), float(b)) for _, a, b in s]
def do_score(qt):
    q, tag = qt; d = S1B / f"{q}_{tag}"
    if not setup(q, tag, SEEDS[tag]): log(f"{q}_{tag}: setup failed"); return
    s = score(d); log(f"{q}_{tag}: train {st.mean(a for a, _ in s):.4f} held-out {st.mean(b for _, b in s):+.4f}")
def do_search(q):
    sc = {}
    for tag in SEEDS:
        s = pat.findall((S1B / f"{q}_{tag}" / "scores.out").read_text()); sc[tag] = st.mean(float(a) for _, a, _ in s)
    tag = max(sc, key=sc.get); d = S1B / f"{q}_best"
    if (d / "legs.png").exists(): return
    setup(q, "best", SEEDS[tag]); (d / "seed_choice.txt").write_text(f"{tag} (training {sc['atl']:.4f} atl, {sc['pac']:.4f} pac)\n")
    with open(d / "run.out", "w") as f:
        subprocess.run(["bash", "lte_run.sh", str(TIMEOUT), str(THREADS)], cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    (d / "result.p").write_bytes((d / "lt.exe.p").read_bytes()); log(f"{q}_best ({tag}): search finished")
if __name__ == "__main__":
    qs = sorted(sample)
    if "--score" in sys.argv:
        log("score start"); 
        with ThreadPoolExecutor(6) as ex: list(ex.map(do_score, [(q, t) for q in qs for t in SEEDS]))
        log("score done")
    if "--search" in sys.argv:
        log("search start")
        with ThreadPoolExecutor(JOBS) as ex: list(ex.map(do_search, qs))
        (S1B / "SEARCH_DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n")); log("search done")
