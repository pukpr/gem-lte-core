#!/usr/bin/env python3
"""Stage 1 of a progressive re-fit of the Kaplan quads with basin seeds from the 6-leg climate-index fit (2026-11-21).
30 of the 89 quads (sample.json: farthest-point spread over the globe, each with its basin). Each quad is the centre of
a five-leg run with its Oct2026_quad5 neighbours, seeded from the 6-leg fit
(legs/...leg6.tna.extrapolate) merged with one leg's own parameters: Atlantic <- AMO (besk 0.2007), Pacific <- NINO4
(besk 0.2306), Indian <- IODE (besk 0.2306). Settings: the 6-leg lt.exe.resp (template.resp, LEG_WINDINGS/LEG_P off,
so the five legs share one parameter set, besk included), 2000-2005 held out, ACCEL and MANIFOLD_PRE on.
Search TIMEOUT s on THREADS threads, JOBS at a time; a quad whose legs.png exists is skipped (resumes).
    setsid nohup python3 quad5_besk_sweep.py > driver.out 2>&1 &      (log: sweep.log; SWEEP_DONE)
"""
import json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE = Path(__file__).resolve().parent; E = HERE.parent
DLOD = E / "Oct2026_dlod"; MAKE = E / "Oct2026_generator/legs/make_leg_run.py"
TIMEOUT, THREADS, JOBS = 600, 5, 3
sample = json.loads((HERE / "sample.json").read_text()); nb = json.loads((E / "Oct2026_quad5/neighbours.json").read_text())
def log(msg):
    with open(HERE / "sweep.log", "a") as f: f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
def run(q):
    d = HERE / q
    if (d / "legs.png").exists(): log(f"{q}: done already, skipped"); return
    args = [sys.executable, str(MAKE), q, "--root", str(HERE), "--seed", str(HERE / f"seed_{sample[q]}.p"), "--template", str(HERE / "template.resp"),
            "--seed-test", "--extra-env", "ACCEL=TRUE MANIFOLD_PRE=TRUE", "--leg", f"{DLOD / q / (q + '.dat')}:C_{q[1:]}"]
    for n in nb[q]: args += ["--leg", f"{DLOD / n['quad'] / (n['quad'] + '.dat')}:{n['side']}_{n['quad'][1:]}"]
    if d.exists(): args.append("--force")
    out = subprocess.run(args, capture_output=True, text=True)
    log(f"{q} ({sample[q]}): set up\n" + out.stdout.strip() + (("\n" + out.stderr.strip()) if out.returncode else ""))
    if out.returncode: return
    t0 = time.time()
    with open(d / "run.out", "w") as f:
        subprocess.run(["bash", "lte_run.sh", str(TIMEOUT), str(THREADS)], cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    if (d / "lt.exe.p").exists(): (d / "result.p").write_bytes((d / "lt.exe.p").read_bytes())
    log(f"{q}: search finished in {time.time() - t0:.0f} s; legs.png {'written' if (d / 'legs.png').exists() else 'MISSING'}")
if __name__ == "__main__":
    for f, src in (("lt.exe", E / "Feb2026/enso_opt"), ("leg_scores2.py", E / "Oct2026_generator/legs/leg_scores2.py"), ("dlod3.dat", E / "Feb2026/dlod3.dat")):
        if not (HERE / f).exists(): (HERE / f).symlink_to(src)
    qs = sys.argv[1:] or sorted(sample)
    log(f"stage 1 start: {len(qs)} quads, {TIMEOUT} s x {THREADS} threads, {JOBS} at a time")
    with ThreadPoolExecutor(JOBS) as ex: list(ex.map(run, qs))
    if len(qs) == len(sample): (HERE / "SWEEP_DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log("stage 1 done")
