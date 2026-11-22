#!/usr/bin/env python3
"""Gradient refinement of every quad's own Oct2026_quad5 fit (2026-11-20): the "near-optimal seed" speed-up test.
For each of the 89 centres:
  1. <centre>/ is a copy of Oct2026_quad5/<centre> (settings, data, neighbours), with lt.exe.p = that run's final
     result.p; TEST_ONLY scores it as it stands (seed_scores.out: train / held-out 2000-2005 per leg);
  2. lpap_refine.py --apply --jac <centre>/quad_jac.npz --lam LAM refines the nonzero lpap lines (trust region on the
     solver's own training score, 8 rounds); the derivative is the
     quad's own (its non-lpap parameters differ from every other quad's); lt.exe.p is replaced only when the
     solver's training composite improves (previous kept as lt.exe.p.pre_refine);
  3. TEST_ONLY scores the final lt.exe.p (run.out, legs.png); result.p = the final lt.exe.p; DONE marks the centre.
JOBS centres at a time; a centre with DONE is skipped, so the sweep resumes.
    setsid nohup python3 quad5_refine_own.py > driver.out 2>&1 &      (log: sweep.log; SWEEP_DONE when finished)
"""
import os, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE = Path(__file__).resolve().parent; E = HERE.parent; OLD = E / "Oct2026_quad5"
REFINE = E / "Oct2026_generator/legs/lpap_refine.py"
JOBS, LAM, ITERS, WORKERS = 2, 1.0, 48, 7
def log(msg):
    with open(HERE / "sweep.log", "a") as f: f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
def score(d, out):
    with open(d / out, "w") as f:
        subprocess.run(["bash", "lte_run.sh", "60", "5"], cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                       env=dict(os.environ, TEST_ONLY="TRUE"))
def run(q):
    d = HERE / q
    if (d / "DONE").exists(): log(f"{q}: done already, skipped"); return
    src = OLD / q
    if not (src / "result.p").exists(): log(f"{q}: no result.p in {src}, skipped"); return
    if d.exists(): d.rename(HERE / f"{q}.partial.{int(time.time())}")
    d.mkdir()
    for f in src.iterdir():
        if f.is_file() and (f.name.startswith(("lt.exe.resp", "lte_run.sh")) or f.suffix == ".dat") and f.name != "dlod3.dat":
            shutil.copy(f, d)
    shutil.copy(src / "result.p", d / "lt.exe.p"); shutil.copy(src / "result.p", d / "seed.p")
    t0 = time.time(); score(d, "seed_scores.out")
    r = subprocess.run([sys.executable, str(REFINE), str(d), "--apply", "--jac", str(d / "quad_jac.npz"), "--lam", str(LAM),
                        "--iters", str(ITERS), "--workers", str(WORKERS)], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="2"))
    (d / "refine.out").write_text(r.stdout + r.stderr)
    log(f"{q}: refined ({time.time() - t0:.0f} s)\n" + "\n".join(l for l in r.stdout.splitlines() if l.startswith(("baseline", "solver, training mean", "applied", "not applied")))
        + (("\n" + r.stderr.strip()[-400:]) if r.returncode else ""))
    score(d, "run.out"); shutil.copy(d / "lt.exe.p", d / "result.p")
    (d / "DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n")); log(f"{q}: final scores written ({time.time() - t0:.0f} s)")
if __name__ == "__main__":
    for f, src in (("lt.exe", E / "Feb2026/enso_opt"), ("leg_scores2.py", E / "Oct2026_generator/legs/leg_scores2.py"), ("dlod3.dat", E / "Feb2026/dlod3.dat")):
        if not (HERE / f).exists(): (HERE / f).symlink_to(src)
    qs = sys.argv[1:] or sorted(p.name for p in OLD.iterdir() if p.is_dir() and p.name.startswith("k") and (p / "result.p").exists())
    log(f"sweep start: {len(qs)} centres from {OLD}/<centre>/result.p, own derivative per centre, lam {LAM}, {ITERS} evaluations, {JOBS} at a time")
    with ThreadPoolExecutor(JOBS) as ex: list(ex.map(run, qs))
    if len(sys.argv) == 1: (HERE / "SWEEP_DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log("sweep done")
