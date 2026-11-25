#!/usr/bin/env python3
"""Hybrid basin-hopping optimizer for single-leg and multi-leg LTE runs (suitable for long / overnight runs).

Alternates:
  Phase A (non-LPAP exploration): optional small perturbation of non-LPAP forcing/Bessel parameters,
           followed by a short stochastic lt.exe search with LOCKT=TRUE (lpap frozen, dLOD preserved).
  Phase B (deterministic polish): lpap_refine.py --full (analytic pre-Bessel Jacobian + Block A lpap +
           Block B per-leg continuous GN + greedy discrete harmonic screening, verified by lt.exe).

Each cycle is accepted strictly on training composite (and optional --min-dlod floor); held-out 1960-1970
CC is logged for monitoring only, never used to select. Best parameters are saved after every accepted cycle
so the run can be stopped at any time or left running overnight.

usage:
  python3 ../basin_hop.py . --cycles 20 --nonlpap-sec 120 --threads 14 --lam 1.0 --rounds 3
"""
import argparse, csv, glob, json, math, os, random, re, shutil, subprocess, sys, time
import numpy as np

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dir", help="run directory (e.g. . or experiments/Oct2026_generator/legs/7indices)")
ap.add_argument("--cycles", type=int, default=10, help="number of outer non-LPAP + deterministic cycles (default 10)")
ap.add_argument("--nonlpap-sec", type=int, default=90, help="seconds of LOCKT=TRUE non-LPAP search per cycle (default 90)")
ap.add_argument("--threads", type=int, default=0, help="threads for LOCKT=TRUE search (default: largest multiple of legs <= cores)")
ap.add_argument("--lam", type=float, default=1.0, help="dLOD penalty passed to lpap_refine.py (default 1.0)")
ap.add_argument("--rounds", type=int, default=3, help="trust-region rounds per cycle in lpap_refine.py (default 3)")
ap.add_argument("--iters", type=int, default=24, help="surrogate evals per cycle in lpap_refine.py (default 24)")
ap.add_argument("--perturb", type=float, default=0.0, help="optional relative log-normal kick to non-LPAP params on stalled cycles (e.g. 0.02)")
ap.add_argument("--min-dlod", type=float, default=0.97, help="reject cycle if dLOD CC falls below this floor (default 0.97)")
ap.add_argument("--seed", type=int, default=42, help="RNG seed for optional non-LPAP perturbations")
a = ap.parse_args()

D = os.path.abspath(a.dir)
HERE = os.path.dirname(os.path.abspath(__file__))
REFINE = os.path.join(HERE, "lpap_refine.py")
W = os.path.join(D, "basin_hop_work")
BEST = os.path.join(W, "best")
os.makedirs(BEST, exist_ok=True)
rng = random.Random(a.seed)

# Parse exported env from lte_run.sh
env = {}
for line in open(os.path.join(D, "lte_run.sh")):
    m = re.match(r"^\s*export\s+([^#\n]*)", line)
    if m:
        for kv in m.group(1).split():
            if "=" in kv and "${" not in kv:
                k, v = kv.split("=", 1); env[k] = v

_nlegs = re.findall(r"^\s*LEGS\s+(\d+)", open(os.path.join(D, "lt.exe.resp")).read(), re.M)
nlegs = int(_nlegs[-1]) if _nlegs and int(_nlegs[-1]) > 1 else 1
SUF = [f".{L}" for L in range(1, nlegs + 1)] if nlegs > 1 else [""]
ncpu = os.cpu_count() or 8
threads = a.threads or max(nlegs, (ncpu // nlegs) * nlegs)
EXE = os.path.abspath(os.path.join(D, "..", "lt.exe"))
if not (os.path.isfile(EXE) and os.path.getsize(EXE) > 1024):
    EXE = os.path.abspath(os.path.join(HERE, "..", "..", "Feb2026", "enso_opt"))

TS = float(env.get("TRAIN_START", "1960"))
TE = float(env.get("TRAIN_END", "1970"))
EXC = env.get("EXCLUDE", "TRUE").upper() in ("TRUE", "T", "1")

def p_files():
    return [os.path.join(D, "lt.exe.p")] + [os.path.join(D, f"lt.exe.p{s}") for s in SUF if s and os.path.exists(os.path.join(D, f"lt.exe.p{s}"))]

def save_best():
    for f in p_files():
        shutil.copy2(f, os.path.join(BEST, os.path.basename(f)))

def restore_best():
    for f in p_files():
        src = os.path.join(BEST, os.path.basename(f))
        if os.path.exists(src):
            shutil.copy2(src, f)

def run_lt(timeout, nproc, lockt=False, test_only=False, log_path=None):
    e = dict(os.environ); e.update(env)
    e.update({
        "DLOD_REF": "TRUE", "MANIFOLD_PRE": "TRUE",
        "LOCKT": "TRUE" if lockt else "FALSE",
        "TEST_ONLY": "TRUE" if test_only else "FALSE",
        "TIMEOUT": str(timeout), "NUMBER_OF_PROCESSORS": str(nproc),
    })
    cmd = f"ulimit -s 65536 2>/dev/null || true; exec {EXE}"
    if log_path:
        with open(log_path, "w") as fh:
            subprocess.run(["bash", "-c", cmd], cwd=D, env=e, stdout=fh, stderr=subprocess.STDOUT)
        return open(log_path).read()
    r = subprocess.run(["bash", "-c", cmd], cwd=D, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return r.stdout

def score_dir(log_text=None):
    trs, hos = [], []
    for s in SUF:
        arr = np.loadtxt(os.path.join(D, f"lte_results{s}.csv"), delimiter=",")
        t, m, o = arr[:, 0], arr[:, 1], arr[:, 2]
        inside = (t >= TS) & (t < TE)
        tr = (~inside if EXC else inside) & (o != 0)
        ho = (inside if EXC else ~inside) & (o != 0)
        trs.append(float(np.corrcoef(m[tr], o[tr])[0, 1]))
        hos.append(float(np.corrcoef(m[ho], o[ho])[0, 1]) if ho.sum() > 2 else float("nan"))
    dlod = float("nan")
    if log_text:
        ms = re.findall(r"([\d.]+):dLOD:", log_text)
        if ms: dlod = float(ms[-1])
    return trs, hos, dlod

NON_LPAP_KEYS = ["ma", "mp", "mb", "shfT", "init", "asym", "impA", "impB", "besk", "offs", "bg", "IR"]

def perturb_non_lpap(scale):
    """Apply small log-normal perturbations to a random subset of non-LPAP scalars per leg."""
    if scale <= 0.0: return
    for idx_l, s in enumerate(SUF):
        fp = os.path.join(D, f"lt.exe.p{s}" if s else "lt.exe.p")
        if not os.path.exists(fp): continue
        p = json.load(open(fp))
        chosen = rng.sample(NON_LPAP_KEYS, k=min(3, len(NON_LPAP_KEYS)))
        for k in chosen:
            if k in p and isinstance(p[k], (int, float)) and p[k] != 0.0:
                p[k] = float(p[k]) * math.exp(rng.gauss(0.0, scale))
        json.dump(p, open(fp, "w"), indent=2)
        if idx_l == 0 and s:
            pm = json.load(open(os.path.join(D, "lt.exe.p")))
            for k in chosen:
                if k in p: pm[k] = p[k]
            json.dump(pm, open(os.path.join(D, "lt.exe.p"), "w"), indent=2)

# ---- Baseline ----
t0 = time.time()
out0 = run_lt(10, nlegs, lockt=True, test_only=True, log_path=os.path.join(W, "baseline.out"))
trs0, hos0, dlod0 = score_dir(out0)
best_tr, best_ho, best_dlod = float(np.mean(trs0)), float(np.mean(hos0)), dlod0
save_best()
print(f"[cycle 0] baseline | train {' '.join(f'{x:.4f}' for x in trs0)} (mean {best_tr:.4f}) | ho mean {best_ho:+.4f} | dLOD {best_dlod:.4f}", flush=True)

hist_path = os.path.join(W, "history.csv")
with open(hist_path, "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["cycle", "elapsed_s", "status", "mean_train", "mean_ho", "dlod"] + [f"tr_{i+1}" for i in range(nlegs)] + [f"ho_{i+1}" for i in range(nlegs)])
    cw.writerow([0, 0.0, "baseline", f"{best_tr:.6f}", f"{best_ho:.6f}", f"{best_dlod:.6f}"] + [f"{x:.5f}" for x in trs0] + [f"{x:.5f}" for x in hos0])

stalls = 0
for cyc in range(1, a.cycles + 1):
    tc = time.time()
    restore_best()
    # Optional kick if previous cycle stalled and --perturb > 0
    if stalls > 0 and a.perturb > 0.0:
        kick = a.perturb * min(2.0, 1.0 + 0.5 * (stalls - 1))
        perturb_non_lpap(kick)
        print(f"[cycle {cyc}] applied non-LPAP perturbation (scale {kick:.3g} after {stalls} stall(s))", flush=True)
    # Phase A: LOCKT=TRUE stochastic search over non-LPAP params
    out_a = run_lt(a.nonlpap_sec, threads, lockt=True, test_only=False, log_path=os.path.join(W, f"cyc{cyc}_nonlpap.out"))
    trs_a, hos_a, dlod_a = score_dir(out_a)
    print(f"[cycle {cyc}] Phase A (LOCKT=TRUE {a.nonlpap_sec}s): train mean {np.mean(trs_a):.4f}, ho mean {np.mean(hos_a):+.4f}, dLOD {dlod_a:.4f}", flush=True)
    # Phase B: deterministic lpap_refine.py --full
    ref_log = os.path.join(W, f"cyc{cyc}_refine.out")
    cmd_b = [sys.executable, "-u", REFINE, D, "--full", "--lam", str(a.lam), "--rounds", str(a.rounds), "--iters", str(a.iters), "--apply"]
    with open(ref_log, "w") as fh:
        subprocess.run(cmd_b, stdout=fh, stderr=subprocess.STDOUT)
    out_b = run_lt(10, nlegs, lockt=True, test_only=True, log_path=os.path.join(W, f"cyc{cyc}_final.out"))
    trs_b, hos_b, dlod_b = score_dir(out_b)
    m_tr, m_ho = float(np.mean(trs_b)), float(np.mean(hos_b))
    ok = (m_tr > best_tr + 1e-6) and (math.isnan(dlod_b) or dlod_b >= a.min_dlod)
    status = "accepted" if ok else "rejected"
    if ok:
        best_tr, best_ho, best_dlod = m_tr, m_ho, dlod_b
        save_best(); stalls = 0
    else:
        restore_best(); stalls += 1
    dt = time.time() - tc
    print(f"[cycle {cyc}] Phase B (--full): train {' '.join(f'{x:.4f}' for x in trs_b)} (mean {m_tr:.4f}) | ho mean {m_ho:+.4f} | dLOD {dlod_b:.4f} -> {status} [{dt:.0f}s]", flush=True)
    with open(hist_path, "a", newline="") as fh:
        csv.writer(fh).writerow([cyc, f"{time.time()-t0:.1f}", status, f"{m_tr:.6f}", f"{m_ho:.6f}", f"{dlod_b:.6f}"] + [f"{x:.5f}" for x in trs_b] + [f"{x:.5f}" for x in hos_b])

# Ensure final files on disk match best and refresh legs.png
restore_best()
out_f = run_lt(10, nlegs, lockt=True, test_only=True, log_path=os.path.join(W, "final.out"))
trs_f, hos_f, dlod_f = score_dir(out_f)
print(f"\n[done] best mean train {np.mean(trs_f):.4f} (start {np.mean(trs0):.4f}, {np.mean(trs_f)-np.mean(trs0):+.4f}) | mean ho {np.mean(hos_f):+.4f} (start {np.mean(hos0):+.4f}) | dLOD {dlod_f:.4f} | total {time.time()-t0:.0f}s", flush=True)
