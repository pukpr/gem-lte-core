#!/usr/bin/env python3
"""Deterministic gradient refinement of a run's lpap (the tidal amplitudes and phases), as a speed-up for seeds that
are already close to optimal. It replaces the long random-walk tail of a search, not the search for a good basin:
from dLOD or the crowd set it finds only the nearest, much weaker valley (see LEGS.md).

usage: lpap_refine.py DIR [--lam 1] [--iters 40] [--rounds 4] [--workers N] [--jac FILE] [--apply]

DIR is a run folder (single leg, or legs with lte_results.<leg>.csv) whose lte_run.sh exports the run's settings.
Nothing in DIR changes unless --apply is given; work files go to DIR/refine_work/.

  1. Baseline: the run is scored as it stands (TEST_ONLY, MANIFOLD_PRE TRUE) in a copy.
  2. Surrogate: each leg's model is rebuilt from the saved manifold exactly as the solver builds it (each leg's own
     Bessel step and windings; least squares of the IR-matched data on 1, F, the windings, trend, acceleration,
     annual, semi-annual and CW terms; the annual impulse; IR on the model; the 9-point filter). The run stops if
     any leg's rebuilt model differs from its saved column 2 by more than 1e-6, i.e. for settings it does not cover.
  3. Derivative: the pre-Bessel manifold is exactly linear in each lpap line's (A cos, A sin). Each leg's derivative
     for all lines comes from perturbed TEST_ONLY runs of the whole run (one run per coefficient gives every leg),
     cached in refine_work/jac.npz; it depends on the non-lpap parameters, settings and dates only.
  4. Refinement: trust-region least squares on all nonzero lines (zero lines stay zero), maximizing the mean training CC of the legs, with a
     penalty LAM * |change of the daily dLOD tidal series|^2 / |series|^2 (the dLOD-calibrated lines as reference),
     so lines move freely where dLOD cannot see them. Steps are scaled by each line's aliased leverage on the manifold.
  5. Trust region on the solver's own score: each round's bounded step is scored by the real solver; a better
     training composite is accepted (every leg re-anchored on its actual manifold, step x1.5), a worse one rejected
     (step x0.5). Needed because the forcing is linear except the integrator's sign-dependent friction term, which
     offsets the manifold when a zero crossing moves. With --apply and a better composite than the start,
     DIR/lt.exe.p is backed up to lt.exe.p.pre_refine and its lpap replaced.
Deterministic: the same inputs give the same result.
"""
import argparse, glob, hashlib, json, os, re, shutil, subprocess, sys, time
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from scipy.optimize import least_squares

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dir"); ap.add_argument("--lam", type=float, default=1.0); ap.add_argument("--iters", type=int, default=48, help="surrogate evaluations in all, split over the rounds")
ap.add_argument("--rounds", type=int, default=8, help="trust-region rounds, each checked by the real solver")
ap.add_argument("--accept", choices=["train", "validate"], default="train",
                help="validate: a round is accepted only if the training composite improves AND the mean validation CC "
                     "(legs' CC over VALIDATE_START on, from the round's lte_results) does not fall by more than --vtol")
ap.add_argument("--vtol", type=float, default=0.0, help="allowed fall of the validation score per accepted round (default 0)")
ap.add_argument("--step", type=float, default=0.3, help="initial trust radius per coefficient, in leverage-scaled units (~0.05 of manifold)")
ap.add_argument("--workers", type=int, default=0, help="parallel solver runs for the derivative (default: cores / legs)"); ap.add_argument("--apply", action="store_true")
ap.add_argument("--jac", default=None, help="shared derivative cache (.npz) for runs whose legs all share the forcing parameters "
                "(no lt.exe.p.<leg> files): one single-leg derivative serves every leg, and every run with the same "
                "parameters, settings and dates (e.g. a sweep from one seed)")
a = ap.parse_args()
D = os.path.abspath(a.dir); W = os.path.join(D, "refine_work"); os.makedirs(W, exist_ok=True)

# ---- run settings ---------------------------------------------------------------------------------------------
env = {}
for line in open(os.path.join(D, "lte_run.sh")):
    m = re.match(r"^\s*export\s+([^#\n]*)", line)
    if m:
        for kv in m.group(1).split():
            if "=" in kv and "${" not in kv: k, v = kv.split("=", 1); env[k] = v
EXE = os.path.abspath(os.path.join(D, "..", "lt.exe"))   # keep the lt.exe name: the binary reads <its name>.p
dl = os.path.realpath(os.path.join(D, "dlod3.dat")) if os.path.exists(os.path.join(D, "dlod3.dat")) else os.path.realpath(os.path.join(D, "..", "dlod3.dat"))
shutil.copy(dl, os.path.join(W, "dlod3.dat"))
_nlegs = re.findall(r"^\s*LEGS\s+(\d+)", open(os.path.join(D, "lt.exe.resp")).read(), re.M) if os.path.exists(os.path.join(D, "lt.exe.resp")) else []
legs = list(range(1, int(_nlegs[-1]) + 1)) if _nlegs and int(_nlegs[-1]) > 1 else []   # LEGS in lt.exe.resp; per-leg files .1 .. .N
SUF = [f".{L}" for L in legs] if legs else [""]
def resp_index(suf):
    for f in ([os.path.join(D, f"lt.exe.resp{suf}")] if suf else []) + [os.path.join(D, "lt.exe.resp")]:
        if os.path.exists(f):
            m = re.findall(r"^\s*CLIMATE_INDEX\s+(\S+)", open(f).read(), re.M)
            if m: return m[-1]
    return env.get("CLIMATE_INDEX")
def test_only(d, extra=None, threads=2):
    e = dict(os.environ); e.update(env); e.update({"TEST_ONLY": "TRUE", "MANIFOLD_PRE": "TRUE", "NUMBER_OF_PROCESSORS": str(threads)})
    #  every run here scores the lpap it is given: DLOD_REF <name> would re-seed lpap from wisdom_lpap at start-up and
    #  overwrite it (FALSE/CAL likewise), so it is always TRUE here (an explicit DLOD_REF in extra still wins)
    e["DLOD_REF"] = "TRUE"
    e.pop("TIMEOUT", None); e.update(extra or {})
    with open(os.path.join(d, "test.out"), "w") as fh:
        subprocess.run(["bash", "-c", f"ulimit -s 65536; exec {EXE}"], cwd=d, env=e, stdout=fh, stderr=subprocess.STDOUT)
def copy_run(dst, lpap=None):
    shutil.rmtree(dst, ignore_errors=True); os.makedirs(dst)
    for f in glob.glob(os.path.join(D, "lt.exe.p")) + glob.glob(os.path.join(D, "lt.exe.p.[0-9]*")) + glob.glob(os.path.join(D, "lt.exe.resp*")):
        if re.search(r"(lt\.exe\.p(\.\d+)?|lt\.exe\.resp(\.\d+)?)$", f): shutil.copy(f, dst)
    for suf in SUF: shutil.copy(os.path.join(D, resp_index(suf)), dst)
    if lpap is not None:
        p = json.load(open(os.path.join(dst, "lt.exe.p"))); p["lpap"] = [[float(x) for x in r] for r in lpap]; json.dump(p, open(os.path.join(dst, "lt.exe.p"), "w"), indent=1)
def run_scores(d):
    s = open(os.path.join(d, "test.out")).read()
    own = {int(m[0]): float(m[1]) for m in re.findall(r"leg (\d+) pair \d+\s+own:\s+(\S+)", s)}
    cc = re.findall(r"^CC\s+(\S+)\s+(\S+)\s+(\d+)", s, re.M); dlod = re.findall(r"([\d.]+):dLOD:", s); yr = re.findall(r"([\d.]+)=Yr", s)
    if not own and cc: own = {1: float(cc[-1][1])}
    return own, float(dlod[-1]) if dlod else float("nan"), float(yr[-1]) if yr else 365.2422

# ---- 1. baseline ----------------------------------------------------------------------------------------------
t0 = time.time(); B = os.path.join(W, "baseline"); copy_run(B); test_only(B, threads=max(2, len(SUF)))
own0, dlod0, YL = run_scores(B)
P0 = json.load(open(os.path.join(B, "lt.exe.p"))); LP0 = np.array(json.load(open(os.path.join(D, "lt.exe.p")))["lpap"])
print(f"baseline: training {' '.join(f'{own0[k]:.4f}' for k in sorted(own0))} (mean {np.mean(list(own0.values())):.4f}), dLOD {dlod0:.4f}  [{time.time()-t0:.0f} s]")

# ---- 2. surrogate --------------------------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Feb2026"))
from lte_forward import filter9point as filter9
F9 = int(float(env.get("F9", "1"))); EXC = env.get("EXCLUDE", "FALSE").upper() in ("TRUE", "T", "1")
class Leg:
    def __init__(s, suf):
        r = np.loadtxt(os.path.join(B, f"lte_results{suf}.csv"), delimiter=","); s.t, s.model, s.data, s.F = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
        w = json.load(open(os.path.join(B, f"lt.exe.windings{suf}.json")))
        s.m, s.ir, s.ks = w["manifold"], float(w["IR"]), [k for k, _, _ in w["k_amp_phase"]]
        p = dict(P0)
        if suf and os.path.exists(os.path.join(B, f"lt.exe.p{suf}")): p.update(json.load(open(os.path.join(B, f"lt.exe.p{suf}"))))
        s.nm = len(p.get("ltep", [0, 0, 0])); s.impc, s.delb, s.cw = float(p.get("impC", 0)), float(p.get("delB", 0)), float(p.get("cw1", 0)) != 0.0
        ts, te = float(env.get("TRAIN_START", s.t[0])), float(env.get("TRAIN_END", s.t[-1]))
        first = np.where(s.t <= ts + 1e-9)[0][-1] if (s.t <= ts + 1e-9).any() else 0; last = np.where(s.t <= te + 1e-9)[0][-1]
        s.tr = np.zeros(len(s.t), bool)
        if EXC: s.tr[:first + 1] = True; s.tr[last:] = True
        else: s.tr[first:last + 1] = True
        s.ref = float(w.get("secular_context", {}).get("accel_ref", s.t[first]))
        cwf = (YL / (27.212220815 / 2.0)) % 1.0
        s.N = np.column_stack([s.t - s.ref, (s.t - s.ref) ** 2] + [g(2 * np.pi * f * s.t) for f in ((1.0, 2.0, cwf) if s.cw else (1.0, 2.0)) for g in (np.sin, np.cos)])
        s.dr = s.data.copy(); s.dr[12:] = s.data[12:] + s.ir * s.data[:-12]
        trunc = np.rint((s.t - np.floor(s.t)) * 12).astype(int); s.imp = np.where(trunc == int(round(abs(s.delb) * 12)) % 12, s.impc, 0.0)
        s.sc = s.tr & (s.data != 0)
    def post(s, F):
        m = s.m; bg = m["bg"]
        k, e2s, e2c = (m["besk"] * (1 - bg), m["offs"], m.get("bg2", bg)) if s.nm == 1 else (m["besk"], m["offs"], bg)
        ph = 2 * np.pi * k * F; return F + m["impA"] * np.sin(ph) + m["impB"] * np.cos(ph) + e2s * m["impA"] * np.sin(2 * ph) + e2c * m["impB"] * np.cos(2 * ph)
    def model_of(s, F):
        Fp = s.post(F); X = np.column_stack([np.ones_like(Fp), Fp] + [g(2 * np.pi * k * Fp) for k in s.ks for g in (np.sin, np.cos)] + [s.N])
        y = X @ np.linalg.lstsq(X[s.tr], s.dr[s.tr], rcond=None)[0] + s.imp
        if s.ir != 0: y = np.concatenate([y[:12], y[12:] - s.ir * y[:-12]])
        for _ in range(F9): y = filter9(y)
        return y
L = [Leg(suf) for suf in SUF]
for i, l in enumerate(L):
    e = float(np.max(np.abs(l.model_of(l.F) - l.model)))
    if e > 1e-6: sys.exit(f"leg {i + 1}: the surrogate does not reproduce the solver's model (max diff {e:.2e}); settings not covered, nothing done")
print(f"surrogate reproduces every leg's model column (max diff <= 1e-6)")

# ---- 3. derivative of each leg's manifold, from the real solver --------------------------------------------------
#  One TEST_ONLY of the whole run (all legs) per perturbed coefficient gives that column for every leg at once; legs
#  differ in their impulse comb and integrator settings, so each leg has its own derivative.
nl = len(LP0)
def ab(lp): return lp[:, 1] * np.cos(lp[:, 2]), lp[:, 1] * np.sin(lp[:, 2])
A0, B0 = ab(LP0); th0 = np.column_stack([A0, B0]).ravel()
def lp_of(th):
    lp = LP0.copy(); A, Bv = th[0::2], th[1::2]; lp[:, 1] = np.hypot(A, Bv); lp[:, 2] = np.arctan2(Bv, A); return lp
O = os.path.join(W, "oracle"); os.makedirs(O, exist_ok=True)
def oracle(lp, tag, extra=None):
    d = os.path.join(O, tag); copy_run(d, lp); test_only(d, extra, threads=max(2, len(SUF)))
    return [np.loadtxt(os.path.join(d, f"lte_results{suf}.csv"), delimiter=",")[:, 3] for suf in SUF], d
shutil.copy(os.path.join(W, "dlod3.dat"), O)
key = hashlib.md5("".join(open(f).read() for f in sorted(glob.glob(os.path.join(D, "lt.exe.p.[0-9]*")) + glob.glob(os.path.join(D, "lt.exe.resp*"))))
                  .encode() + json.dumps({k: v for k, v in json.load(open(os.path.join(D, "lt.exe.p"))).items() if k != "lpap"}, sort_keys=True).encode()
                  + json.dumps([float(r[0]) for r in LP0]).encode() + json.dumps(env, sort_keys=True).encode()).hexdigest()
JF = os.path.join(W, "jac.npz"); workers = a.workers or max(1, (os.cpu_count() or 4) // max(2, len(SUF)))
def shared_jac(path):
    """One single-leg derivative for every leg: exact when no leg has its own forcing parameters."""
    if glob.glob(os.path.join(D, "lt.exe.p.[0-9]*")):
        sys.exit("--jac needs legs without their own lt.exe.p.<leg> files (shared forcing parameters)")
    resp1 = [l for l in open(os.path.join(D, "lt.exe.resp")).read().splitlines() if not re.match(r"^\s*-?(LEGS|LEG_WINDINGS|LEG_P|CLIMATE_INDEX|SELECT_TOP|SELECT_TIE)\b", l)]
    pk = {k: v for k, v in json.load(open(os.path.join(D, "lt.exe.p"))).items() if k != "lpap"}
    k1 = hashlib.md5((json.dumps(pk, sort_keys=True) + "\n".join(resp1) + json.dumps(env, sort_keys=True)
                      + json.dumps([float(r[0]) for r in LP0])).encode()).hexdigest()
    dates = [np.loadtxt(os.path.join(B, f"lte_results{suf}.csv"), delimiter=",")[:, 0] for suf in SUF]
    if os.path.exists(path):
        z = np.load(path)
        if str(z["key"]) == k1:
            ix = {round(x, 4): i for i, x in enumerate(z["t"])}
            if all(all(round(x, 4) in ix for x in t_) for t_ in dates):
                print(f"derivative: shared cache {path}")
                return [z["J"][[ix[round(x, 4)] for x in t_]] for t_ in dates], z["ref"]
    t1 = time.time(); lsuf = max(range(len(SUF)), key=lambda i: len(dates[i]))
    sd = os.path.join(W, "oracle1"); os.makedirs(sd, exist_ok=True); shutil.copy(os.path.join(W, "dlod3.dat"), sd)
    rsp = resp1 + [f"CLIMATE_INDEX                  {resp_index(SUF[lsuf])}"]
    def one(lp, tag, extra=None):
        d = os.path.join(sd, tag); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
        pp = dict(pk); pp["lpap"] = [[float(x) for x in r] for r in lp]; json.dump(pp, open(os.path.join(d, "lt.exe.p"), "w"), indent=1)
        open(os.path.join(d, "lt.exe.resp"), "w").write("\n".join(rsp) + "\n"); shutil.copy(os.path.join(D, resp_index(SUF[lsuf])), d)
        test_only(d, extra, threads=1); r = np.loadtxt(os.path.join(d, "lte_results.csv"), delimiter=","); return r[:, 0], r[:, 3], d
    h = 0.01; lps = [LP0]; tags = ["base"]
    for j in range(nl):
        for w_ in (0, 1):
            th = th0.copy(); th[2 * j + w_] += h; lps.append(lp_of(th)); tags.append(f"j{j}_{w_}")
    with ThreadPoolExecutor(a.workers or (os.cpu_count() or 4)) as ex: res = list(ex.map(lambda i: one(lps[i], tags[i]), range(len(lps))))
    t_, Fb = res[0][0], res[0][1]; Jm = np.column_stack([(r[1] - Fb) / h for r in res[1:]])
    _, _, dd = one(LP0, "dlod_ref", {"DLOD_REF": "FALSE"}); ref = np.array(json.load(open(os.path.join(dd, "lt.exe.p")))["lpap"])
    for i, suf in enumerate(SUF):   # the shared derivative must reproduce every leg's saved manifold at the seed
        ix = {round(x, 4): k for k, x in enumerate(t_)}; e = float(np.max(np.abs(Fb[[ix[round(x, 4)] for x in dates[i]]] - L[i].F)))
        if e > 1e-6: sys.exit(f"leg {i + 1}: its manifold differs from the shared one by {e:.2e}; --jac is not valid for this run")
    np.savez(path, key=k1, t=t_, J=Jm, ref=ref); print(f"derivative: shared, {Jm.shape[1]} columns from {len(lps)} single-leg runs, saved to {path} [{time.time()-t1:.0f} s]")
    ix = {round(x, 4): i for i, x in enumerate(t_)}
    return [Jm[[ix[round(x, 4)] for x in tt]] for tt in dates], ref
if a.jac:
    JL, LPREF = shared_jac(os.path.abspath(a.jac))
elif os.path.exists(JF) and str(np.load(JF)["key"]) == key:
    jz = np.load(JF); JL = [jz[f"J{i}"] for i in range(len(SUF))]; LPREF = jz["ref"]; print("derivative: cached")
else:
    t1 = time.time(); h = 0.01; lps = [LP0]; tags = ["base"]
    for j in range(nl):
        for w_ in (0, 1):
            th = th0.copy(); th[2 * j + w_] += h; lps.append(lp_of(th)); tags.append(f"j{j}_{w_}")
    with ThreadPoolExecutor(workers) as ex: res = list(ex.map(lambda i: oracle(lps[i], tags[i])[0], range(len(lps))))
    JL = [np.column_stack([(r[i] - res[0][i]) / h for r in res[1:]]) for i in range(len(SUF))]
    _, dd = oracle(LP0, "dlod_ref", {"DLOD_REF": "FALSE"}); LPREF = np.array(json.load(open(os.path.join(dd, "lt.exe.p")))["lpap"])
    np.savez(JF, key=key, ref=LPREF, **{f"J{i}": JL[i] for i in range(len(SUF))})
    print(f"derivative: {2 * nl} columns for {len(SUF)} legs from {len(lps)} solver runs [{time.time()-t1:.0f} s]")
for l, Jl in zip(L, JL): l.J = Jl

# ---- 4. refinement -------------------------------------------------------------------------------------------
Ar, Br = ab(LPREF); thr = np.column_stack([Ar, Br]).ravel()
td = np.loadtxt(os.path.join(B, "dlod_compare.csv"), delimiter=",")[::7, 0]; wv = 2 * np.pi * YL / np.abs(LP0[:, 0])
PHI = np.zeros((len(td), 2 * nl)); PHI[:, 0::2] = np.cos(np.outer(td, wv)); PHI[:, 1::2] = -np.sin(np.outer(td, wv)); NP = np.linalg.norm(PHI @ thr)
J = np.vstack(JL)
XS = 0.05 / np.maximum(np.repeat(np.sqrt(np.mean(J[:, 0::2] ** 2 + J[:, 1::2] ** 2, axis=0) / 2), 2), 1e-9)
FREE = np.repeat(np.abs(LP0[:, 1]) > 0, 2)   # zero-amplitude lines stay at zero (the search never moves them either)
def full(x): th = th_b.copy(); th[FREE] = x; return th
def ccs(th):
    out = []
    for l in L:
        y = l.model_of(l.F + l.J @ (th - th_b)); out.append(float(np.corrcoef(y[l.sc], l.data[l.sc])[0, 1]))
    return out
def resid(x):
    th = full(x); r = []
    for l in L:
        y = l.model_of(l.F + l.J @ (th - th_b))[l.sc]; d = l.data[l.sc]
        r.append(((y - y.mean()) / y.std() - (d - d.mean()) / d.std()) / np.sqrt(2 * len(d) * len(L)))   # sum of squares = mean(1 - CC)
    return np.concatenate(r + [np.sqrt(a.lam) * (PHI @ (th - thr)) / NP])
dlod_cc = lambda th: float(np.corrcoef(PHI @ th, PHI @ thr)[0, 1])

# ---- 4-5. refinement as a trust region on the real solver's score ------------------------------------------------
#  The forcing is linear in lpap except for the integrator's sign-dependent friction term (mp): when the manifold
#  crosses zero (e.g. in the run-up from IDATE), moving lpap moves the crossing and offsets the whole manifold, which
#  the derivative cannot see. So each round takes a bounded step (each coefficient within STEP x its leverage-scaled
#  unit of the current best), the real solver scores it, and: better -> accepted, every leg re-anchored on its actual
#  manifold, step x1.5; worse -> back to the best, step x0.5. Ends after --rounds rounds or when the step is tiny.
VS = float(env.get("VALIDATE_START", "nan"))
if a.accept == "validate" and not np.isfinite(VS): sys.exit("--accept validate needs VALIDATE_START in lte_run.sh")
def vscore(dd):   # mean over legs of CC(model, data) on t >= VALIDATE_START (data rows that are present)
    out = []
    for suf in SUF:
        r = np.loadtxt(os.path.join(dd, f"lte_results{suf}.csv"), delimiter=","); m = (r[:, 0] >= VS) & (r[:, 2] != 0)
        if m.sum() > 24: out.append(float(np.corrcoef(r[m, 1], r[m, 2])[0, 1]))
    return float(np.mean(out)) if out else float("nan")
v0 = vscore(B) if a.accept == "validate" else float("nan")
if a.accept == "validate": print(f"validation (t >= {VS:g}): baseline {v0:.4f}")
m0 = np.mean(list(own0.values())); best = (m0, th0.copy(), dlod0, own0); th_b = th0.copy(); step = a.step; best_v = v0
per = max(4, a.iters // a.rounds); anchor = [l.F.copy() for l in L]
for rnd in range(1, a.rounds + 1):
    for l, F_ in zip(L, anchor): l.F = F_.copy()
    th_b = best[1].copy(); t2 = time.time(); c0 = ccs(th_b); x0 = th_b[FREE]; lo, hi = x0 - step * XS[FREE], x0 + step * XS[FREE]
    sol = least_squares(resid, x0, method="trf", x_scale=XS[FREE], bounds=(lo, hi), max_nfev=per, diff_step=1e-6); th1 = full(sol.x); c1 = ccs(th1)
    R = os.path.join(W, f"round{rnd}"); copy_run(R, lp_of(th1)); test_only(R, threads=max(2, len(SUF))); own1, dlod1, _ = run_scores(R)
    m1 = np.mean(list(own1.values())); v1 = vscore(R) if a.accept == "validate" else float("nan")
    ok = m1 > best[0] and (a.accept == "train" or v1 >= best_v - a.vtol)
    print(f"round {rnd} (step {step:.3g}): surrogate {np.mean(c0):.4f} -> {np.mean(c1):.4f} | solver {m1:.4f} vs best {best[0]:.4f} "
          + (f" | validation {v1:.4f} vs {best_v:.4f}" if a.accept == "validate" else "") +
          f" {'accepted' if ok else 'rejected'}, dLOD {dlod1:.4f} [{sol.nfev} evaluations, {time.time()-t2:.0f} s]")
    if ok:
        best = (m1, th1.copy(), dlod1, own1); best_v = v1; step *= 1.5
        anchor = [np.loadtxt(os.path.join(R, f"lte_results{suf}.csv"), delimiter=",")[:, 3] for suf in SUF]
    else:
        step *= 0.5
    if step < 0.01: break
m1, th1, dlod1, own1 = best; LP1 = lp_of(th1)
print(f"solver, training per leg: " + " ".join(f"{own0[k]:.4f}->{own1[k]:.4f}" for k in sorted(own0)))
print(f"solver, training mean {m0:.4f} -> {m1:.4f}; dLOD {dlod0:.4f} -> {dlod1:.4f}  [total {time.time()-t0:.0f} s]")
json.dump({"lpap": LP1.tolist()}, open(os.path.join(W, "refined_lpap.json"), "w"))
if a.apply:
    if m1 > m0:
        shutil.copy(os.path.join(D, "lt.exe.p"), os.path.join(D, "lt.exe.p.pre_refine"))
        p = json.load(open(os.path.join(D, "lt.exe.p"))); p["lpap"] = LP1.tolist(); json.dump(p, open(os.path.join(D, "lt.exe.p"), "w"), indent=2)
        print(f"applied: {D}/lt.exe.p updated (previous saved as lt.exe.p.pre_refine)")
    else:
        print("not applied: the solver's training score did not improve")
else:
    print(f"refined lpap: {W}/refined_lpap.json (use --apply to write it into lt.exe.p)")
