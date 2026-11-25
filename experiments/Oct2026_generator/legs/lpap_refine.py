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
ap.add_argument("--fd-jac", action="store_true", help="force finite-difference solver probes instead of the analytic pre-Bessel Jacobian")
ap.add_argument("--per-leg", action="store_true", help="alternate Block B per-leg continuous GN (Bessel, ltep, IR) with Block A (lpap)")
ap.add_argument("--leg-step", type=float, default=0.05, help="initial relative trust radius for per-leg continuous parameters (default 0.05)")
ap.add_argument("--leg-lam", type=float, default=1e-3, help="L2 drift penalty on per-leg continuous parameters (default 1e-3)")
ap.add_argument("--harm-screen", action="store_true", help="deterministic greedy 1-swap harmonic screening over 2..MAXH (replaces FLIP)")
ap.add_argument("--harm-expand", type=int, default=0, help="append N new greedy harmonic slots to NH from 2..MAXH (updates NH in lt.exe.resp on --apply)")
ap.add_argument("--harm-max-swaps", type=int, default=2, help="max harmonic swaps per leg per round (0 = until local optimum, default 2)")
ap.add_argument("--harm-min-gain", type=float, default=1e-4, help="minimum per-leg training CC gain to accept a harmonic swap/expansion")
ap.add_argument("--homotopy", action="store_true", help="coarse-to-fine frequency continuation (|k|<=0.35 -> 0.9 -> 2.5 -> all) before full rounds to avoid cycle-skipping from cold seeds")
ap.add_argument("--full", action="store_true", help="shorthand for --per-leg --harm-screen")
ap.add_argument("--jac", default=None, help="shared derivative cache (.npz) for runs whose legs all share the forcing parameters "
                "(no lt.exe.p.<leg> files): one single-leg derivative serves every leg, and every run with the same "
                "parameters, settings and dates (e.g. a sweep from one seed)")
a = ap.parse_args()
if a.full: a.per_leg = True; a.harm_screen = True
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
def _update_resp_nh(path, new_nh):
    txt = open(path).read()
    nh_str = '"' + " ".join(["1"] * new_nh) + '"'
    out, done = [], False
    for line in txt.splitlines():
        if re.match(r"^\s*NH\s+", line):
            out.append(f"NH                             {nh_str}"); done = True
        else:
            out.append(line)
    if not done: out.append(f"NH                             {nh_str}")
    open(path, "w").write("\n".join(out) + "\n")
def copy_run(dst, lpap=None, leg_patches=None, new_nh=None):
    shutil.rmtree(dst, ignore_errors=True); os.makedirs(dst)
    for f in glob.glob(os.path.join(D, "lt.exe.p")) + glob.glob(os.path.join(D, "lt.exe.p.[0-9]*")) + glob.glob(os.path.join(D, "lt.exe.resp*")):
        if re.search(r"(lt\.exe\.p(\.\d+)?|lt\.exe\.resp(\.\d+)?)$", f): shutil.copy(f, dst)
    for suf in SUF: shutil.copy(os.path.join(D, resp_index(suf)), dst)
    p = json.load(open(os.path.join(dst, "lt.exe.p"))); changed_main = False
    if lpap is not None:
        p["lpap"] = [[float(x) for x in r] for r in lpap]; changed_main = True
    if leg_patches:
        for idx_l, (suf, patch) in enumerate(zip(SUF, leg_patches)):
            if not patch: continue
            if suf and os.path.exists(os.path.join(dst, f"lt.exe.p{suf}")):
                pl = json.load(open(os.path.join(dst, f"lt.exe.p{suf}")))
                pl.update(patch); json.dump(pl, open(os.path.join(dst, f"lt.exe.p{suf}"), "w"), indent=2)
            if idx_l == 0:
                p.update(patch); changed_main = True
    if changed_main:
        json.dump(p, open(os.path.join(dst, "lt.exe.p"), "w"), indent=1)
    if new_nh is not None:
        _update_resp_nh(os.path.join(dst, "lt.exe.resp"), new_nh)
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
        s.suf = suf
        r = np.loadtxt(os.path.join(B, f"lte_results{suf}.csv"), delimiter=","); s.t, s.model, s.data, s.F = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
        w = json.load(open(os.path.join(B, f"lt.exe.windings{suf}.json")))
        s.m, s.ir, s.ks = dict(w["manifold"]), float(w["IR"]), [k for k, _, _ in w["k_amp_phase"]]
        p = dict(P0)
        if suf and os.path.exists(os.path.join(B, f"lt.exe.p{suf}")): p.update(json.load(open(os.path.join(B, f"lt.exe.p{suf}"))))
        s.ltep = [float(x) for x in p.get("ltep", [])]
        s.harm = [int(round(float(h))) for h in p.get("harm", [])]
        s.nm = len(s.ltep); s.impc, s.delb, s.cw = float(p.get("impC", 0)), float(p.get("delB", 0)), float(p.get("cw1", 0)) != 0.0
        ts, te = float(env.get("TRAIN_START", s.t[0])), float(env.get("TRAIN_END", s.t[-1]))
        first = np.where(s.t <= ts + 1e-9)[0][-1] if (s.t <= ts + 1e-9).any() else 0; last = np.where(s.t <= te + 1e-9)[0][-1]
        s.tr = np.zeros(len(s.t), bool)
        if EXC: s.tr[:first + 1] = True; s.tr[last:] = True
        else: s.tr[first:last + 1] = True
        s.ref = float(w.get("secular_context", {}).get("accel_ref", s.t[first]))
        cwf = (YL / (27.212220815 / 2.0)) % 1.0
        s.N = np.column_stack([s.t - s.ref, (s.t - s.ref) ** 2] + [g(2 * np.pi * f * s.t) for f in ((1.0, 2.0, cwf) if s.cw else (1.0, 2.0)) for g in (np.sin, np.cos)])
        s.sync_dr()
        trunc = np.rint((s.t - np.floor(s.t)) * 12).astype(int); s.imp = np.where(trunc == int(round(abs(s.delb) * 12)) % 12, s.impc, 0.0)
        s.sc = s.tr & (s.data != 0)
        s.orig_snap = s.snapshot()
    def sync_ks(s):
        if s.nm > 0: s.ks = list(s.ltep) + [s.ltep[s.nm - 1] * h for h in s.harm]
    def sync_dr(s):
        s.dr = s.data.copy(); s.dr[12:] = s.data[12:] + s.ir * s.data[:-12]
    def snapshot(s):
        return (dict(s.m), list(s.ltep), list(s.harm), float(s.ir))
    def restore(s, snap):
        s.m, s.ltep, s.harm, s.ir = dict(snap[0]), list(snap[1]), list(snap[2]), float(snap[3])
        s.sync_ks(); s.sync_dr()
    def to_patch(s, include_bessel=True):
        out = {"ltep": [float(x) for x in s.ltep], "harm": [float(h) for h in s.harm], "IR": float(s.ir)}
        if include_bessel:
            for k in ("impA", "impB", "besk", "offs", "bg"):
                if k in s.m: out[k] = float(s.m[k])
            if s.nm == 1 and "bg2" in s.m: out["bg2"] = float(s.m["bg2"])
        return out
    def post(s, F):
        m = s.m; bg = m["bg"]
        k, e2s, e2c = (m["besk"] * (1 - bg), m["offs"], m.get("bg2", bg)) if s.nm == 1 else (m["besk"], m["offs"], bg)
        ph = 2 * np.pi * k * F; return F + m["impA"] * np.sin(ph) + m["impB"] * np.cos(ph) + e2s * m["impA"] * np.sin(2 * ph) + e2c * m["impB"] * np.cos(2 * ph)
    def model_of(s, F, kmax=None):
        Fp = s.post(F)
        ks = s.ks if kmax is None else [k for k in s.ks if abs(k) <= kmax]
        X = np.column_stack([np.ones_like(Fp), Fp] + [g(2 * np.pi * k * Fp) for k in ks for g in (np.sin, np.cos)] + [s.N])
        y = X @ np.linalg.lstsq(X[s.tr], s.dr[s.tr], rcond=None)[0] + s.imp
        if s.ir != 0: y = np.concatenate([y[:12], y[12:] - s.ir * y[:-12]])
        for _ in range(F9): y = filter9(y)
        return y
    def cc_of(s, F):
        y = s.model_of(F); return float(np.corrcoef(y[s.sc], s.data[s.sc])[0, 1])
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
def leg_opt(suf, key, default):
    """Read a config key with Ada precedence: leg .resp > main .resp > env > default."""
    for f in ([os.path.join(D, f"lt.exe.resp{suf}")] if suf else []) + [os.path.join(D, "lt.exe.resp")]:
        if os.path.exists(f):
            m = re.findall(r"^\s*" + re.escape(key) + r"\s+(\S+)", open(f).read(), re.M | re.I)
            if m: return m[-1]
    return env.get(key, default)
def _as_bool(v): return str(v).strip().upper() in ("TRUE", "T", "1")
def analytic_jac():
    """Exact linear tangent of Tide_Sum_Diff -> Impulse_Amplify -> Shape_Impulses -> IIR wrt (A_j cos P_j, A_j sin P_j).
    Verified by one single-column solver spot-check (run concurrently with dlod_ref); returns None to fall back to FD."""
    t1 = time.time()
    tides_mode = str(leg_opt("", "TIDES", "")).strip().upper()
    if tides_mode in ("ZONAL", "BLEND") or _as_bool(leg_opt("", "TIDE_LEGACY", "FALSE")):
        return None
    per = LP0[:, 0].astype(float)
    JL_a = []
    for suf, l in zip(SUF, L):
        p = dict(P0)
        if suf and os.path.exists(os.path.join(B, f"lt.exe.p{suf}")):
            p.update(json.load(open(os.path.join(B, f"lt.exe.p{suf}"))))
        spy = float(leg_opt(suf, "SAMPLING", "12.0"))
        monthly = abs(spy - 12.0) < 1e-6; step_scale = 12.0 / spy
        idate = float(leg_opt(suf, "IDATE", "0.0"))
        init_raw = float(leg_opt(suf, "INIT_DATE", "-1.0"))
        strict = _as_bool(leg_opt(suf, "STRICT_IDATE", "FALSE")) or init_raw >= 0.0
        jerk = float(leg_opt(suf, "JERK", "0.0"))
        lod = float(leg_opt(suf, "LOD_TERM", "0.0")); lod_floor = float(leg_opt(suf, "LOD_FLOOR", "0.0"))
        imp_only = _as_bool(leg_opt(suf, "IMPULSE", "TRUE"))
        mod12 = _as_bool(leg_opt(suf, "IMPULSE_MOD12", "FALSE"))
        imp_exact = _as_bool(leg_opt(suf, "IMPULSE_EXACT", "FALSE"))
        aliased = _as_bool(leg_opt(suf, "ALIAS", "FALSE"))
        dates = l.t
        if strict and idate < dates[0]:
            n_extra = int(np.ceil((dates[0] - idate) * spy))
            t_ext = np.concatenate([dates[0] - np.arange(n_extra, 0, -1) / spy, dates])
        else:
            n_extra = 0; t_ext = dates.copy()
        delA, delB, asym, delC = float(p["delA"]), float(p["delB"]), float(p["asym"]), float(p.get("delC", 0.0))
        spy_i = int(round(spy))
        def imp_base(tv):
            tr = np.rint((tv - np.floor(tv)) * spy).astype(int)
            dp = int(round(abs(delB) * spy))
            if mod12: dp %= spy_i
            dp2 = (dp + spy_i // 2) % spy_i
            v = np.zeros_like(tv)
            if imp_only:
                v[tr == dp] = delA; v[tr == dp2] = asym
            else:
                if delA < 0.0: v[tr == dp2 - 1] = -delA
                v[tr == dp2] = asym
                if delA > 0.0: v[tr == dp2 + 1] = delA
            return v
        def imp_fn(tv):
            v = imp_base(tv)
            if delC != 0.0: v = v - abs(delC) * imp_base(tv - np.copysign(0.25, delC))
            return v
        imp = imp_fn(t_ext)
        t_tide = t_ext.copy()
        if imp_exact and imp_only:
            tau = abs(delB) - np.floor(abs(delB))
            nz = imp != 0.0
            f = t_ext[nz] - np.floor(t_ext[nz])
            d1 = tau - f; d1 -= np.round(d1)
            d2 = tau + 0.5 - f; d2 -= np.round(d2)
            if delC != 0.0:
                d1 = np.where(np.abs(d1) <= np.abs(d2), d1, d2)
                d2 = tau + 0.25 - f; d2 -= np.round(d2 * 2.0) / 2.0
            t_tide[nz] = t_ext[nz] + np.where(np.abs(d1) <= np.abs(d2), d1, d2)
        freqs = YL / per
        if aliased: freqs = np.remainder(freqs + 0.5, 1.0) - 0.5
        lod_f = np.zeros_like(freqs) if lod == 0.0 else (lod / freqs if lod_floor == 0.0 else lod / (freqs * np.sqrt(1.0 + (lod_floor / freqs) ** 2)))
        k_int = float(p["shfT"]) * freqs - lod_f
        ang = 2.0 * np.pi * np.outer(t_tide, freqs)
        C, S = np.cos(ang), np.sin(ang)
        uA = S + C * k_int; uB = C - S * k_int   # Cos_Phase = False (Tide_Legacy = False)
        if jerk != 0.0:
            nw = 1.0 - jerk; dw = jerk * 13.66083077 / per
            uA, uB = nw * uA + dw * uB, -dw * uA + nw * uB
        U = np.empty((len(t_ext), 2 * nl))
        U[:, 0::2] = uA * imp[:, None]; U[:, 1::2] = uB * imp[:, None]
        mB = float(p.get("mb", 0.0))
        if mB != 0.0:
            B_ = min(abs(mB), 0.9999); R2 = B_ if monthly else B_ ** step_scale
            G = np.empty_like(U); G[0] = (1.0 - R2) * U[0]
            for k in range(1, len(t_ext)): G[k] = (1.0 - R2) * U[k] + R2 * G[k - 1]
            U = G
        der = 1.0 - float(p["ma"])
        lag_a = der if monthly else 1.0 - (1.0 - der) * step_scale
        mem = min(max(lag_a, 0.0), 1.0)
        start = idate - 0.1 / spy if strict else idate
        idx = min(max(int(np.searchsorted(t_ext, start, side="right")), 0), len(t_ext) - 1)
        F = np.zeros_like(U)
        for k in range(idx + 1, len(t_ext)): F[k] = U[k] + mem * F[k - 1]
        for k in range(idx, 0, -1): F[k - 1] = (F[k] - U[k]) / mem if mem > 0.0 else -U[k - 1] + mem * F[k]
        JL_a.append(F[n_extra:])
    # Spot-check 1 column (first nonzero line, Re) in parallel with dlod_ref (small h=1e-4 avoids artificial IIR sign flips)
    nz_lines = np.where(np.abs(LP0[:, 1]) > 0)[0]
    j_chk = int(nz_lines[0]) if len(nz_lines) else 0; h = 1e-4
    th_chk = th0.copy(); th_chk[2 * j_chk] += h
    with ThreadPoolExecutor(2) as ex:
        f_chk = ex.submit(oracle, lp_of(th_chk), "analytic_chk")
        f_ref = ex.submit(oracle, LP0, "dlod_ref", {"DLOD_REF": "FALSE"})
        F_chk, _ = f_chk.result(); _, dd = f_ref.result()
    lpref = np.array(json.load(open(os.path.join(dd, "lt.exe.p")))["lpap"])
    max_rel = 0.0
    for i, (l, Ja) in enumerate(zip(L, JL_a)):
        col_fd = (F_chk[i] - l.F) / h; col_an = Ja[:, 2 * j_chk]
        denom = max(float(np.max(np.abs(col_fd))), 1e-12)
        max_rel = max(max_rel, float(np.max(np.abs(col_an - col_fd))) / denom)
    if max_rel > 1e-4:
        print(f"derivative: analytic spot-check failed (max rel diff {max_rel:.2e} > 1e-4); falling back to finite differences")
        return None
    np.savez(JF, key=key, ref=lpref, **{f"J{i}": JL_a[i] for i in range(len(SUF))})
    print(f"derivative: analytic, {2 * nl} columns for {len(SUF)} legs (spot-check max rel diff {max_rel:.2e}) [{time.time()-t1:.0f} s]")
    return JL_a, lpref
if a.jac:
    JL, LPREF = shared_jac(os.path.abspath(a.jac))
elif os.path.exists(JF) and str(np.load(JF)["key"]) == key:
    jz = np.load(JF); JL = [jz[f"J{i}"] for i in range(len(SUF))]; LPREF = jz["ref"]; print("derivative: cached")
else:
    an = None if a.fd_jac else analytic_jac()
    if an is not None:
        JL, LPREF = an
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

# ---- 4b. Block B (per-leg continuous GN) & discrete harmonic screening/expansion -----------------------------
LEG_P_MODE = str(leg_opt("", "LEG_P", "NONE")).strip().upper()
LEG_WIND = _as_bool(leg_opt("", "LEG_WINDINGS", "FALSE")) or len(SUF) == 1
HAS_PER_LEG_BESSEL = len(SUF) == 1 or LEG_P_MODE in ("BESSEL", "PARAMS") or _as_bool(leg_opt("", "LEG_BESSEL", "FALSE"))
MAXH = int(float(leg_opt("", "MAXH", "30")))
DO_LEG = a.per_leg or a.harm_screen or a.harm_expand > 0
if DO_LEG and len(SUF) > 1 and not LEG_WIND:
    sys.exit("--per-leg / --harm-screen / --harm-expand on multi-leg runs requires LEG_WINDINGS TRUE")

def screen_harmonics(l, leg_idx, F_cur, maxh, max_swaps, min_gain, expand_n=0):
    """Greedy deterministic harmonic expansion + 1-swap screening over 2..maxh (replaces stochastic FLIP)."""
    if l.nm == 0 or not l.harm: return False
    Fp = l.post(F_cur); k_base = l.ltep[l.nm - 1]
    if k_base == 0.0: return False
    X_fixed = np.column_stack([np.ones_like(Fp), Fp] + [g(2 * np.pi * k * Fp) for k in l.ltep for g in (np.sin, np.cos)] + [l.N])
    cand = {}
    for h in range(2, maxh + 1):
        kh = h * k_base
        if any(abs(kh - k0) < 1e-3 for k0 in l.ltep): continue
        ang = 2 * np.pi * kh * Fp
        cand[h] = np.column_stack([np.sin(ang), np.cos(ang)])
    d_sc = l.data[l.sc]
    def eval_h(hlist):
        X = np.hstack([X_fixed] + [cand[h] for h in hlist])
        y = X @ np.linalg.lstsq(X[l.tr], l.dr[l.tr], rcond=None)[0] + l.imp
        if l.ir != 0: y = np.concatenate([y[:12], y[12:] - l.ir * y[:-12]])
        for _ in range(F9): y = filter9(y)
        return float(np.corrcoef(y[l.sc], d_sc)[0, 1])
    h_cur = list(l.harm); h_start = list(h_cur)
    # Ensure existing harmonics are in cand even if outside 2..maxh
    for h in h_cur:
        if h not in cand:
            ang = 2 * np.pi * (h * k_base) * Fp
            cand[h] = np.column_stack([np.sin(ang), np.cos(ang)])
    cc_cur = eval_h(h_cur); cc_start = cc_cur; changed = False
    # Optional HARM expansion: greedily append up to expand_n new harmonic slots
    for _ in range(expand_n):
        used = set(h_cur); best_h, best_cc = None, cc_cur
        for h in range(2, maxh + 1):
            if h in used or h not in cand: continue
            c = eval_h(h_cur + [h])
            if c > best_cc: best_cc, best_h = c, h
        if best_h is not None and (best_cc - cc_cur) >= min_gain:
            h_cur.append(best_h); cc_cur = best_cc; changed = True
        else:
            break
    # Greedy 1-swap screening over existing slots
    swaps = 0
    while max_swaps == 0 or swaps < max_swaps:
        used = set(h_cur); best_move, best_cc = None, cc_cur
        for pos in range(len(h_cur)):
            old_h = h_cur[pos]
            for h in range(2, maxh + 1):
                if h in used or h not in cand: continue
                trial = list(h_cur); trial[pos] = h
                c = eval_h(trial)
                if c > best_cc: best_cc, best_move = c, (pos, h)
        if best_move is not None and (best_cc - cc_cur) >= min_gain:
            pos, h = best_move; h_cur[pos] = h; cc_cur = best_cc; swaps += 1; changed = True
        else:
            break
    if changed:
        l.harm = h_cur; l.sync_ks()
        print(f"  leg {leg_idx + 1} harm: {h_start} -> {h_cur} ({cc_start:.4f} -> {cc_cur:.4f}, {cc_cur - cc_start:+.4f})")
    return changed

def refine_per_leg(l, F_cur, leg_step, leg_lam, max_nfev=12):
    """Block B: bounded trust-region GN on one leg's post-manifold continuous params (Bessel, ltep, IR)."""
    keys = []
    if HAS_PER_LEG_BESSEL:
        for k in ("impA", "impB", "besk", "offs", "bg"):
            if abs(float(l.m.get(k, 0.0))) > 0.0: keys.append(("m", k))
        if l.nm == 1 and abs(float(l.m.get("bg2", 0.0))) > 0.0: keys.append(("m", "bg2"))
    for j in range(l.nm):
        if abs(l.ltep[j]) > 0.0: keys.append(("ltep", j))
    if abs(l.ir) > 0.0: keys.append(("ir", 0))
    if not keys: return
    def get_vec(snap=None):
        m_, lt_, _, ir_ = snap if snap is not None else (l.m, l.ltep, l.harm, l.ir)
        v = []
        for kind, k in keys:
            if kind == "m": v.append(float(m_[k]))
            elif kind == "ltep": v.append(float(lt_[k]))
            else: v.append(float(ir_))
        return np.array(v, dtype=float)
    def set_vec(v):
        for (kind, k), val in zip(keys, v):
            if kind == "m": l.m[k] = float(val)
            elif kind == "ltep": l.ltep[k] = float(val)
            else: l.ir = float(val)
        l.sync_ks(); l.sync_dr()
    z_orig = get_vec(l.orig_snap); z0 = get_vec()
    zs = np.maximum(np.abs(z_orig), 1e-4)
    lo, hi = z0 - leg_step * zs, z0 + leg_step * zs
    min_lt = float(leg_opt(l.suf, "MIN_LT", "0.01")) * 1.0001
    for i_k, (kind, _) in enumerate(keys):
        if kind == "ltep":
            if z0[i_k] >= 0.0: lo[i_k] = max(lo[i_k], min_lt); hi[i_k] = max(hi[i_k], lo[i_k] + 1e-6)
            else: hi[i_k] = min(hi[i_k], -min_lt); lo[i_k] = min(lo[i_k], hi[i_k] - 1e-6)
    d_sc = l.data[l.sc]; d_std = (d_sc - d_sc.mean()) / d_sc.std(); inv_n = 1.0 / np.sqrt(2 * len(d_sc))
    cc0 = l.cc_of(F_cur); snap0 = l.snapshot()
    def lresid(z):
        set_vec(z); y = l.model_of(F_cur)[l.sc]
        ys = (y - y.mean()) / max(float(y.std()), 1e-12)
        r_cc = (ys - d_std) * inv_n
        if leg_lam > 0.0:
            return np.concatenate([r_cc, np.sqrt(leg_lam / len(z)) * (z - z_orig) / zs])
        return r_cc
    sol = least_squares(lresid, z0, method="trf", x_scale=zs, bounds=(lo, hi), max_nfev=max_nfev, diff_step=1e-6)
    set_vec(sol.x)
    if l.cc_of(F_cur) <= cc0:
        l.restore(snap0)

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
m0 = np.mean(list(own0.values()))
best_snaps = [l.snapshot() for l in L]
best = (m0, th0.copy(), dlod0, own0, best_snaps); th_b = th0.copy(); step = a.step; leg_step = a.leg_step; best_v = v0
per = max(4, a.iters // a.rounds); anchor = [l.F.copy() for l in L]; orig_nh = len(L[0].harm)
if a.homotopy:
    # Coarse-to-fine frequency continuation (FWI homotopy): align slow/base windings first so fast harmonics don't cycle-skip
    max_k_all = max(max((abs(k) for k in l.ks), default=0.0) for l in L)
    stages = [(km, min(a.lam, 0.25), st) for km, st in ((0.35, 3.0), (0.90, 2.5), (2.50, 2.0)) if km < max_k_all - 1e-6]
    if stages:
        th_h = best[1].copy(); t_h = time.time()
        for idx_s, (km, lam_h, st_h) in enumerate(stages, 1):
            def resid_h(x):
                th = full(x); r = []
                for l in L:
                    y = l.model_of(l.F + l.J @ (th - th_b), kmax=km)[l.sc]; d = l.data[l.sc]
                    r.append(((y - y.mean()) / max(float(y.std()), 1e-12) - (d - d.mean()) / d.std()) / np.sqrt(2 * len(d) * len(L)))
                return np.concatenate(r + [np.sqrt(lam_h) * (PHI @ (th - thr)) / NP])
            x0_h = th_h[FREE]; lo_h, hi_h = x0_h - st_h * XS[FREE], x0_h + st_h * XS[FREE]
            sol_h = least_squares(resid_h, x0_h, method="trf", x_scale=XS[FREE], bounds=(lo_h, hi_h), max_nfev=max(12, per * 2), diff_step=1e-6)
            th_h = full(sol_h.x)
            print(f"homotopy stage {idx_s} (|k|<={km:g}, lam={lam_h:g}): surrogate full-spectrum mean {np.mean(ccs(th_h)):.4f}, dLOD {dlod_cc(th_h):.4f} [{sol_h.nfev} evals]")
        R0 = os.path.join(W, "homotopy"); copy_run(R0, lp_of(th_h)); test_only(R0, threads=max(2, len(SUF))); own_h, dlod_h, _ = run_scores(R0)
        m_h = np.mean(list(own_h.values())) if len(own_h) == len(SUF) else float("nan")
        v_h = vscore(R0) if (a.accept == "validate" and len(own_h) == len(SUF)) else float("nan")
        ok_h = len(own_h) == len(SUF) and m_h > best[0] and (a.accept == "train" or v_h >= best_v - a.vtol)
        print(f"homotopy solver: {m_h:.4f} vs baseline {best[0]:.4f} {'accepted' if ok_h else 'rejected'}, dLOD {dlod_h:.4f} [{time.time()-t_h:.0f} s]")
        if ok_h:
            best = (m_h, th_h.copy(), dlod_h, own_h, best[4]); best_v = v_h
            anchor = [np.loadtxt(os.path.join(R0, f"lte_results{suf}.csv"), delimiter=",")[:, 3] for suf in SUF]
for rnd in range(1, a.rounds + 1):
    for l, F_, sn in zip(L, anchor, best[4]): l.F = F_.copy(); l.restore(sn)
    th_b = best[1].copy(); t2 = time.time(); c0 = ccs(th_b)
    # Optional harmonic expansion (round 1 only) & 1-swap screening
    if a.harm_screen or (rnd == 1 and a.harm_expand > 0):
        exp_n = a.harm_expand if rnd == 1 else 0
        # Keep NH identical across legs when expanding: expand each leg by exact same count if possible
        for idx_l, l in enumerate(L):
            screen_harmonics(l, idx_l, l.F, MAXH, a.harm_max_swaps if a.harm_screen else 0, a.harm_min_gain, expand_n=exp_n)
        if exp_n > 0:
            min_nh = min(len(l.harm) for l in L)
            for l in L:
                if len(l.harm) > min_nh: l.harm = l.harm[:min_nh]; l.sync_ks()
    # Block A: shared lpap trust-region GN
    x0 = th_b[FREE]; lo, hi = x0 - step * XS[FREE], x0 + step * XS[FREE]
    sol = least_squares(resid, x0, method="trf", x_scale=XS[FREE], bounds=(lo, hi), max_nfev=per, diff_step=1e-6); th1 = full(sol.x)
    # Block B: per-leg continuous GN (independent per leg, holding F1 = F + J @ dth fixed)
    if a.per_leg:
        for l in L:
            F1 = l.F + l.J @ (th1 - th_b)
            refine_per_leg(l, F1, leg_step, a.leg_lam, max_nfev=max(8, per))
    c1 = ccs(th1); snaps1 = [l.snapshot() for l in L]
    patches1 = [l.to_patch(include_bessel=HAS_PER_LEG_BESSEL) for l in L] if DO_LEG else None
    new_nh = len(L[0].harm) if (a.harm_expand > 0 and len(L[0].harm) != orig_nh) else None
    R = os.path.join(W, f"round{rnd}"); copy_run(R, lp_of(th1), leg_patches=patches1, new_nh=new_nh); test_only(R, threads=max(2, len(SUF))); own1, dlod1, _ = run_scores(R)
    m1 = np.mean(list(own1.values())) if len(own1) == len(SUF) else float("nan")
    v1 = vscore(R) if (a.accept == "validate" and len(own1) == len(SUF)) else float("nan")
    ok = len(own1) == len(SUF) and m1 > best[0] and (a.accept == "train" or v1 >= best_v - a.vtol)
    print(f"round {rnd} (step {step:.3g}" + (f", leg_step {leg_step:.3g}" if a.per_leg else "") + f"): surrogate {np.mean(c0):.4f} -> {np.mean(c1):.4f} | solver {m1:.4f} vs best {best[0]:.4f} "
          + (f" | validation {v1:.4f} vs {best_v:.4f}" if a.accept == "validate" else "") +
          f" {'accepted' if ok else 'rejected'}, dLOD {dlod1:.4f} [{sol.nfev} evaluations, {time.time()-t2:.0f} s]")
    if ok:
        best = (m1, th1.copy(), dlod1, own1, snaps1); best_v = v1; step *= 1.5; leg_step = min(leg_step * 1.5, 0.25)
        anchor = [np.loadtxt(os.path.join(R, f"lte_results{suf}.csv"), delimiter=",")[:, 3] for suf in SUF]
    else:
        step *= 0.5; leg_step *= 0.5
    if step < 0.01: break
m1, th1, dlod1, own1, snaps_best = best; LP1 = lp_of(th1)
for l, sn in zip(L, snaps_best): l.restore(sn)
final_nh = len(L[0].harm) if (a.harm_expand > 0 and len(L[0].harm) != orig_nh) else None
print(f"solver, training per leg: " + " ".join(f"{own0[k]:.4f}->{own1[k]:.4f}" for k in sorted(own0)))
print(f"solver, training mean {m0:.4f} -> {m1:.4f}; dLOD {dlod0:.4f} -> {dlod1:.4f}  [total {time.time()-t0:.0f} s]")
json.dump({"lpap": LP1.tolist(), "legs": [l.to_patch(include_bessel=HAS_PER_LEG_BESSEL) for l in L] if DO_LEG else None}, open(os.path.join(W, "refined_lpap.json"), "w"))
if a.apply:
    if m1 > m0:
        shutil.copy(os.path.join(D, "lt.exe.p"), os.path.join(D, "lt.exe.p.pre_refine"))
        p = json.load(open(os.path.join(D, "lt.exe.p"))); p["lpap"] = LP1.tolist()
        if DO_LEG:
            patches = [l.to_patch(include_bessel=HAS_PER_LEG_BESSEL) for l in L]
            p.update(patches[0])
            for suf, patch in zip(SUF, patches):
                fp = os.path.join(D, f"lt.exe.p{suf}")
                if suf and os.path.exists(fp):
                    shutil.copy(fp, fp + ".pre_refine")
                    pl = json.load(open(fp)); pl.update(patch); json.dump(pl, open(fp, "w"), indent=2)
        json.dump(p, open(os.path.join(D, "lt.exe.p"), "w"), indent=2)
        if final_nh is not None:
            rp = os.path.join(D, "lt.exe.resp"); shutil.copy(rp, rp + ".pre_refine"); _update_resp_nh(rp, final_nh)
            print(f"applied: updated NH {orig_nh} -> {final_nh} in {rp}")
        print(f"applied: {D}/lt.exe.p" + (" and per-leg files" if DO_LEG and len(SUF) > 1 else "") + " updated (backups *.pre_refine)")
    else:
        print("not applied: the solver's training score did not improve")
else:
    print(f"refined lpap: {W}/refined_lpap.json (use --apply to write it into lt.exe.p)")
