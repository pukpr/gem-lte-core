#!/usr/bin/env python3
"""Tests of the two-leg worker groups in lt.exe (LEGS 2, lt.exe.resp.1 / lt.exe.resp.2), 2026-10-12.
 T0  no legs: new binary against the pre-change binary, 5-day and monthly (must be identical)
 T1  legs with identical settings on one series: both leg files equal the single run
 T2  one series, different interval per leg (leg 1 trains 1992-2005, leg 2 trains outside it): each leg file equals
     the corresponding single run; the composite printed is the mean of the two single scores
 T3  5-day leg + monthly leg (different SAMPLING and index file): each leg file equals its single run
 T4  a short paired search on one series with two intervals: ends by itself, composite rises, the saved parameter
     file reproduces both saved leg files
All from soi_daily/rr2/start.p (the t_daily fit)."""
import json, os, re, shutil, subprocess, sys, time, numpy as np
from pathlib import Path
E = Path("/home/paul/eval/gem-lte-core/experiments"); H = Path(__file__).resolve().parent; R = E / "Oct2026_generator/soi_daily/rr2"; NEW, OLD = E / "Feb2026/enso_opt", E / "Feb2026/enso_opt.pre-legs"
ENV0 = {k: v for k, v in re.findall(r"^export (\w+)=(\S+)", (R / "daily/lte_run.sh").read_text(), re.M)}
for k in ("CLIMATE_INDEX", "SAMPLING", "EXCLUDE", "TRAIN_START", "TRAIN_END", "TEST_ONLY", "TIMEOUT", "NUMBER_OF_PROCESSORS"): ENV0.pop(k, None)
RESP0 = re.sub(r"^-?(SAMPLING)\s.*\n", "", (R / "daily/lte_run.sh").read_text() * 0 + (R / "daily/lt.exe.resp").read_text(), flags=re.M)
D5 = dict(CLIMATE_INDEX="soi_lag_5.dat", SAMPLING="73.0"); DM = dict(CLIMATE_INDEX="soi_lag_m.dat", SAMPLING="12.0"); IN = dict(TRAIN_START="1992", TRAIN_END="2005", EXCLUDE="FALSE"); OUT = dict(TRAIN_START="1992", TRAIN_END="2005", EXCLUDE="TRUE")
def resp_text(d): return "".join(f"{k:30s} {v}\n" for k, v in d.items())
def run(name, binary, common, legs=None, test=True, timeout=60, threads=2, pfile=None):
    d = H / name; shutil.rmtree(d, ignore_errors=True); d.mkdir()
    for f in ("soi_lag_5.dat",): shutil.copy(R / "daily" / f, d / f)
    shutil.copy(R / "monthly/soi_lag_m.dat", d / "soi_lag_m.dat"); shutil.copy(R / "daily/dlod_ref.dat", d / "dlod_ref.dat"); shutil.copy(pfile or (R / "start.p"), d / "lt.exe.p"); os.symlink(binary, d / "lt.exe"); os.symlink(E / "Feb2026/dlod3.dat", H / "dlod3.dat") if not (H / "dlod3.dat").exists() else None
    (d / "lt.exe.resp").write_text(RESP0 + resp_text(common) + (f"LEGS                           {len(legs)}\n" if legs else ""))
    for i, L in enumerate(legs or [], 1): (d / f"lt.exe.resp.{i}").write_text(resp_text(L))
    env = {**os.environ, **ENV0, "TEST_ONLY": "TRUE" if test else "FALSE", "TIMEOUT": str(timeout), "NUMBER_OF_PROCESSORS": str(threads)}; t0 = time.time()
    try: p = subprocess.run("ulimit -s 65536; ./lt.exe", shell=True, cwd=d, env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=timeout + 240); out = p.stdout + p.stderr
    except subprocess.TimeoutExpired as e: out = (e.stdout or b"").decode(errors="ignore") + "\n*** DID NOT EXIT ***"
    (d / "run.out").write_text(out); return d, out, time.time() - t0
def res(d, tag=""):
    a = np.loadtxt(d / f"lte_results{tag}.csv", delimiter=","); return a
def cmp(a, b): return "IDENTICAL" if a.shape == b.shape and np.array_equal(a, b) else (f"max |diff| {np.abs(a - b).max():.2e}" if a.shape == b.shape else f"SHAPES {a.shape} {b.shape}")
cc = lambda x, y: float(np.corrcoef(x, y)[0, 1])
def scores(out): return [tuple(float(x) for x in m) for m in re.findall(r"^DTW\s+([-\d.]+)\s+([-\d.]+)\s+\d+", out, re.M)]
def status(out): return [float(m) for m in re.findall(r"^COMPOSITE\s+([-\d.]+)", out, re.M)] or [float(m) for m in re.findall(r"^Status:\s+\d+\s+\d+\s+([-\d.]+)", out, re.M)]
which = sys.argv[1:] or ["T0", "T1", "T2", "T3", "T4"]; S = {}
def single(key, common):
    if key not in S: d, out, _ = run("s_" + key, OLD, common); S[key] = (res(d), scores(out)[-1] if scores(out) else None)
    return S[key]
if "T0" in which:
    for key, common in (("5in", {**D5, **IN}), ("m_out", {**DM, **OUT})):
        a, sc = single(key, common); d, out, _ = run("t0_" + key, NEW, common); print(f"T0 {key:6s} new binary without legs vs old binary: {cmp(res(d), a)}; scores old {sc}, new {scores(out)[-1]}", flush=True)
if "T1" in which:
    a, sc = single("5in", {**D5, **IN}); d, out, dt = run("t1", NEW, {**D5, **IN}, legs=[{}, {}]); print(f"T1 identical legs: leg 1 vs single {cmp(res(d, '.1'), a)}; leg 2 vs single {cmp(res(d, '.2'), a)}; single score {sc}; composite shown {status(out)[-1:]}; exited in {dt:.0f} s", flush=True)
if "T2" in which:
    a, sa = single("5in", {**D5, **IN}); b, sb = single("5out", {**D5, **OUT}); d, out, dt = run("t2", NEW, D5, legs=[IN, OUT])
    print(f"T2 one series, two intervals: leg 1 (trains 1992-2005) vs single {cmp(res(d, '.1'), a)}; leg 2 (trains outside) vs single {cmp(res(d, '.2'), b)}")
    print(f"   single training scores {sa[0]:.5f} (inside) and {sb[1]:.5f} (outside), mean {0.5 * (sa[0] + sb[1]):.5f}; composite shown {status(out)[-1:]}; exited in {dt:.0f} s", flush=True)
if "T3" in which:
    a, sa = single("5in", {**D5, **IN}); b, sb = single("m_out", {**DM, **OUT}); d, out, dt = run("t3", NEW, {}, legs=[{**D5, **IN}, {**DM, **OUT}])
    print(f"T3 5-day leg + monthly leg: leg 1 vs single 5-day {cmp(res(d, '.1'), a)}; leg 2 vs single monthly {cmp(res(d, '.2'), b)}")
    print(f"   single training scores {sa[0]:.5f} (inside) and {sb[1]:.5f} (outside), mean {0.5 * (sa[0] + sb[1]):.5f}; composite shown {status(out)[-1:]}; exited in {dt:.0f} s", flush=True)
if "T5" in which:
    LATE = dict(TRAIN_START="2005", TRAIN_END="2026", EXCLUDE="FALSE"); a, sa = single("5in", {**D5, **IN}); b, sb = single("5out", {**D5, **OUT}); c, sc = single("5late", {**D5, **LATE}); d, out, dt = run("t5", NEW, D5, legs=[IN, OUT, LATE], threads=3)
    print(f"T5 three legs, one series, three intervals: leg 1 {cmp(res(d, '.1'), a)}; leg 2 {cmp(res(d, '.2'), b)}; leg 3 {cmp(res(d, '.3'), c)}")
    print(f"   single training scores {sa[0]:.5f}, {sb[1]:.5f}, {sc[0]:.5f}, mean {(sa[0] + sb[1] + sc[0]) / 3:.5f}; composite shown {status(out)[-1:]}; exited in {dt:.0f} s", flush=True)
    d, out, dt = run("t5s", NEW, D5, legs=[IN, OUT, LATE], test=False, timeout=60, threads=9); st = status(out); print(f"   three-leg search, 9 tasks (3 groups), 60 s: exited in {dt:.0f} s{' *** HUNG ***' if 'DID NOT EXIT' in out else ''}; composite at save {st[-1:]}; leg files written: {[f.name for f in sorted(d.glob('lte_results.*.csv'))]}; errors {len(re.findall('Solution err', out))}")
if "T6" in which:
    import json as _j
    a, sa = single("5in", {**D5, **IN}); d, out, dt = run("t6", NEW, {**D5, **IN, "LEG_COMMON": "30"}, legs=[{}, {}]); r1, r2 = res(d, ".1"), res(d, ".2")
    print(f"T6 common response, two identical legs: leader vs single {cmp(r1, a)}; follower vs leader: max |model diff| {np.abs(r1[:, 1] - r2[:, 1]).max():.1e} (model sd {r1[:, 1].std():.2f}); exited in {dt:.0f} s")
    d, out, dt = run("t7", NEW, {**D5, "LEG_COMMON": "30"}, legs=[IN, OUT]); w1, w2 = (_j.load(open(d / f"lt.exe.windings.{i}.json"))["k_amp_phase"] for i in (1, 2)); kb = abs(w1[1][0]); hi = [i for i, x in enumerate(w1) if i > 1 and round(abs(x[0]) / kb) > 30]; lo = [i for i in range(len(w1)) if i not in hi]
    r2 = res(d, ".2"); t_ = r2[:, 0]; w_ = (t_ >= 1992) & (t_ < 2005); b_, sb = single("5out", {**D5, **OUT})
    print(f"T7 common response, leg 1 trains 1992-2005, leg 2 outside: harmonics above 30 identical in both legs: {all(abs(w1[i][1] - w2[i][1]) < 1e-12 and abs(w1[i][2] - w2[i][2]) < 1e-12 for i in hi)} ({len(hi)} of them); lower terms differ: {any(abs(w1[i][1] - w2[i][1]) > 1e-6 for i in lo)}")
    print(f"   leg 2 outside-window CC with the leader's fast part imposed {cc(r2[~w_, 1], r2[~w_, 2]):+.3f}; with its own fast part (no LEG_COMMON) {cc(b_[~w_, 1], b_[~w_, 2]):+.3f}; composite {status(out)[-1:]}")
    d, out, dt = run("t8", NEW, {**D5, "LEG_COMMON": "30"}, legs=[IN, OUT], test=False, timeout=60, threads=8); st = status(out); print(f"   60 s search with LEG_COMMON, 4 pairs: exited in {dt:.0f} s{' *** HUNG ***' if 'DID NOT EXIT' in out else ''}; composite at save {st[-1:]}; errors {len(re.findall('Solution err', out))}")
if "T9" in which:
    FAST = {"LEG_SCORE": "FAST"}; d, out, dt = run("t9", NEW, {**D5, **IN, "LEG_COMMON": "30", "METRIC": "CC"}, legs=[{}, FAST]); r1 = res(d, ".1"); r2 = res(d, ".2"); st = status(out)
    print(f"T9 fast-part score, two identical legs, METRIC CC from the environment override: leg files equal {np.abs(r1[:, 1] - r2[:, 1]).max() < 1e-9}; composite {st[-1:]} (leader's own score {scores(out)[-1][0] if scores(out) else '?'})")
    d, out, dt = run("t9b", NEW, {**D5, "LEG_COMMON": "30"}, legs=[IN, {**OUT, **FAST}]); print(f"   leg 1 trains 1992-2005, leg 2 outside scored on the fast part: composite {status(out)[-1:]}; exited in {dt:.0f} s; errors {len(re.findall('Solution err', out))}")
    d, out, dt = run("t9c", NEW, {**D5, "LEG_COMMON": "30"}, legs=[IN, {**OUT, **FAST}], test=False, timeout=60, threads=8); print(f"   60 s search: exited in {dt:.0f} s{' *** HUNG ***' if 'DID NOT EXIT' in out else ''}; composite at save {status(out)[-1:]}; errors {len(re.findall('Solution err', out))}")
if "T4" in which:
    d, out, dt = run("t4", NEW, D5, legs=[IN, OUT], test=False, timeout=90, threads=8); st = status(out); a1, a2 = res(d, ".1"), res(d, ".2")
    print(f"T4 paired search, 8 tasks (4 pairs), 90 s: exited in {dt:.0f} s{' *** HUNG ***' if 'DID NOT EXIT' in out else ''}; composite {st[0]:.5f} -> {st[-1]:.5f} over {len(st)} reports; errors in output: {len(re.findall('Solution err', out))}")
    d2, out2, _ = run("t4_check", NEW, D5, legs=[IN, OUT], pfile=d / "lt.exe.p"); print(f"   saved parameter file re-run: leg 1 {cmp(res(d2, '.1'), a1)}; leg 2 {cmp(res(d2, '.2'), a2)}; composite on re-run {status(out2)[-1:]}")
    t = a1[:, 0]; w = (t >= 1992) & (t < 2005); print(f"   CC of the saved fit: leg 1 inside 1992-2005 {cc(a1[w, 1], a1[w, 2]):+.3f}; leg 2 outside {cc(a2[~w, 1], a2[~w, 2]):+.3f} (start: +0.720 with outside fitted separately)")
