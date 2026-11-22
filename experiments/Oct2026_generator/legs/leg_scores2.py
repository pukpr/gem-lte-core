#!/usr/bin/env python3
"""Two legs on two indices with 2000-2005 held out: train and held-out correlations per leg, dLOD-free summary of the
shared windings, the per-leg amplitudes by harmonic, and a plot. usage: leg_scores2.py [dir] [name1] [name2]
-post (or --post): also draw each leg's post-Bessel manifold on the bottom panel, dashed in the leg's colour. It is
computed from column 4, which must then hold the pre-Bessel manifold (MANIFOLD_PRE TRUE), with the leg's own impA, impB,
besk, offs, bg, bg2 (lt.exe.p.<leg>, else lt.exe.p), as the solver's Bessel step does. LAYER and FRICTION runs are not covered.
Regions follow the run's EXCLUDE (environment, else the last uncommented export in DIR/lte_run.sh; default TRUE):
  EXCLUDE TRUE:  the --held interval (else TRAIN_START..TRAIN_END, else 2000-2005) is held out, the rest is training.
  EXCLUDE FALSE: the inverse. TRAIN_START..TRAIN_END (else the --held interval) is training, the rest is held out;
                 with VALIDATE TRUE and VALIDATE_START, the part from VALIDATE_START on is reported as validation."""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H0, H1 = 2000.0, 2005.0
X0 = None; OUTPNG = "legs.png"   # X0 None: start at the decade of the earliest leg's first sample
if "--from" in sys.argv: i = sys.argv.index("--from"); X0 = float(sys.argv[i + 1]); del sys.argv[i:i + 2]
if "--out" in sys.argv: i = sys.argv.index("--out"); OUTPNG = sys.argv[i + 1]; del sys.argv[i:i + 2]
POST = any(x in sys.argv for x in ("-post", "--post")); sys.argv = [x for x in sys.argv if x not in ("-post", "--post")]
HELD = None
if "--held" in sys.argv: i = sys.argv.index("--held"); HELD = (float(sys.argv[i + 1]), float(sys.argv[i + 2])); del sys.argv[i:i + 3]
d = sys.argv[1] if len(sys.argv) > 1 else "."
def run_setting(key):   # the environment (set when lte_run.sh runs this), else the last uncommented export of KEY in lte_run.sh
    import os, re
    if key in os.environ: return os.environ[key]
    try: v = re.findall(r"^\s*export\s+[^#\n]*?\b" + key + r"=([^\s#]+)", open(f"{d}/lte_run.sh").read(), re.M)
    except OSError: v = []
    return v[-1] if v else None
def run_year(key):
    v = run_setting(key)
    try: return float(v)
    except (TypeError, ValueError): return None
EXC = (run_setting("EXCLUDE") or "TRUE").upper() not in ("FALSE", "0", "NO")
TS, TE = run_year("TRAIN_START"), run_year("TRAIN_END")
VS = run_year("VALIDATE_START") if (run_setting("VALIDATE") or "").upper() == "TRUE" else None
if EXC: H0, H1 = HELD or ((TS, TE) if TS is not None and TE is not None else (H0, H1))
else:
    H0, H1 = (TS if TS is not None else -1e9, TE if TE is not None else 1e9) if (TS is not None or TE is not None) else (HELD or (H0, H1))
    if HELD and HELD != (H0, H1): print(f"EXCLUDE FALSE: training is TRAIN_START..TRAIN_END ({H0:g}-{H1:g}); --held {HELD[0]:g} {HELD[1]:g} is not used")
def regions(t):   # train, held-out and validation masks, and the held-out label
    inside = (t >= H0) & (t < H1)
    if EXC: return ~inside, inside, np.zeros_like(inside), f"{H0:g}-{H1:g}"
    va = (t >= VS) if VS is not None else np.zeros_like(inside); ho = ~inside & ~va
    lab = (f"{H1:g}-{VS:g}" if VS is not None else f"{H1:g}-" if H0 <= t[0] else f"before {H0:g}" if H1 >= t[-1] else f"outside {H0:g}-{H1:g}")
    return inside, ho, va, lab
print(f"regions: EXCLUDE {'TRUE' if EXC else 'FALSE'}, " + (f"held out {H0:g}-{H1:g}, training the rest" if EXC else f"training {H0:g}-{H1:g}, held out the rest" + (f", validation from {VS:g}" if VS is not None else "")))
names = sys.argv[2:] if len(sys.argv) > 3 else ["leg 1", "leg 2"]; NL = len(names); cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v: (v - v.mean()) / v.std()
fig, ax = plt.subplots(NL + 1, 1, figsize=(17, 4 * (NL + 1)), facecolor="#fcfcfb"); F = []   # last panel: every leg's manifold together
if X0 is None:
    def _start(k):
        a = np.loadtxt(f"{d}/lte_results.{k}.csv", delimiter=","); nz = a[a[:, 2] != 0, 0]; return float(nz[0]) if len(nz) else float(a[0, 0])
    X0 = 10.0 * np.floor(min(_start(k) for k in range(1, NL + 1)) / 10.0)
def shade(A, t0, t1):   # held out: grey; validation: lighter grey; training: unshaded
    if EXC: A.axvspan(H0, H1, color="#e4e3df", alpha=0.7, lw=0); return
    e = VS if VS is not None else 1e9
    if H0 > t0: A.axvspan(t0, H0, color="#e4e3df", alpha=0.7, lw=0)
    if H1 < min(e, t1): A.axvspan(H1, min(e, t1), color="#e4e3df", alpha=0.7, lw=0)
    if VS is not None and VS < t1: A.axvspan(VS, t1, color="#efeeeb", alpha=0.7, lw=0)
for i, name in enumerate(names):
    a = np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=","); a = a[np.abs(a[:, 2]) > 0] if (a[:60, 2] == 0).all() else a; t = a[:, 0]; tr, g, va, hl = regions(t); g2 = (t >= 2000) & (t < 2005); F.append(a[:, [0, 3]])
    dq = lambda v: v - np.polyval(np.polyfit(t[tr], v[tr], 2), t); trend_only = cc(np.polyval(np.polyfit(t[tr], a[tr, 2], 2), t)[tr], a[tr, 2])
    vtxt = f" | validate {VS:g}- {cc(a[va, 1], a[va, 2]):+.3f}" if va.sum() > 2 else ""
    print(f"{name}: {t[0]:.1f}-{t[-1]:.1f} | train {cc(a[tr, 1], a[tr, 2]):+.3f} (trend alone {trend_only:+.3f}; with quadratic trend removed from both {cc(dq(a[:, 1])[tr], dq(a[:, 2])[tr]):+.3f}) | held out {hl} {cc(a[g, 1], a[g, 2]):+.3f}{vtxt} | 2000-05 {cc(a[g2, 1], a[g2, 2]):+.3f} | by 36-yr era: " + " ".join(f"{cc(a[(t >= y) & (t < y + 36), 1], a[(t >= y) & (t < y + 36), 2]):+.2f}" for y in range(1880, 2020, 36) if ((t >= y) & (t < y + 36)).sum() > 100))
    A = ax[i]; A.set_facecolor("#fcfcfb"); shade(A, t[0], t[-1]); A.plot(t, z(a[:, 2]), color="#2a78d6", lw=1.1, label=f"{name} (standardized)"); A.plot(t, z(a[:, 1]), color="#eb6834", lw=1.1, label="model (standardized)")
    A.set_title(f"{name}: train {cc(a[tr, 1], a[tr, 2]):+.3f}, held out {hl} {cc(a[g, 1], a[g, 2]):+.3f}" + (f", validate {VS:g}- {cc(a[va, 1], a[va, 2]):+.3f}" if va.sum() > 2 else "") + ("" if EXC else f"   (EXCLUDE FALSE: training {H0:g}-{H1:g})"), loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=2); A.set_xlim(X0, max(2023, float(t[-1])))
    for s in ("top", "right"): A.spines[s].set_visible(False)
print("manifold of each leg against leg 1: " + ", ".join(f"{cc(np.interp(F[i][:, 0], F[0][:, 0], F[0][:, 1])[F[i][:, 0] <= F[0][-1, 0]], F[i][F[i][:, 0] <= F[0][-1, 0], 1]):.5f}" for i in range(1, NL)))
def bessel_post(k, m):   # the solver's Bessel step (gem-lte-primitives-solution.adb), with leg k's own parameters
    p = json.load(open(f"{d}/lt.exe.p"));
    try: p.update(json.load(open(f"{d}/lt.exe.p.{k}")))
    except OSError: pass
    g = lambda key, dv=0.0: float(p.get(key, dv)); eS, eC, bg = g("impA"), g("impB"), g("bg"); nm = len(p.get("ltep", [0, 0, 0]))
    kk, e2s, e2c = (g("besk") * (1.0 - bg), g("offs"), g("bg2", bg)) if nm == 1 else (g("besk"), g("offs"), bg)
    w = 2 * np.pi * kk * m; return m + eS * np.sin(w) + eC * np.cos(w) + e2s * eS * np.sin(2 * w) + e2c * eC * np.cos(2 * w)
if POST:   # column 4 must be pre-Bessel: MANIFOLD_PRE on in the run script or resp file, or (several legs) identical across legs
    import re, glob
    said = any(re.search(r"^\s*(export\s+.*\bMANIFOLD_PRE=(TRUE|true|1)\b|MANIFOLD_PRE\s+(TRUE|true))", open(f).read(), re.M) for f in glob.glob(f"{d}/lte_run.sh") + glob.glob(f"{d}/lt.exe.resp"))
    same = NL > 1 and all(cc(np.diff(np.interp(F[i][F[i][:, 0] <= F[0][-1, 0], 0], F[0][:, 0], F[0][:, 1])), np.diff(F[i][F[i][:, 0] <= F[0][-1, 0], 1])) > 0.9999 for i in range(1, NL))
    if not (said or same): print("-post: column 4 already looks post-Bessel (MANIFOLD_PRE not set, legs differ), so nothing extra is drawn"); POST = False
    else: FP = [np.column_stack([f[:, 0], bessel_post(i + 1, f[:, 1])]) for i, f in enumerate(F)]; print("post-Bessel manifold of each leg against leg 1: " + ", ".join(f"{cc(np.interp(FP[i][:, 0], FP[0][:, 0], FP[0][:, 1])[FP[i][:, 0] <= FP[0][-1, 0]], FP[i][FP[i][:, 0] <= FP[0][-1, 0], 1]):.5f}" for i in range(1, NL)))
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]   # categorical order, fixed per leg
A = ax[NL]; A.set_facecolor("#fcfcfb"); x1 = max(2023, max(float(f[-1, 0]) for f in F)); shade(A, min(float(f[0, 0]) for f in F), x1)
for i in range(NL):   # column 4 of lte_results.<leg>.csv, raw (same units on every leg), so offsets and amplitudes show as they are
    m = F[i][:, 0] <= F[0][-1, 0]; r = 1.0 if i == 0 else cc(np.interp(F[i][m, 0], F[0][:, 0], F[0][:, 1]), F[i][m, 1])
    A.plot(F[i][:, 0], F[i][:, 1], color=SERIES[i % len(SERIES)], lw=1.1, alpha=0.9, label=names[i] if i == 0 else f"{names[i]} (CC vs {names[0]} {r:.4f})")
for i in range(NL if POST else 0):   # dashed on top of every solid line (the pre-Bessel lines coincide), same colour, kept out of the legend
    A.plot(FP[i][:, 0], FP[i][:, 1], color=SERIES[i % len(SERIES)], lw=1.0, alpha=0.85, ls=(0, (4, 2.5)), label="_nolegend_")
if POST: A.plot([], [], color="#52514d", lw=1.0, ls=(0, (4, 2.5)), label="dashed: post-Bessel, same colours")   # one legend key for all the dashed lines
lo, hi = A.get_ylim(); A.set_ylim(lo, hi + 0.18 * (hi - lo))   # headroom for the legend
if NL <= 4:   # end labels, pushed apart where the lines finish close together, each with a leader to its line's end
    ys = [float(f[-1, 1]) for f in F]; gap = 0.06 * (hi - lo); y = {}
    for i in sorted(range(NL), key=lambda j: ys[j]): y[i] = ys[i] if not y else max(ys[i], max(y.values()) + gap)
    for i in range(NL): A.annotate(names[i], (F[i][-1, 0], ys[i]), xytext=(F[i][-1, 0] + 1.5, y[i]), textcoords="data", va="center", fontsize=9, color="#52514d", arrowprops=dict(arrowstyle="-", color="#a3a29c", lw=0.6, shrinkA=0, shrinkB=1))
A.set_title(f"manifolds of all {NL} legs (lte_results column 4, not standardized" + ("; solid pre-Bessel, dashed post-Bessel)" if POST else ")"), loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=min(NL, 4), loc="upper left"); A.set_xlim(X0, x1 + (6 if NL <= 4 else 0))
for s_ in ("top", "right"): A.spines[s_].set_visible(False)
W = [json.load(open(f"{d}/lt.exe.windings.{i}.json")) for i in range(1, NL + 1)]; K = [k for k, _, _ in W[0]["k_amp_phase"]]; base = K[1] if len(K) > 1 else K[0]   # no harmonics: slow winding only
KL = [[k for k, _, _ in w["k_amp_phase"]] for w in W]
if any(len(k) != len(K) or not np.allclose(k, K, rtol=0, atol=1e-12) for k in KL[1:]):   # LEG_WINDINGS: each leg its own ltep/harm
    print("per-leg windings (LEG_WINDINGS): wavenumber amplitude (phase in degrees)")
    for i in range(NL): print(f"   {names[i]:>8s}  " + "   ".join(f"{k:+.4f} {a_:6.3f} ({np.degrees(p_):+5.0f})" for k, a_, p_ in W[i]["k_amp_phase"]))
else: print(f"shared windings: slow {K[0]:+.4f}, base {base:.4f}; amplitude (phase in degrees) by harmonic:")
for j in (np.argsort(np.abs(K)) if not any(len(k) != len(K) or not np.allclose(k, K, rtol=0, atol=1e-12) for k in KL[1:]) else []): print(f"   {('slow' if j == 0 else f'{K[j] / K[0]:.0f}x'):>5s}  " + "   ".join((f"{names[i]} {W[i]['k_amp_phase'][j][1]:6.3f} ({np.degrees(W[i]['k_amp_phase'][j][2]):+5.0f})" if j < len(W[i]['k_amp_phase']) else f"{names[i]}      -       ") for i in range(NL)))   # a leg with fewer windings (shorter per-leg NH) shows -
fig.savefig(f"{d}/{OUTPNG}", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print(f"plot: {d}/{OUTPNG}")
