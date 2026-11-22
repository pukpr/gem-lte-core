#!/usr/bin/env python3
"""LTE modulation of every leg on one graph: the model (lte_results column 2) against the manifold (column 4), folded
into one Brillouin zone. Started from harmonic_spectrum.py. usage: lte_modulation.py DIR [NAME ...] [--fit H] [--offset Y] [--out FILE]

Each NAME picks a leg: the one whose lt.exe.resp.<leg> CLIMATE_INDEX file it names, case-insensitive and without .dat
(nino4 -> nino4.dat), or a leg number. No names: every leg. Options may go anywhere on the line.

For each leg <i> (lte_results.<i>.csv, lt.exe.windings.<i>.json, lt.exe.p.<i>):
  * IR: the model is IR-decompensated with the exact recursive inverse, u(n) = model(n) + IR*u(n-12) (as in
    harmonic_spectrum.py and Feb2026/plot.py). Skipped when |IR| >= 1, where the inverse diverges.
  * Manifold: the windings act on the post-Bessel manifold. If column 4 is pre-Bessel (MANIFOLD_PRE on in the run
    script or lt.exe.resp), the leg's own Bessel step is applied first, with its impA, impB, besk, offs, bg and bg2
    from the windings JSON (the solver's formula; LAYER and FRICTION runs are not covered).
  * Terms that are not a function of the folded manifold are removed by least squares: 1, F, t, (t - accel_ref)^2,
    and sin/cos at 1/yr, 2/yr and the CW alias (frac(year / (draconic/2)) ~ 0.844/yr).
  * x: the manifold folded by the leg's base winding k (its last ltep, the BZONE convention), frac(k F), in [0, 1).
  * y: what is left, divided by its standard deviation, so legs in different units share the axis.
  * The points are sorted by x and joined; lines are offset vertically by OFFSET (default 0 = overlaid) with --offset.
    Windings that are not integer multiples of k (the slow ltep ones) do not fold onto one curve, so the joined
    line also shows a band around it. --fit H also draws, bold, the least-squares fit of sin/cos(2 pi h x),
    h = 1..H, which is the part that does fold onto one curve.

--limbs: a test of whether the modulation is a rule (the same F gives the same response) rather than a description.
  The manifold rises to a maximum and falls back, so the falling limb revisits manifold values first seen on the
  rising one. For each leg, ONE harmonic h of k (amplitude and phase: 2 numbers) is fitted to the DATA (column 3,
  IR-matched as in harmonic_spectrum.py: data(n) + IR*data(n-12)) on one limb, together with terms that do not involve
  F (1, t, t^2, sin/cos at 1/yr, 2/yr and the CW alias), and scored on the other limb: the CC between the predicted
  sinusoid and the data there, after the same non-F terms are removed from it. Both directions are reported.
    --turn YEAR   where the limbs meet (default: the date of the largest column-4 value of each leg)
    --harm H      use harmonic H of each leg's own k (default: the best of 1..--hmax, 9, judged on both limbs,
                  so the score is mildly selected; give --harm or --kh to avoid that)
    --kh K        use the wavenumber K itself on every leg, whatever its k (e.g. 1.155, shared by NINO4, PDO, IODE)
    --null N      reference: the same test with N random wavenumbers in [0.05, 3] (default 200; 0 = none),
                  reported as mean, 95th percentile and maximum (--seed S, default 0)
  The manifold and k come from the run's own fit, so the test is not independent of the search that produced them;
  a manifold from tides constrained by dLOD alone (no climate index fitted) gives the clean version.
"""
import glob, json, os, re, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

def opt(flag, n=1, conv=str, default=None):
    if flag in sys.argv:
        i = sys.argv.index(flag); v = [conv(x) for x in sys.argv[i + 1:i + 1 + n]]; del sys.argv[i:i + 1 + n]; return v[0] if n == 1 else v
    return default
LIMBS = "--limbs" in sys.argv; sys.argv = [x for x in sys.argv if x != "--limbs"]
TURN = opt("--turn", conv=float); KH = opt("--kh", conv=float); HARM = opt("--harm", conv=int); HMAX = opt("--hmax", conv=int, default=9)
NNULL = opt("--null", conv=int, default=200); SEED = opt("--seed", conv=int, default=0)
OUT = opt("--out", default="lte_modulation.png"); OFFSET = opt("--offset", conv=float, default=0.0); FIT = opt("--fit", conv=int, default=0)
d, names = sys.argv[1], sys.argv[2:]
def leg_index(k):   # the leg's CLIMATE_INDEX file name without .dat, from lt.exe.resp.<k>
    try: return re.search(r"^\s*CLIMATE_INDEX\s+(\S+)", open(f"{d}/lt.exe.resp.{k}").read(), re.M).group(1).removesuffix(".dat")
    except (OSError, AttributeError): return None
avail = {}; k = 1
while os.path.exists(f"{d}/lte_results.{k}.csv"): avail[k] = leg_index(k); k += 1
if not names: names = [avail[k] or f"leg {k}" for k in avail]   # no names: every leg
legs = []   # each name picks the leg whose index file it names (case-insensitive), or a leg number
for n in names:
    hit = [k for k, ix in avail.items() if ix and ix.lower() == n.lower()] or ([int(n)] if n.isdigit() and int(n) in avail else [])
    if not hit: sys.exit(f"{n}: no leg with that index; legs here: " + ", ".join(f"{k} {ix}" for k, ix in avail.items()))
    legs.append(hit[0])
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]   # same order as leg_scores2.py
CW = (365.24219 / (27.212221 / 2.0)) % 1.0   # CW alias, cycles per year

pre = any(re.search(r"^\s*(export\s+.*\bMANIFOLD_PRE=(TRUE|true|1)\b|MANIFOLD_PRE\s+(TRUE|true))", open(f).read(), re.M)
          for f in glob.glob(f"{d}/lte_run.sh") + glob.glob(f"{d}/lt.exe.resp"))

def decomp(v, ir):   # exact inverse of model(n) = orig(n) - IR*orig(n-12)
    u = v.copy()
    if ir != 0 and abs(ir) < 1:
        for i in range(12, len(v)): u[i] = v[i] + ir * u[i - 12]
    return u

def leg_p(i):
    p = json.load(open(f"{d}/lt.exe.p"))
    try: p.update(json.load(open(f"{d}/lt.exe.p.{i}")))
    except OSError: pass
    return p

def limbs(name, L, t, data, F, F4, k, ir, rng):
    """Cross-limb test (see --limbs): one harmonic of k fitted to the data on one limb, scored on the other."""
    tp = TURN if TURN is not None else float(t[np.argmax(F4)]); up, dn = t < tp, t >= tp
    if min(up.sum(), dn.sum()) < 120: print(f"{name:8s} leg {L}  turn {tp:.1f}: a limb is shorter than 10 years, skipped"); return
    dm = data.copy()
    if ir != 0: dm[12:] = data[12:] + ir * data[:-12]   # IR-matched data
    N = np.column_stack([np.ones_like(t), t - t[0], (t - t[0]) ** 2] + [g(2 * np.pi * f * t) for f in (1.0, 2.0, CW) for g in (np.sin, np.cos)])
    def score(kh):
        Z = np.column_stack([np.sin(2 * np.pi * kh * F), np.cos(2 * np.pi * kh * F)]); NZ = np.column_stack([N, Z]); out = []
        for A_, B_ in ((up, dn), (dn, up)):
            c = np.linalg.lstsq(NZ[A_], dm[A_], rcond=None)[0]
            r = dm[B_] - N[B_] @ np.linalg.lstsq(N[B_], dm[B_], rcond=None)[0]; out.append(float(np.corrcoef(Z[B_] @ c[-2:], r)[0, 1]))
        return out
    if KH is not None: kh, tag = KH, f"k_h {KH:+.4f} (fixed)"
    else:
        h = HARM or max(range(1, HMAX + 1), key=lambda j: sum(score(j * k))); kh = h * k; tag = f"h={h}{'' if HARM else '*'} (h*k {kh:+.4f})"
    s_ = score(kh)
    nul = np.array([sum(score(kr)) / 2 for kr in rng.uniform(0.05, 3.0, NNULL)]) if NNULL > 0 else None
    print(f"{name:8s} leg {L}  turn {tp:.0f}  {tag}  rise->fall {s_[0]:+.2f}  fall->rise {s_[1]:+.2f}"
          + (f"  | random k: mean {nul.mean():+.2f}, 95th pct {np.percentile(nul, 95):+.2f}, max {nul.max():+.2f}" if nul is not None else ""))

rng = np.random.default_rng(SEED)
fig, A = plt.subplots(figsize=(13, 7.5), facecolor=SURF); A.set_facecolor(SURF)
print(f"{d}: column 4 is {'pre-Bessel (MANIFOLD_PRE), converted to post-Bessel' if pre else 'post-Bessel'}")
if LIMBS:
    print("--limbs: one sinusoid of the manifold fitted to the data on one limb, CC on the other"
          + ("" if HARM or KH is not None else f"  (* harmonic chosen from 1..{HMAX} on both limbs)"))
for i, (name, L) in enumerate(zip(names, legs)):   # i: colour slot, L: leg number
    a = np.loadtxt(f"{d}/lte_results.{L}.csv", delimiter=","); t, model, data, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]; F4 = F.copy()
    have = np.abs(data) > 0 if (data[:60] == 0).all() else np.ones_like(t, dtype=bool)   # --limbs: skip leading rows without data
    w = json.load(open(f"{d}/lt.exe.windings.{L}.json")); m = w["manifold"]; nm = len(leg_p(L)["ltep"]); ir = float(w["IR"])
    if pre:   # the solver's Bessel step (gem-lte-primitives-solution.adb), this leg's own parameters
        eS, eC, bg = m["impA"], m["impB"], m["bg"]
        kk, e2s, e2c = (m["besk"] * (1.0 - bg), m["offs"], m.get("bg2", bg)) if nm == 1 else (m["besk"], m["offs"], bg)
        ph = 2 * np.pi * kk * F; F = F + eS * np.sin(ph) + eC * np.cos(ph) + e2s * eS * np.sin(2 * ph) + e2c * eC * np.cos(2 * ph)
    k = float(w["k_amp_phase"][nm - 1][0])   # base winding = last ltep
    u = decomp(model, ir); ref = float(w.get("secular_context", {}).get("accel_ref", t[0]))
    X = np.column_stack([np.ones_like(t), F, t - t[0], (t - ref) ** 2] + [g(2 * np.pi * f * t) for f in (1.0, 2.0, CW) for g in (np.sin, np.cos)])
    y = u - X @ np.linalg.lstsq(X, u, rcond=None)[0]; y = y / y.std()
    x = np.mod(k * F, 1.0); o = np.argsort(x, kind="stable"); c = SERIES[i % len(SERIES)]; off = -OFFSET * i
    A.plot(x[o], y[o] + off, color=c, lw=0.7, alpha=0.55 if FIT else 0.8, label=f"{name} (leg {L}; k = {k:+.4f}, IR {ir:+.3f})")
    if LIMBS: limbs(name, L, t[have], data[have], F[have], F4[have], k, ir, rng)
    if FIT:
        B = np.column_stack([np.ones_like(x)] + [g(2 * np.pi * h * x) for h in range(1, FIT + 1) for g in (np.sin, np.cos)]); cf = np.linalg.lstsq(B, y, rcond=None)[0]
        xs = np.linspace(0, 1, 801); Bs = np.column_stack([np.ones_like(xs)] + [g(2 * np.pi * h * xs) for h in range(1, FIT + 1) for g in (np.sin, np.cos)])
        A.plot(xs, Bs @ cf + off, color=c, lw=2.2, label="_nolegend_"); r2 = 1 - np.var(y - B @ cf) / np.var(y)
        print(f"{name:8s} leg {L}  k {k:+.5f}  IR {ir:+.3f}  n {len(x)}  share of the modulation on one zone curve (h<={FIT}): {r2:.2f}")
    else:
        print(f"{name:8s} leg {L}  k {k:+.5f}  IR {ir:+.3f}  n {len(x)}")
if FIT: A.plot([], [], color=MUTED, lw=2.2, label=f"bold: zone-periodic fit, harmonics 1..{FIT}")
A.set_xlim(0, 1); A.set_xlabel("manifold folded into one Brillouin zone, frac(k F)" + ("  (post-Bessel F)" if pre else ""), color=MUTED)
A.set_ylabel("LTE modulation: model, IR removed, non-zone terms removed (standardized)" + (f", legs offset by {OFFSET:g}" if OFFSET else ""), color=MUTED)
A.set_title(f"{os.path.basename(os.path.abspath(d))}: LTE modulation of each leg within one zone, points sorted by folded manifold", loc="left", fontsize=12, color=INK)
A.grid(True, color=GRID, lw=0.6); A.set_axisbelow(True); A.tick_params(colors=MUTED)
for s in ("top", "right"): A.spines[s].set_visible(False)
for s in ("left", "bottom"): A.spines[s].set_color(GRID)
A.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.0, 1.0)); fig.tight_layout()
out = OUT if os.path.isabs(OUT) else f"{d}/{OUT}"; fig.savefig(out, dpi=100, facecolor=SURF); print("plot:", out)
