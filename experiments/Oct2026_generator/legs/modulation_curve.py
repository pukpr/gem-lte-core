#!/usr/bin/env python3
"""Complete LTE modulation of each leg as a function of the pre-Bessel forcing F, from lt.exe.windings.<leg>.json.

For F on a grid (default 0 .. 35, the range of the pre-Bessel manifold), each leg's curve is
    G(F) = sum over all NM + NH windings of  A * sin(2 pi k F_post + phase),    F_post = Bessel(F)
with the leg's own winding wavenumbers, amplitudes and phases (k_amp_phase) and its own Bessel step from the JSON's
"manifold" block, exactly as the solver applies it:
    NM >= 2:  F_post = F + impA sin(2 pi besk F) + impB cos(2 pi besk F) + offs impA sin(4 pi besk F) + bg impB cos(4 pi besk F)
    NM = 1:   the same with besk (1 - bg) as the wavenumber and bg2 in place of bg
(LAYER and FRICTION runs use other steps and are not covered.) The JSON amplitudes and phases are the solver's own
regression coefficients. What is left out: level, trend, acceleration, the annual/semi-annual/CW terms (functions of
time, not of F), and IR and the 9-point filter (applied along time). --linear adds the k0 * F_post term.
Top panel: G(F), one line per leg. Bottom panel: each leg's Bessel correction, F_post - F.

--check: one panel per leg comparing the curve with the solver's own output: lte_results column 2 (the model) with
IR undone (exact recursive inverse), less k0 F_post and the time-dependent terms (level, trend, acceleration,
annual, semi-annual, CW and the annual impulse; least squares on the model minus the curve, so they cannot absorb
the modulation), against column 4 (pre-Bessel F, MANIFOLD_PRE), sorted. The 9-point time filter on the model
cannot be undone; the same filter applied to the curve along time shows how much of any mismatch it explains.

--std sets the scale of each curve:
  data   (default with more than one leg) divided by the leg's index standard deviation (lte_results column 3), so the
         curves are in units of each index's sigma and a leg's amplitude shows how much of its index it carries
  curve  divided by the curve's own rms on the plotted range: equal sizes, shapes and phases only (--norm is the same)
  raw    (default for a single leg) each leg in its index's own units

usage: modulation_curve.py [DIR] [NAME ...] [--fmin 0] [--fmax 35] [--std data|curve|raw] [--offset Y] [--linear] [--check] [--out FILE]
Each NAME picks a leg by its lt.exe.resp.<leg> CLIMATE_INDEX file (case-insensitive, without .dat) or by number;
no names: every leg. A single-leg run uses lt.exe.windings.json.
"""
import argparse, glob, json, os, re, sys
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dir", nargs="?", default="."); ap.add_argument("names", nargs="*")
ap.add_argument("--fmin", type=float, default=0.0); ap.add_argument("--fmax", type=float, default=35.0)
ap.add_argument("--n", type=int, default=14001, help="grid points (default %(default)s)")
ap.add_argument("--std", choices=["data", "curve", "raw"], default=None,
                help="scale: data = / index std (default with several legs), curve = / curve rms, raw = own units (default, one leg)")
ap.add_argument("--norm", action="store_true", help="same as --std curve")
ap.add_argument("--offset", type=float, default=0.0, help="shift each successive leg down by this much")
ap.add_argument("--linear", action="store_true", help="include the k0 * F_post term")
ap.add_argument("--check", action="store_true", help="per-leg panels: the curve against lte_results column 2 vs column 4 (see above)")
ap.add_argument("--out", default=None, help="default modulation_curve.png (modulation_check.png with --check)")
a = ap.parse_args()
d = a.dir
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

def leg_index(k):
    try: return re.search(r"^\s*CLIMATE_INDEX\s+(\S+)", open(f"{d}/lt.exe.resp.{k}").read(), re.M).group(1).removesuffix(".dat")
    except (OSError, AttributeError): return None
avail = {}; k = 1
while os.path.exists(f"{d}/lt.exe.windings.{k}.json"): avail[k] = leg_index(k); k += 1
single = not avail
if single:
    if not os.path.exists(f"{d}/lt.exe.windings.json"): sys.exit(f"no lt.exe.windings*.json in {d}")
    legs, labels = [None], [a.names[0] if a.names else os.path.basename(os.path.abspath(d))]
else:
    labels = a.names or [avail[k] or f"leg {k}" for k in avail]; legs = []
    for n in labels:
        hit = [k for k, ix in avail.items() if ix and ix.lower() == n.lower()] or ([int(n)] if n.isdigit() and int(n) in avail else [])
        if not hit: sys.exit(f"{n}: no leg with that index; legs here: " + ", ".join(f"{k} {ix}" for k, ix in avail.items()))
        legs.append(hit[0])

def nm_of(L):
    p = json.load(open(f"{d}/lt.exe.p"))
    if L is not None and os.path.exists(f"{d}/lt.exe.p.{L}"): p.update(json.load(open(f"{d}/lt.exe.p.{L}")))
    return len(p.get("ltep", [0, 0, 0]))
def bessel(F, m, nm):
    bg = m["bg"]
    kk, e2s, e2c = (m["besk"] * (1 - bg), m["offs"], m.get("bg2", bg)) if nm == 1 else (m["besk"], m["offs"], bg)
    ph = 2 * np.pi * kk * F
    return F + m["impA"] * np.sin(ph) + m["impB"] * np.cos(ph) + e2s * m["impA"] * np.sin(2 * ph) + e2c * m["impB"] * np.cos(2 * ph)

if a.out is None: a.out = "modulation_check.png" if a.check else "modulation_curve.png"
F = np.linspace(a.fmin, a.fmax, a.n)
mode = "curve" if a.norm else (a.std or ("data" if len(legs) > 1 else "raw"))
if a.check:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Feb2026"))
    from lte_forward import filter9point
    pre = any(re.search(r"^\s*(export\s+.*\bMANIFOLD_PRE=(TRUE|true|1)\b|MANIFOLD_PRE\s+(TRUE|true))", open(f).read(), re.M)
              for f in glob.glob(f"{d}/lte_run.sh") + glob.glob(f"{d}/lt.exe.resp"))
    if not pre: sys.exit("--check needs column 4 to be the pre-Bessel manifold (MANIFOLD_PRE TRUE in lte_run.sh or lt.exe.resp)")
    def env_of(key, default):
        v = re.findall(r"^\s*export\s+[^#\n]*?\b" + key + r"=([^\s#]+)", open(f"{d}/lte_run.sh").read(), re.M) if os.path.exists(f"{d}/lte_run.sh") else []
        return v[-1] if v else default
    F9 = int(float(env_of("F9", "1")))
    n_leg = len(legs); fig, axs = plt.subplots(n_leg, 1, figsize=(15, 2.6 * n_leg + 1), facecolor=SURF, sharex=True, squeeze=False)
    for i, (L, name) in enumerate(zip(legs, labels)):
        sfx = "" if L is None else f".{L}"
        w = json.load(open(f"{d}/lt.exe.windings{sfx}.json")); m, nm = w["manifold"], nm_of(L)
        r = np.loadtxt(f"{d}/lte_results{sfx}.csv", delimiter=","); t, model, Fd = r[:, 0], r[:, 1], r[:, 3]
        ir = float(w["IR"]); u = model.copy()
        if ir != 0 and abs(ir) < 1:
            for j in range(12, len(u)): u[j] = model[j] + ir * u[j - 12]
        Fp = bessel(Fd, m, nm); Gd = sum(Am * np.sin(2 * np.pi * kw * Fp + P) for kw, Am, P in w["k_amp_phase"])
        p_ = json.load(open(f"{d}/lt.exe.p"))
        if L is not None and os.path.exists(f"{d}/lt.exe.p.{L}"): p_.update(json.load(open(f"{d}/lt.exe.p.{L}")))
        ref = float(w.get("secular_context", {}).get("accel_ref", t[0])); cwf = (365.2422484 / (27.212220815 / 2.0)) % 1.0
        trunc = np.rint((t - np.floor(t)) * 12).astype(int); imp = (trunc == int(round(abs(float(p_.get("delB", 0))) * 12)) % 12).astype(float)
        X = np.column_stack([np.ones_like(t), t - ref, (t - ref) ** 2, imp] + [g(2 * np.pi * f * t) for f in (1.0, 2.0, cwf) for g in (np.sin, np.cos)])
        base = u - w["k0"] * Fp - Gd; y = u - w["k0"] * Fp - X @ np.linalg.lstsq(X, base, rcond=None)[0]   # data-derived modulation
        Gs = Gd.copy()
        for _ in range(F9): Gs = filter9point(Gs)                    # the curve through the same time filter
        cc_raw, cc_f9 = np.corrcoef(y, Gd)[0, 1], np.corrcoef(y, Gs)[0, 1]
        e_raw, e_f9 = np.std(y - Gd) / np.std(Gd), np.std(y - Gs) / np.std(Gd)
        o = np.argsort(Fd, kind="stable"); A = axs[i, 0]; A.set_facecolor(SURF); c = SERIES[i % len(SERIES)]
        A.plot(F, sum(Am * np.sin(2 * np.pi * kw * bessel(F, m, nm) + P) for kw, Am, P in w["k_amp_phase"]), color=c, lw=0.9, alpha=0.55, label="curve from lt.exe.windings json")
        A.plot(Fd[o], y[o], color="#000000", lw=0.8, ls=":", label="lte_results: column 2 (IR and time terms removed) vs column 4, sorted")
        A.set_title(f"{name}" + (f" (leg {L})" if L else "") + f": data vs curve CC {cc_raw:+.4f}, rms misfit {e_raw:.2f} of the curve's rms; "
                    f"with the curve also 9-point filtered CC {cc_f9:+.4f}, misfit {e_f9:.3f}", loc="left", fontsize=9.5, color=INK)
        print(f"{name:8s} {'leg ' + str(L) if L else '':6s} data F {Fd.min():.2f}..{Fd.max():.2f} | data vs curve CC {cc_raw:+.4f} (misfit {e_raw:.2f}); "
              f"vs the filtered curve CC {cc_f9:+.4f} (misfit {e_f9:.3f})")
        A.grid(True, color=GRID, lw=0.6); A.set_axisbelow(True); A.tick_params(colors=MUTED); A.set_xlim(a.fmin, a.fmax)
        for s_ in ("top", "right"): A.spines[s_].set_visible(False)
        if i == 0: A.legend(frameon=False, fontsize=8.5, loc="upper right", ncol=2)
    axs[-1, 0].set_xlabel("pre-Bessel forcing F (lte_results column 4)", color=MUTED)
    fig.suptitle(f"{os.path.basename(os.path.abspath(d))}: modulation curves checked against the solver's model", x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = a.out if os.path.isabs(a.out) else os.path.join(d, a.out); fig.savefig(out, dpi=110, facecolor=SURF); print("plot:", out); sys.exit()
F = np.linspace(a.fmin, a.fmax, a.n)
fig, (A, Bx) = plt.subplots(2, 1, figsize=(15, 9), facecolor=SURF, sharex=True, gridspec_kw={"height_ratios": [3, 1.2]})
for ax in (A, Bx): ax.set_facecolor(SURF)
for i, (L, name) in enumerate(zip(legs, labels)):
    w = json.load(open(f"{d}/lt.exe.windings.json" if L is None else f"{d}/lt.exe.windings.{L}.json"))
    m, nm = w["manifold"], nm_of(L); Fp = bessel(F, m, nm)
    G = sum(Am * np.sin(2 * np.pi * kw * Fp + P) for kw, Am, P in w["k_amp_phase"])
    if a.linear: G = G + w["k0"] * Fp
    if mode == "curve": scale = G.std()
    elif mode == "data":
        dd = np.loadtxt(f"{d}/lte_results.csv" if L is None else f"{d}/lte_results.{L}.csv", delimiter=",")[:, 2]
        dd = dd[dd != 0] if (dd[:60] == 0).all() else dd; scale = dd.std()
    else: scale = 1.0
    c = SERIES[i % len(SERIES)]
    nk = len(w["k_amp_phase"]); dom = max(w["k_amp_phase"], key=lambda r: r[1])
    A.plot(F, G / scale - a.offset * i, color=c, lw=1.1, label=f"{name}" + (f" (leg {L})" if L else "") + f": {nk} windings, largest k {dom[0]:+.3f} (amp {dom[1]:.3f}); besk {m['besk']:.4f}")
    Bx.plot(F, Fp - F, color=c, lw=1.0)
    print(f"{name:8s} {'leg ' + str(L) if L else '':6s} NM {nm}, {nk} windings, rms {G.std():.4f}, scale ({mode}) {scale:.4f}; Bessel: impA {m['impA']:+.3f} impB {m['impB']:+.3f} besk {m['besk']:.4f} offs {m['offs']:+.4f} bg {m['bg']:+.3f}; F_post - F in [{(Fp - F).min():+.2f}, {(Fp - F).max():+.2f}]")
A.axhline(0, color=MUTED, lw=0.5); Bx.axhline(0, color=MUTED, lw=0.5)
A.set_ylabel("modulation G(F)" + {"curve": " / its own rms", "data": " / index std (units of each index's sigma)", "raw": " (each leg in its own units)"}[mode] + (f", legs offset by {a.offset:g}" if a.offset else ""), color=MUTED)
Bx.set_ylabel("Bessel correction\nF_post - F", color=MUTED); Bx.set_xlabel("pre-Bessel forcing F (lte_results column 4 with MANIFOLD_PRE)", color=MUTED)
A.set_title(f"{os.path.basename(os.path.abspath(d))}\ncomplete modulation of each leg: all NM+NH windings on the leg's Bessel-corrected F"
            + (", with k0 F" if a.linear else ""), loc="left", fontsize=11.5, color=INK)
for ax in (A, Bx):
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(colors=MUTED); ax.set_xlim(a.fmin, a.fmax)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
h, lab = A.get_legend_handles_labels()
fig.legend(h, lab, frameon=False, fontsize=8.5, loc="lower center", ncol=2)
fig.tight_layout(rect=(0, 0.04 + 0.022 * ((len(lab) + 1) // 2), 1, 1))
out = a.out if os.path.isabs(a.out) else os.path.join(d, a.out); fig.savefig(out, dpi=110, facecolor=SURF); print("plot:", out)
