#!/usr/bin/env python3
"""Composed sum of a fit's tidal factors (lpap in lt.exe.p), sampled daily over a span of years (default 3), with the
dLOD-calibrated tidal factors composed the same way on the same chart.

The sum is the solver's own tidal series (Tide_Sum):
    y(t) = sum_i A_i * sin(2 pi (YL / P_i) t + phi_i),   t in decimal years,   YL = 365.2422484 + YEAR + year
with YEAR from lt.exe.resp (or the environment) and year from lt.exe.p, so the curves sit at their true calendar
dates. (lte_gui.py's composed-sum plot counts time in days from 1 January of an arbitrary start year instead, so its
phases are not the solver's.) Composed this way, the dLOD-calibrated lines reproduce the measured dLOD record at
CC 0.97, as the solver reports. Daily sampling matches the dLOD record's own; the diurnal lines alias, as there.
The impulse comb, integrator and shfT term are not applied: this is the tidal sum that feeds them.

usage: plot_lpap_year.py [DIR] [--year 2000] [--years 3] [--dlod FILE | --no-dlod] [--measured] [--csv FILE] [--out FILE] [--no-show]
  DIR          run folder with lt.exe.p (default: current folder)
  --year Y     start of the span, decimal year (default 2000); --years N its length in years (default 3)
  --dlod FILE  dLOD-calibrated lpap (default: experiments/wisdom_lpap/dlod_calibrated.lt.exe.p)
  --measured   also draw the measured dLOD (Feb2026/dlod3.dat, 1962-2019) for that year
"""
import argparse, json, os, re, sys
import numpy as np
import matplotlib

HERE = os.path.dirname(os.path.abspath(__file__))
DLOD_DEFAULT = os.path.normpath(os.path.join(HERE, "..", "wisdom_lpap", "dlod_calibrated.lt.exe.p"))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dir", nargs="?", default=".")
ap.add_argument("--year", type=float, default=2000.0, help="start of the plotted span, decimal year (default %(default)s)")
ap.add_argument("--years", type=float, default=3.0, help="length of the plotted span in years (default %(default)s)")
ap.add_argument("--dlod", default=DLOD_DEFAULT, help="dLOD-calibrated lpap file (default: %(default)s)")
ap.add_argument("--no-dlod", action="store_true", help="omit the dLOD-calibrated curve")
ap.add_argument("--measured", action="store_true", help="also draw the measured dLOD record for that year")
ap.add_argument("--csv", default=None, help="write the daily series to this CSV")
ap.add_argument("--out", default="lpap_year.png", help="output image (default %(default)s; 'none' to skip)")
ap.add_argument("--no-show", action="store_true")
a = ap.parse_args()
if a.no_show or not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

def lpap_of(path):
    p = json.load(open(path)); return p, np.array(p["lpap"] if isinstance(p, dict) else p, float)

def resp_year(d):   # YEAR: lt.exe.resp (later lines win, '-' comments a line out), else the environment
    v = None
    f = os.path.join(d, "lt.exe.resp")
    if os.path.exists(f):
        for line in open(f):
            m = re.match(r"^\s*YEAR\s+(\S+)", line)
            if m: v = float(m.group(1))
    return v if v is not None else float(os.environ.get("YEAR", 0.0))

def compose(lp, t, yl):
    y = np.zeros_like(t)
    for P, A, ph in lp[:, :3]:
        if P != 0: y += A * np.sin(2 * np.pi * (yl / P) * t + ph)
    return y

fit_path = os.path.join(a.dir, "lt.exe.p")
if not os.path.exists(fit_path): sys.exit(f"no lt.exe.p in {a.dir}")
P, LP = lpap_of(fit_path)
YEAR = resp_year(a.dir); yl = 365.2422484 + YEAR + float(P.get("year", 0.0) if isinstance(P, dict) else 0.0)
n = int(round(a.years * 365.25)) + 1
t = a.year + np.arange(n) / 365.25                     # daily, as decimal years
days = t                                               # x-axis: decimal year
_sd = 1 if a.years >= 2 else 2 if a.years >= 0.2 else 3
span = f"{a.year:.{_sd}f}-{a.year + a.years:.{_sd}f}"
y_fit = compose(LP, t, yl)
series = {"fit": y_fit}
name = os.path.basename(os.path.abspath(a.dir))

fig, ax = plt.subplots(figsize=(14, 5.5), facecolor="#fcfcfb"); ax.set_facecolor("#fcfcfb")
ax.plot(days, y_fit, color="#2a78d6", lw=1.3, label=f"{name} lt.exe.p ({int((LP[:, 1] != 0).sum())} lines)")
title = f"Composed tidal sum, daily, {span} — YL {yl:.5f} d"
if not a.no_dlod:
    if not os.path.exists(a.dlod): sys.exit(f"dLOD-calibrated lpap not found: {a.dlod} (use --no-dlod)")
    _, LD = lpap_of(a.dlod); y_d = compose(LD, t, yl); series["dlod_calibrated"] = y_d
    cc = float(np.corrcoef(y_fit, y_d)[0, 1])
    ax.plot(days, y_d, color="#eb6834", lw=1.1, alpha=0.9, label="dLOD calibrated (" + os.path.basename(a.dlod) + ")")
    title += f"   |   CC fit vs dLOD {cc:+.3f}, rms ratio {np.std(y_fit) / np.std(y_d):.2f}"
    print(f"{span}: CC fit vs dLOD-calibrated {cc:+.4f}; rms fit {np.std(y_fit):.4f}, dLOD {np.std(y_d):.4f}")
if a.measured:
    m = np.loadtxt(os.path.join(HERE, "dlod3.dat")); sel = (m[:, 0] >= t[0]) & (m[:, 0] <= t[-1])
    if sel.any():
        dm = m[sel, 0]
        ax.plot(dm, m[sel, 1], color="#000000", lw=1.0, ls=":", label="measured dLOD (dlod3.dat)")
        series_m = np.interp(days, dm, m[sel, 1]); series["measured_dlod"] = series_m
        inside = (t >= dm[0]) & (t <= dm[-1])
        print(f"{span}: CC fit vs measured dLOD {np.corrcoef(y_fit[inside], series_m[inside])[0, 1]:+.4f}")
    else:
        print(f"measured dLOD covers 1962-2019; nothing in {span}")
ax.axhline(0, color="#52514e", lw=0.5)
from matplotlib.ticker import FormatStrFormatter
_dec = 1 if a.years >= 2 else 2 if a.years >= 0.2 else 3
ax.xaxis.set_major_formatter(FormatStrFormatter(f"%.{_dec}f"))   # plain decimal years, never an offset or exponent
ax.set_xlim(t[0], t[-1]); ax.set_xlabel("decimal year (daily samples)"); ax.set_ylabel("tidal sum (lpap units; dLOD calibration in ms)")
ax.set_title(title, loc="left", fontsize=11); ax.grid(True, ls=":", lw=0.5, alpha=0.7)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3); fig.tight_layout()
if a.csv:
    import csv
    with open(a.csv, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["decimal_year"] + list(series))
        for i in range(n): w.writerow([f"{t[i]:.6f}"] + [f"{series[k][i]:.6g}" for k in series])
    print("csv:", a.csv)
if a.out.lower() != "none":
    out = a.out if os.path.isabs(a.out) else os.path.join(a.dir, a.out); fig.savefig(out, dpi=120, facecolor="#fcfcfb"); print("plot:", out)
if not a.no_show and os.environ.get("DISPLAY"):
    plt.show()
