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

usage: plot_lpap_year.py [DIR] [--year 2000] [--years 3] [--dlod FILE | --no-dlod] [--measured] [--no-seasonal] [--csv FILE] [--out FILE] [--no-show]
  DIR           run folder with lt.exe.p (default: current folder)
  --year Y      start of the span, decimal year (default 2000); --years N its length in years (default 3)
  --dlod FILE   dLOD-calibrated lpap (default: experiments/wisdom_lpap/dlod_calibrated.lt.exe.p)
  --measured    also draw the measured dLOD (Feb2026/dlod3.dat, 1962-2019) for that year
  --no-seasonal omit annual (~365.24 d) and semi-annual (~182.62 d, Ssa) lines from the composed sums
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
ap.add_argument("--no-seasonal", "--no-annual", dest="no_seasonal", action="store_true",
                help="omit annual (~365.24 d) and semi-annual (~182.62 d, Ssa) lines from the composed sums")
ap.add_argument("--no-long", action="store_true",
                help="omit multi-year (|P| > 365.25 d) nodal/perigee terms from the composed sums")
ap.add_argument("--poly-tier", default=None,
                help="override POLY_TIER (auto-detected from lt.exe.resp when TIDES=POLY)")
ap.add_argument("--csv", default=None, help="write the daily series to this CSV")
ap.add_argument("--out", default="lpap_year.png", help="output image (default %(default)s; 'none' to skip)")
ap.add_argument("--no-show", action="store_true")
a = ap.parse_args()
if a.no_show or not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

def lpap_of(path):
    p = json.load(open(path)); return p, np.array(p["lpap"] if isinstance(p, dict) else p, float)

def resp_val(d, key, default=None):   # later lines win, '-' comments a line out
    v = None
    f = os.path.join(d, "lt.exe.resp")
    if os.path.exists(f):
        pat = re.compile(rf"^\s*{re.escape(key)}\s+(\S+)", re.I)
        for line in open(f):
            m = pat.match(line)
            if m: v = m.group(1)
    return v if v is not None else os.environ.get(key, default)

def resp_year(d):
    return float(resp_val(d, "YEAR", 0.0))

def is_seasonal(P):
    ap_ = abs(P)
    return abs(ap_ - 365.2422) < 20.0 or abs(ap_ - 182.6211) < 12.0

# Exact Meeus-free GEM.Tidal_Poly catalog matching src/gem-tidal_poly.adb
_POLY_CATALOG = [
    ( 2,  0,  0,  0, 0),  #  1: Mf(2L)      13.6608 d
    ( 2,  0, -1,  0, 0),  #  2: Mf_nod      13.6334 d
    ( 2,  0, -2,  0, 0),  #  3: 2F          13.6061 d
    ( 0,  1,  0,  0, 0),  #  4: Mm(Mp)      27.5545 d
    ( 0, -1,  0,  2, 0),  #  5: MSm(ev)     31.8119 d
    ( 2,  1,  0,  0, 0),  #  6: Mt           9.1329 d
    ( 2,  1, -1,  0, 0),  #  7: Mt_nod       9.1207 d
    ( 0,  0,  0,  2, 0),  #  8: MSf(2D)     14.7653 d
    ( 2,  0,  0, -2, 0),  #  9: Ssa(2Ls)   182.6211 d
    ( 0,  1,  0,  2, 0),  # 10: 2D+Mp        9.6137 d
    ( 1,  0,  0,  0, 0),  # 11: T(L)        27.3216 d
    ( 1,  0, -1,  0, 0),  # 12: D(F)        27.2122 d
    ( 0,  1, -1,  0, 0),  # 13: Mp-Om       27.4432 d
    ( 2, -1, -1,  0, 0),  # 14: 2L-Mp-Om    26.9850 d
    ( 2, -1,  0,  2, 0),  # 15: Mtm          9.5569 d
    ( 1,  0,  0, -1, 1),  # 16: K1           0.9973 d
    (-1,  0,  0, -1, 1),  # 17: O1           1.0758 d
    (-1,  0,  0,  1, 1),  # 18: P1           1.0027 d
    ( 0,  0,  0, -2, 2),  # 19: M2           0.5175 d
    ( 0, -1,  0, -2, 2),  # 20: N2           0.5274 d
    ( 2,  0,  0, -2, 2),  # 21: K2           0.4986 d
    ( 0,  0, -1,  0, 0),  # 22: N(18.6y)  6798.38 d
    ( 0,  0, -2,  0, 0),  # 23: 2N(9.3y)  3399.19 d
    ( 1, -1,  0,  0, 0),  # 24: p(8.85y)  3231.50 d
    ( 2, -2,  0,  0, 0),  # 25: 2p(4.43y) 1615.75 d
    ( 1, -1, -1,  0, 0),  # 26: p+N(6.0y) 2190.24 d
    ( 0, -1, -1,  2, 0),  # 27: MSm_nod     31.6639 d
    ( 0,  0, -1,  2, 0),  # 28: MSf_nod+    14.7333 d
    ( 0,  0,  1,  2, 0),  # 29: MSf_nod-    14.7975 d
    ( 2, -1, -1,  2, 0),  # 30: Mtm_nod      9.5435 d
    ( 2,  2,  0,  0, 0),  # 31: Mqm          6.8594 d
    ( 2,  0,  0,  2, 0),  # 32: Msqm         7.0958 d
    ( 0,  2,  0,  0, 0),  # 33: 2Mp         13.7773 d
    ( 2,  2, -1,  0, 0),  # 34: Mqm_nod      6.8525 d
    ( 2,  0, -1,  2, 0),  # 35: Msqm_nod     7.0884 d
    ( 2, -2, -2,  0, 0),  # 36: 2p+2N(3y)  1095.18 d
    ( 2, -2, -1,  0, 0),  # 37: 2p+N(3.6y) 1305.66 d
    ( 1, -1,  1,  0, 0),  # 38: p-N(16.9y) 6167.21 d
    (-2,  2, -1,  0, 0),  # 39: -2p+N(5.8) -2120.9 d
    ( 1, -1, -2,  0, 0),  # 40: p+2N(4.5y) 1656.32 d
    ( 0,  1,  1,  0, 0),  # 41: Mp+Om       27.6668 d
    ( 2,  1, -2,  0, 0),  # 42: Mt_2n        9.1085 d
    ( 0,  3,  0,  0, 0),  # 43: 3Mp          9.1848 d
    ( 2,  1,  0,  2, 0),  # 44: M5_1         5.6427 d
    ( 2, -1,  0,  0, 0),  # 45: 2L-Mp       27.0927 d
    ( 2,  1, -1,  2, 0),  # 46: M5_1n        5.6380 d
    ( 2,  3,  0,  0, 0),  # 47: 2L+3Mp       5.4922 d
    (-1,  0,  0,  3, 0),  # 48: 3D-L        15.3874 d
    ( 2,  1,  0, -2, 0),  # 49: 2L+Mp-2D    23.9421 d
    ( 3,  0,  0, -1, 1),  # 50: OO1          0.9294 d
    ( 0,  0,  0, -1, 1),  # 51: M1           1.0351 d
    ( 1,  1,  0, -1, 1),  # 52: J1           0.9624 d
    (-1, -1,  0, -1, 1),  # 53: Q1           1.1195 d
]
_POLY_TIERS = {
    "Z9": list(range(0, 9)), "Z18": list(range(0, 9)),
    "Z10": list(range(0, 10)), "Z20": list(range(0, 10)),
    "Z12": list(range(0, 10)) + [30, 31],
    "Z15": list(range(0, 15)), "Z30": list(range(0, 15)),
    "Z10_D6": list(range(0, 10)) + list(range(15, 21)), "F32": list(range(0, 10)) + list(range(15, 21)),
    "Z15_D6": list(range(0, 21)),
    "Z10_N5": list(range(0, 10)) + list(range(21, 26)),
    "Z10_D6_N5": list(range(0, 10)) + list(range(15, 26)),
    "Z10_S5_N5": list(range(0, 10)) + [14, 26, 27, 28, 29] + list(range(21, 26)),
    "Z12_S5_N5": list(range(0, 10)) + [30, 31, 14, 26, 27, 28, 29] + list(range(21, 26)),
    "Z10_S5_D6_N5": list(range(0, 10)) + [14, 26, 27, 28, 29] + list(range(15, 26)),
    "Z12_S5_D6_N5": list(range(0, 10)) + [30, 31, 14, 26, 27, 28, 29] + list(range(15, 26)),
    "Z35": list(range(0, 35)),
    "Z40": list(range(0, 40)),
    "Z52": list(range(0, 26)) + list(range(27, 53)),
    "FULL": list(range(0, 26)) + list(range(27, 53)),
}
_R_L  = 481267.88123421 / (360.0 * 36525.0)
_R_Mp = 477198.8675055  / (360.0 * 36525.0)
_R_F  = 483202.0175233  / (360.0 * 36525.0)
_R_D  = 445267.1114034  / (360.0 * 36525.0)
_R_Om = _R_L - _R_F

def compose_poly(lp, t, yl, tier, drop_seasonal=False, drop_long=False):
    idxs = _POLY_TIERS[tier.upper()]
    jd = 2448074.761995 + (t - 1990.5) * yl
    T = (jd - 2451545.0) / 36525.0
    d2r = np.pi / 180.0
    L  = (218.3164477 + 481267.88123421 * T) * d2r
    Mp = (134.9633964 + 477198.8675055  * T) * d2r
    F  = ( 93.2720950 + 483202.0175233  * T) * d2r
    D  = (297.8501921 + 445267.1114034  * T) * d2r
    Om = L - F
    th_d = 2.0 * np.pi * ((jd - 2451545.0) - np.floor(jd - 2451545.0))
    y = np.zeros_like(t)
    kept = 0
    for k, cat_i in enumerate(idxs):
        if k >= len(lp): break
        A, ph = float(lp[k, 1]), float(lp[k, 2])
        nL, nMp, nOm, nD, nDay = _POLY_CATALOG[cat_i]
        rate = nL * _R_L + nMp * _R_Mp + nOm * _R_Om + nD * _R_D + nDay
        per = 1.0 / rate if abs(rate) > 1e-15 else 0.0
        if drop_seasonal and is_seasonal(per): continue
        if drop_long and abs(per) > 365.25: continue
        if A != 0.0: kept += 1
        theta = nL * L + nMp * Mp + nOm * Om + nD * D + nDay * th_d
        y += A * np.sin(theta + ph)
    return y, kept

def compose(lp, t, yl, drop_seasonal=False, drop_long=False):
    y = np.zeros_like(t)
    kept = 0
    for P, A, ph in lp[:, :3]:
        if P == 0:
            continue
        if drop_seasonal and is_seasonal(P):
            continue
        if drop_long and abs(P) > 365.25:
            continue
        if A != 0:
            kept += 1
        y += A * np.sin(2 * np.pi * (yl / P) * t + ph)
    return y, kept

fit_path = os.path.join(a.dir, "lt.exe.p")
if not os.path.exists(fit_path): sys.exit(f"no lt.exe.p in {a.dir}")
P, LP = lpap_of(fit_path)
YEAR = resp_year(a.dir); yl = 365.2422484 + YEAR + float(P.get("year", 0.0) if isinstance(P, dict) else 0.0)
tides_mode = str(resp_val(a.dir, "TIDES", "TABLE")).upper()
poly_tier = a.poly_tier or (str(resp_val(a.dir, "POLY_TIER", "Z10")).upper() if tides_mode == "POLY" else None)
n = int(round(a.years * 365.25)) + 1
t = a.year + np.arange(n) / 365.25                     # daily, as decimal years
days = t                                               # x-axis: decimal year
_sd = 1 if a.years >= 2 else 2 if a.years >= 0.2 else 3
span = f"{a.year:.{_sd}f}-{a.year + a.years:.{_sd}f}"
if poly_tier:
    y_fit, n_fit = compose_poly(LP, t, yl, poly_tier, drop_seasonal=a.no_seasonal, drop_long=a.no_long)
else:
    y_fit, n_fit = compose(LP, t, yl, drop_seasonal=a.no_seasonal, drop_long=a.no_long)
series = {"fit": y_fit}
name = os.path.basename(os.path.abspath(a.dir))
flags_tag = (", no ann/semi" if a.no_seasonal else "") + (", no >1y" if a.no_long else "")
mode_tag = f"POLY {poly_tier}, " if poly_tier else ""

fig, ax = plt.subplots(figsize=(14, 5.5), facecolor="#fcfcfb"); ax.set_facecolor("#fcfcfb")
ax.plot(days, y_fit, color="red", lw=1.3, zorder=2, label=f"{name} lt.exe.p ({mode_tag}{n_fit} lines{flags_tag})")
title = f"Composed tidal sum, daily, {span} — YL {yl:.5f} d" + (f" [{mode_tag.rstrip(', ')}]" if poly_tier else "") + (f" [{flags_tag.lstrip(', ')}]" if flags_tag else "")
if not a.no_dlod:
    if not os.path.exists(a.dlod): sys.exit(f"dLOD-calibrated lpap not found: {a.dlod} (use --no-dlod)")
    _, LD = lpap_of(a.dlod); y_d, _ = compose(LD, t, yl, drop_seasonal=a.no_seasonal, drop_long=a.no_long); series["dlod_calibrated"] = y_d
    cc = float(np.corrcoef(y_fit, y_d)[0, 1])
    ax.plot(days, y_d, color="blue", lw=1.2, ls=":", alpha=0.9, zorder=3, label="dLOD calibrated (" + os.path.basename(a.dlod) + ")")
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
