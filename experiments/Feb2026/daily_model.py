#!/usr/bin/env python3
"""Daily time series of a fitted LTE model, from a run directory's saved files.

Reads lte_results.csv (monthly: time, model, data, manifold), lt.exe.windings.json (the regression: level, k0,
trend, accel, IR and the winding amplitudes/phases), lt.exe.p (the annual and semi-annual regression terms) and,
for the filter setting, lte_run.sh.

The manifold (column 4) is a staircase: it changes in one month of the year, where the annual impulse lands
(two if a semi-annual impulse is active), and only drifts slightly otherwise. The monthly file cannot say on which
day inside that month the step happens. This script ASSUMES a date:

    --impulse-date MM-DD      calendar date of the (largest) annual step; default: the middle of the step interval

and builds the daily manifold as: constant before the date, the next monthly value after it, inside each step
interval; linear interpolation of the small drift in every other interval. A second impulse month, if present,
is placed the same number of days into its own interval.

The response is then evaluated day by day exactly as lt.exe does month by month:

    m(t) = level + k0 F(t) + trend t + accel (t - accel_ref)^2 + sum_j A_j sin(2 pi k_j F(t) + phi_j)
           + ann1 sin(2 pi t) + ann2 cos(2 pi t) + sem1 sin(4 pi t) + sem2 cos(4 pi t)      (from lt.exe.p)
    model(t) = m(t) - IR m(t - 1 year)                       (no lag term in the first year)

and, as a separate column, smoothed with the daily equivalent of the F9 filter the fit used (each pass of the
monthly 1-2-1 filter = a triangle two months wide).

Time convention: decimal year = year + (day of year - 1) / (days in that year); the monthly samples are taken at
their own decimal times (year + month/12).

Checks printed: the monthly reconstruction against column 2 (must agree to rounding), and the daily series sampled
at the monthly times.

usage: daily_model.py run_dir [--impulse-date MM-DD] [--start YEAR] [--end YEAR] [--out file.csv] [--plot [YEARS]]
"""
import argparse, datetime as dt, json, os, re, sys
import numpy as np

ap = argparse.ArgumentParser(description="daily LTE model series from a run directory")
ap.add_argument("run_dir")
ap.add_argument("--impulse-date", default=None, help="MM-DD of the annual manifold step (default: middle of the step interval)")
ap.add_argument("--start", type=float, default=None)
ap.add_argument("--end", type=float, default=None)
ap.add_argument("--out", default=None, help="output CSV (default: run_dir/daily_model.csv)")
ap.add_argument("--plot", nargs="?", const=12.0, type=float, default=None, help="also write run_dir/daily_model.png showing the last N years (default 12)")
ap.add_argument("--f9", type=int, default=None, help="filter passes for the smoothed column (default: F9 from lte_run.sh, else 0)")
args = ap.parse_args()
d = args.run_dir

a = np.loadtxt(os.path.join(d, "lte_results.csv"), delimiter=",")
t, mod, dat, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
w = json.load(open(os.path.join(d, "lt.exe.windings.json")))
sc = w.get("secular_context", {})
IR, level, k0, trend, accel = w["IR"], w["level"], w["k0"], w["trend"], w.get("accel", 0.0)
accel_ref = sc.get("accel_ref", t[0])
KAP = [(float(k), float(am), float(ph)) for k, am, ph in w["k_amp_phase"]]
pj = json.load(open(os.path.join(d, "lt.exe.p")))
ANN = [float(pj.get(k, 0.0)) for k in ("ann1", "ann2", "sem1", "sem2")]
for key, why in (("aero", "the aerosol term"), ("alpha", "the ALPHA gain"), ("layer", "the LAYER manifold"), ("annual_dither", "the annual dither")):
    v = w.get(key, sc.get(key, 0.0))
    if v not in (0, 0.0, False, None):
        print(f"  WARNING: this fit uses {why} ({key} = {v}), which this script does not reproduce", file=sys.stderr)
f9 = args.f9
if f9 is None:
    f9 = 0
    sh = os.path.join(d, "lte_run.sh")
    if os.path.exists(sh):
        m_ = re.search(r"^export F9=(-?\d+)", open(sh).read(), re.M)
        if m_:
            f9 = int(m_.group(1))

def response(tt, FF):
    m = level + k0 * FF + trend * tt + accel * (tt - accel_ref) ** 2
    for k, am, ph in KAP:
        m = m + am * np.sin(2 * np.pi * k * FF + ph)
    return m + ANN[0] * np.sin(2 * np.pi * tt) + ANN[1] * np.cos(2 * np.pi * tt) + ANN[2] * np.sin(4 * np.pi * tt) + ANN[3] * np.cos(4 * np.pi * tt)

# --- 1. the monthly model, reconstructed, as a check on the formula -------------------------------------------
mm = response(t, F)
mon = mm.copy()
mon[12:] = mm[12:] - IR * mm[:-12]
for _ in range(max(f9, 0)):
    s = mon.copy()
    s[1:-1] = 0.25 * mon[:-2] + 0.5 * mon[1:-1] + 0.25 * mon[2:]
    mon = s
cc = lambda x, y: float(np.corrcoef(x, y)[0, 1])
print(f"{os.path.abspath(d)}")
print(f"  monthly record {t[0]:.3f} to {t[-1]:.3f} ({len(t)} months); {len(KAP)} windings, IR {IR:+.4f}, F9 {f9}")
print(f"  check 1, monthly reconstruction against lte_results.csv column 2: max |difference| {np.abs(mon - mod)[13:-1].max():.2e}, CC {cc(mon[13:-1], mod[13:-1]):.6f}")

# --- 2. where the manifold steps ------------------------------------------------------------------------------
dF = np.diff(F)
moy = np.round((t[:-1] - np.floor(t[:-1] + 1e-9)) * 12).astype(int) % 12       # month-of-year index of each interval's start
size = np.array([np.abs(dF[moy == m]).mean() if (moy == m).any() else 0.0 for m in range(12)])
order = np.argsort(-size)
base = np.median(size)
steps = [int(m) for m in order[:2] if size[m] > 5 * max(base, 1e-12) and size[m] > 0.2 * size[order[0]]]
names = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
print("  manifold step interval(s): " + "; ".join(f"{names[m]} -> {names[(m + 1) % 12]} samples (mean |step| {size[m]:.2f})" for m in steps)
      + f"; mean |change| in the other months {np.mean([size[m] for m in range(12) if m not in steps]):.4f}")

def decyear(date):
    y0 = dt.date(date.year, 1, 1)
    return date.year + (date - y0).days / ((dt.date(date.year + 1, 1, 1) - y0).days)

main = steps[0]
if args.impulse_date:
    mo, da = (int(x) for x in args.impulse_date.split("-"))
    frac_year = decyear(dt.date(2001, mo, da)) - 2001                       # a non-leap reference year
    lo, hi = main / 12.0, (main + 1) / 12.0
    frac = (frac_year - lo) / (hi - lo)
    if not 0.0 <= frac <= 1.0:
        a0 = dt.date(2001, 1, 1) + dt.timedelta(days=round(lo * 365)); a1 = dt.date(2001, 1, 1) + dt.timedelta(days=round(hi * 365) - 1)
        sys.exit(f"  --impulse-date {args.impulse_date} is outside the interval in which this fit's manifold steps "
                 f"({a0:%m-%d} to {a1:%m-%d}); the monthly samples on either side fix it to that range")
else:
    frac = 0.5
day0 = dt.date(2001, 1, 1) + dt.timedelta(days=round((main + frac) / 12.0 * 365))
print(f"  assumed impulse date: {day0:%m-%d} (fraction {frac:.2f} of the step interval)"
      + (f"; second impulse {(dt.date(2001, 1, 1) + dt.timedelta(days=round(((steps[1] + frac) / 12.0) * 365))):%m-%d}" if len(steps) > 1 else ""))

# --- 3. daily grid and daily manifold -------------------------------------------------------------------------
y0 = int(np.floor(args.start if args.start is not None else t[0])); y1 = args.end if args.end is not None else t[-1]
first = dt.date(y0, 1, 1); last_year = int(np.floor(y1))
last = dt.date(last_year, 1, 1) + dt.timedelta(days=int(round((y1 - last_year) * 365.25)))
days = [first + dt.timedelta(days=i) for i in range((last - first).days + 1)]
td = np.array([decyear(x) for x in days])
keep = (td >= t[0]) & (td <= t[-1])
days = [x for x, k_ in zip(days, keep) if k_]; td = td[keep]
i = np.clip(np.searchsorted(t, td, side="right") - 1, 0, len(t) - 2)       # monthly interval of each day
u = (td - t[i]) / (t[i + 1] - t[i])
Fd = F[i] + u * (F[i + 1] - F[i])                                           # drift: linear
is_step = np.isin(moy[i], steps)
Fd[is_step] = np.where(u[is_step] < frac, F[i][is_step], F[i + 1][is_step])  # impulse: a step on the assumed date

# --- 4. response, lag term, filter ----------------------------------------------------------------------------
md = response(td, Fd)
lag = np.interp(td - 1.0, td, md)
raw = np.where(td - 1.0 >= td[0], md - IR * lag, md)
smooth = raw.copy()
if f9 > 0:
    n = int(round(365.2425 / 12))
    box = np.ones(n) / n
    pad = 2 * n * f9
    x = np.concatenate([np.full(pad, raw[0]), raw, np.full(pad, raw[-1])])
    for _ in range(2 * f9):                                                 # each 1-2-1 pass = two one-month boxcars
        x = np.convolve(x, box, mode="same")
    smooth = x[pad:-pad]
at = lambda v: np.interp(t, td, v)
ok = slice(14, -2)
print(f"  check 2, daily series sampled at the monthly times against column 2: unfiltered CC {cc(at(raw)[ok], mod[ok]):.4f}; filtered CC {cc(at(smooth)[ok], mod[ok]):.4f}, "
      f"max |difference| {np.abs(at(smooth) - mod)[ok].max():.3f} (model sd {mod.std():.3f})")
print(f"  daily series: {len(td)} days, {days[0]} to {days[-1]}; sd unfiltered {raw.std():.3f}, filtered {smooth.std():.3f}; "
      f"largest one-day change in the unfiltered series {np.abs(np.diff(raw)).max():.3f}")

out = args.out or os.path.join(d, "daily_model.csv")
with open(out, "w") as fh:
    fh.write("date,decimal_year,manifold,model_unfiltered,model_filtered\n")
    for x, tt, ff, r_, s_ in zip(days, td, Fd, raw, smooth):
        fh.write(f"{x.isoformat()},{tt:.6f},{ff:.6f},{r_:.6f},{s_:.6f}\n")
print(f"  wrote {out}")

if args.plot is not None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, MUTED, SURF, GRID, BL, OR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"
    t1 = td[-1]; t0_ = t1 - args.plot
    sel = td >= t0_; selm = t >= t0_
    fig, ax = plt.subplots(3, 1, figsize=(15, 11), facecolor=SURF, sharex=True, gridspec_kw=dict(hspace=0.22))
    for A in ax:
        A.set_facecolor(SURF); A.tick_params(colors=MUTED, labelsize=9); A.grid(axis="y", color=GRID, lw=0.6)
        for s in ("top", "right"): A.spines[s].set_visible(False)
        for s in ("left", "bottom"): A.spines[s].set_color(GRID)
    ax[0].plot(td[sel], Fd[sel], color=BL, lw=1.3, label="daily manifold (step on the assumed date)")
    ax[0].plot(t[selm], F[selm], ls="none", marker="o", ms=3, color=INK, label="monthly manifold (lte_results.csv)")
    ax[0].set_title(f"Manifold, impulse assumed on {day0:%d %b}", loc="left", fontsize=11, color=INK); ax[0].legend(frameon=False, fontsize=9, ncol=2)
    ax[1].plot(td[sel], raw[sel], color=OR, lw=0.9, label="daily model, unfiltered")
    ax[1].plot(td[sel], smooth[sel], color=BL, lw=1.6, label=f"daily model, filtered as the fit (F9 {f9})")
    ax[1].set_title("Daily model", loc="left", fontsize=11, color=INK); ax[1].legend(frameon=False, fontsize=9, ncol=2)
    ax[2].plot(t[selm], dat[selm], color=MUTED, lw=1.0, label="data, monthly (as saved)")
    ax[2].plot(t[selm], mod[selm], color=INK, lw=1.0, ls="--", label="model, monthly (column 2)")
    ax[2].plot(td[sel], smooth[sel], color=BL, lw=1.4, label="daily model, filtered")
    ax[2].set_title("Against the monthly model and data", loc="left", fontsize=11, color=INK); ax[2].legend(frameon=False, fontsize=9, ncol=3); ax[2].set_xlabel("Year", color=MUTED)
    fig.suptitle(os.path.basename(os.path.abspath(d)) + ": daily model series", x=0.125, ha="left", fontsize=13, color=INK, y=0.94)
    png = os.path.join(d, "daily_model.png")
    fig.savefig(png, dpi=105, bbox_inches="tight", facecolor=SURF)
    print(f"  plot: {png}")
