#!/usr/bin/env python3
"""Rebuild monthly GLOBAL SST from the 89 grid quads, three stages.

Standalone: numpy + matplotlib + netCDF4 only, no project imports.

  1. DATA.   Column 3 (data) of every quad's lte_results.csv, prorated by
             each quad's ocean area, summed to a global monthly series.
             CC against the real global mean SST.
  2. MODEL.  Column 2 (model), same prorating. CC against the real
             global mean SST.
  3. TREND.  The secular part of each quad's model, from
             lt.exe.windings.json, same prorating: the global linear
             temperature rise and acceleration, plotted over the span.

Prorating weights. The quads are the 20x20 degree boxes that
Feb2026/build_k_sst.py cut out of the Kaplan SST anomaly grid
(k_sst/sst.mean.anom.nc, boxes >= 30% ocean). Each quad's weight is its
own OCEAN area on that grid: sum of cos(lat) over its cells, each cell
counted by the fraction of post-1950 months it has data. So a
half-land coastal box counts about half as much as an open-ocean box
at the same latitude. Weights are renormalized each month over the
quads that have a value that month. Without the .nc file (--nc
missing), falls back to cos(latitude of the box centre).

The real global mean SST ("the data" in stages 1 and 2) is the
cos(lat)-weighted mean of every valid ocean cell on that same Kaplan
grid, i.e. the same source the quads came from, including the
high-latitude/coastal cells no quad covers. Also reported against the
"covered" mean (only cells inside a quad) to separate reconstruction
error from coverage error.

Secular term (stage 3). Each quad's MLR (Regression_Factors) fits
    trend * t + accel * (t - accel_ref)**2          (t = decimal year)
and LTE evaluates the model with the same accel_ref. Plotted per quad,
measured from its record start t0:
    S(t) = trend*(t - t0) + accel*((t - accel_ref)**2 - (t0 - accel_ref)**2)
so its warming rate at year t is trend + 2*accel*(t - accel_ref).

lt.exe.windings.json from builds after 2026-09-29 carries a
"secular_context" block with accel_ref (returned by the MLR itself), the
fit span, and the role of TRAIN_START/TRAIN_END (under EXCLUDE=TRUE they
bound the excluded TEST interval; the fit starts at the record start).
Older JSON files have two known defects, both fixed in the Ada the same
day: trend/accel were saved as abs() values (so a negative MLR accel was
reported positive), and LTE had evaluated accel around TRAIN_START instead
of the MLR's own reference. For those files accel_ref is inferred as the
first date the MLR was fitted on, and they are counted and flagged in the
stage 3 output; re-run the quad with the current lt.exe to replace them.

Curvature budget (stage 3). The lowest winding
amp*sin(2*pi*k*F(t)+phase), the other windings and k0*F(t) also bend
slowly over a 73-year record, so the model's total curvature is shared
among them and accel. Stage 3 prints the area-weighted quadratic
coefficient of each component next to the whole column-2 model and the
column-3 data.

CCs are reported raw AND with a quadratic (trend + acceleration) removed
from both series. The shared warming trend alone makes any two SST
series correlate, so the detrended CC is the honest measure of whether
the month-to-month and interannual structure is reproduced.

Usage:
    python3 global_sst_from_quads.py
    python3 global_sst_from_quads.py --root . --nc ../Feb2026/k_sst/sst.mean.anom.nc \
        --out global_sst_from_quads.png --smooth 12
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_NC = HERE.parent / "Feb2026" / "k_sst" / "sst.mean.anom.nc"
NAME_RE = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")
BOX_DEG = 20.0
MIN_YEAR = 1950.0

# Reference palette (dataviz skill, light mode): categorical slots 1-3,
# ink and chrome.
C_DATA, C_RECON, C_MODEL = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, AXIS, SURFACE = ("#0b0b0b", "#52514e", "#898781",
                                         "#e1e0d9", "#c3c2b7", "#fcfcfb")


# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------

def month_key(date: float) -> int:
    """Project date convention is year + (month-1)/12; key = months since
    year 0, robust to the 1950.083333 vs 1950.0833333 rounding seen
    across files."""
    return int(round(date * 12.0))


def quad_centre(name: str) -> tuple[float, float]:
    m = NAME_RE.match(name)
    lat = float(m.group(2)) * (1 if m.group(1) == "N" else -1)
    lon = float(m.group(4)) * (1 if m.group(3) == "E" else -1)
    return lat, lon


def load_quads(root: Path) -> dict:
    """name -> dict(dates, model (col 2), data (col 3), manifold (col 4),
    windings, t_ref, exclude, pre_fix)."""
    quads = {}
    for d in sorted(root.iterdir()):
        if not (d.is_dir() and NAME_RE.match(d.name)):
            continue
        csv = d / "lte_results.csv"
        if not csv.is_file():
            print(f"  skip {d.name}: no lte_results.csv", file=sys.stderr)
            continue
        a = np.loadtxt(csv, delimiter=",", ndmin=2)
        wj = d / "lt.exe.windings.json"
        windings = json.loads(wj.read_text()) if wj.is_file() else None
        ctx = (windings or {}).get("secular_context")
        if ctx and "accel_ref" in ctx:  # written by lt.exe itself
            refs = dict(t_ref=ctx["accel_ref"], exclude=ctx["exclude"],
                        pre_fix=False)
        else:
            refs = inferred_accel_ref(d, a[:, 0])
        quads[d.name] = dict(dates=a[:, 0], model=a[:, 1], data=a[:, 2],
                             manifold=a[:, 3], windings=windings, **refs)
    return quads


def lt_setting(cell_dir: Path, key: str) -> str | None:
    """A setting as lt.exe resolves it: resp file first, then the
    environment its lte_run.sh exports."""
    resp = cell_dir / "lt.exe.resp"
    if resp.is_file():
        for line in resp.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] == key:
                return parts[1].strip('"')
    sh = cell_dir / "lte_run.sh"
    if sh.is_file():
        m = re.search(rf"^export {key}=['\"]?([^'\"\s]+)", sh.read_text(), re.M)
        if m:
            return m.group(1)
    return None


def inferred_accel_ref(cell_dir: Path, dates: np.ndarray) -> dict:
    """For windings JSON written before secular_context existed: the MLR's
    accel_ref is the first date of the array it was fitted on -- the
    record start under EXCLUDE (TRAIN_START..TRAIN_END is then the
    excluded test gap), otherwise the last record date <= TRAIN_START."""
    exclude = (lt_setting(cell_dir, "EXCLUDE") or "FALSE").upper() == "TRUE"
    ts = lt_setting(cell_dir, "TRAIN_START")
    if exclude or ts is None:
        t_ref = float(dates[0])
    else:
        le = dates[dates <= float(ts) + 1e-6]
        t_ref = float(le[-1]) if len(le) else float(dates[0])
    return dict(t_ref=t_ref, exclude=exclude, pre_fix=True)


def load_kaplan(nc_path: Path):
    """Returns (month keys, anomaly field (t, lat, lon) masked, lat, lon
    in 0..360), post-1950 only."""
    import netCDF4 as nc
    ds = nc.Dataset(nc_path)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    tv = ds.variables["time"]
    cal = nc.num2date(tv[:], tv.units)
    dates = np.array([c.year + (c.month - 1) / 12.0 for c in cal])
    keep = dates >= MIN_YEAR
    sst = ds.variables["sst"][keep]
    return np.array([month_key(x) for x in dates[keep]]), sst, lat, lon


def box_cells(name: str, lat: np.ndarray, lon: np.ndarray):
    """Index arrays of the grid cells inside a quad's box (same edges as
    build_k_sst.extract_regions: lat from -90, lon from 0, 20 deg)."""
    lat_c, lon_c = quad_centre(name)
    lon_c360 = lon_c % 360.0
    lat_sel = np.nonzero((lat >= lat_c - BOX_DEG / 2) & (lat < lat_c + BOX_DEG / 2))[0]
    lon_sel = np.nonzero((lon >= lon_c360 - BOX_DEG / 2) & (lon < lon_c360 + BOX_DEG / 2))[0]
    return lat_sel, lon_sel


def global_mean(sst, lat, cell_mask=None) -> np.ndarray:
    """cos(lat)-weighted mean over valid cells each month (optionally only
    cells where cell_mask is True)."""
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones(sst.shape[2])[None, :]
    if cell_mask is not None:
        w = w * cell_mask
    valid = ~np.ma.getmaskarray(sst)
    vals = sst.filled(0.0)
    num = (vals * w[None] * valid).sum(axis=(1, 2))
    den = (w[None] * valid).sum(axis=(1, 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


# ---------------------------------------------------------------------
# Prorating
# ---------------------------------------------------------------------

def ocean_area_weights(names, nc_field):
    """Per-quad ocean area on the Kaplan grid, plus a mask of all cells
    covered by any quad (for the 'covered' reference mean)."""
    keys, sst, lat, lon = nc_field
    valid_frac = (~np.ma.getmaskarray(sst)).mean(axis=0)       # (lat, lon)
    area = np.cos(np.deg2rad(lat))[:, None] * valid_frac
    covered = np.zeros_like(area, dtype=bool)
    weights = {}
    for n in names:
        ls, os_ = box_cells(n, lat, lon)
        weights[n] = float(area[np.ix_(ls, os_)].sum())
        covered[np.ix_(ls, os_)] = True
    return weights, covered


def coslat_weights(names):
    return {n: float(np.cos(np.deg2rad(quad_centre(n)[0]))) for n in names}


def prorate(quads, weights, column: str, keys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Area-weighted sum of one column across quads on the month grid
    `keys`; renormalized each month over the quads present. Returns
    (series, fraction of total weight present)."""
    num = np.zeros(len(keys))
    den = np.zeros(len(keys))
    pos = {k: i for i, k in enumerate(keys)}
    for n, q in quads.items():
        w = weights[n]
        for date, v in zip(q["dates"], q[column]):
            i = pos.get(month_key(date))
            if i is not None and np.isfinite(v):
                num[i] += w * v
                den[i] += w
    total = sum(weights[n] for n in quads)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan), den / total


def secular_curve(q, t: np.ndarray) -> np.ndarray:
    w = q["windings"]
    t0, tr = float(q["dates"][0]), q["t_ref"]
    return w["trend"] * (t - t0) + w["accel"] * ((t - tr) ** 2 - (t0 - tr) ** 2)


# ---------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------

def quad_detrend(t, y):
    c = np.polyfit(t - t.mean(), y, 2)
    return y - np.polyval(c, t - t.mean())


def running_mean(y, n):
    if n <= 1:
        return y
    k = np.ones(n) / n
    out = np.full_like(y, np.nan)
    half = n // 2
    out[half:len(y) - (n - 1 - half)] = np.convolve(y, k, mode="valid")
    return out


def cc_report(label, t, x, ref):
    ok = np.isfinite(x) & np.isfinite(ref)
    t, x, ref = t[ok], x[ok], ref[ok]
    raw = np.corrcoef(x, ref)[0, 1]
    det = np.corrcoef(quad_detrend(t, x), quad_detrend(t, ref))[0, 1]
    print(f"    {label:34s} CC raw {raw:+.3f}   detrended {det:+.3f}   (n={ok.sum()})")
    return raw, det


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", type=Path, default=HERE,
                    help="directory holding the kNxxx_Exxx quad folders")
    ap.add_argument("--nc", type=Path, default=DEFAULT_NC,
                    help="Kaplan SST anomaly grid the quads were cut from")
    ap.add_argument("--out", type=Path, default=HERE / "global_sst_from_quads.png")
    ap.add_argument("--smooth", type=int, default=12,
                    help="running-mean months for the plot and smoothed CCs")
    args = ap.parse_args()

    quads = load_quads(args.root)
    if not quads:
        ap.error(f"no quads with lte_results.csv under {args.root}")
    names = sorted(quads)
    print(f"{len(names)} quads loaded from {args.root}")

    have_nc = args.nc.is_file()
    if have_nc:
        field = load_kaplan(args.nc)
        keys = field[0]
        weights, covered = ocean_area_weights(names, field)
        glob = global_mean(field[1], field[2])
        glob_cov = global_mean(field[1], field[2], covered)
        print(f"weights: ocean area per quad from {args.nc.name}")
    else:
        print(f"WARNING: {args.nc} not found -- cos(lat) weights, "
              f"no real global SST to compare against", file=sys.stderr)
        weights = coslat_weights(names)
        all_keys = sorted({month_key(d) for q in quads.values() for d in q["dates"]})
        keys, glob, glob_cov = np.array(all_keys), None, None
    t = keys / 12.0

    wsum = sum(weights.values())
    top = sorted(weights.items(), key=lambda kv: -kv[1])
    print(f"  largest weights: " + ", ".join(f"{n} {w / wsum:.3f}" for n, w in top[:3]))
    print(f"  smallest weights: " + ", ".join(f"{n} {w / wsum:.4f}" for n, w in top[-3:]))

    # ---- Stage 1: data (column 3) ------------------------------------
    recon, frac = prorate(quads, weights, "data", keys)
    span = np.isfinite(recon)
    print(f"\nStage 1 -- prorated DATA (col 3), {t[span][0]:.2f}-{t[span][-1]:.2f}, "
          f"min weight present {np.nanmin(frac[span]):.3f}")
    # ---- Stage 2: model (column 2) -----------------------------------
    model, _ = prorate(quads, weights, "model", keys)
    print("Stage 2 -- prorated MODEL (col 2)")
    if glob is not None:
        sm = args.smooth
        for lbl, ref in (("vs global mean SST", glob),
                         ("vs mean over quad-covered cells", glob_cov)):
            print(f"  {lbl}:")
            cc_report("stage 1 data recon, monthly", t, recon, ref)
            cc_report(f"stage 1 data recon, {sm}-mo mean", t,
                      running_mean(recon, sm), running_mean(ref, sm))
            cc_report("stage 2 model, monthly", t, model, ref)
            cc_report(f"stage 2 model, {sm}-mo mean", t,
                      running_mean(model, sm), running_mean(ref, sm))
    print("  model vs data recon (both prorated):")
    cc_report("stage 2 vs stage 1, monthly", t, model, recon)

    # ---- Stage 3: secular trend + acceleration -----------------------
    with_w = [n for n in names if quads[n]["windings"] is not None]
    missing = sorted(set(names) - set(with_w))
    wsum3 = sum(weights[n] for n in with_w)
    tt = t[span]
    sec = sum(weights[n] * secular_curve(quads[n], tt) for n in with_w) / wsum3
    rate = np.gradient(sec, tt)
    g_trend = sum(weights[n] * quads[n]["windings"]["trend"] for n in with_w) / wsum3
    g_accel = sum(weights[n] * quads[n]["windings"]["accel"] for n in with_w) / wsum3
    n_excl = sum(quads[n]["exclude"] for n in with_w)
    pre_fix = [n for n in with_w if quads[n]["pre_fix"]]
    refs = ", ".join(f"{x:g}" for x in sorted({round(quads[n]["t_ref"], 3)
                                               for n in with_w}))
    print(f"\nStage 3 -- prorated secular term from lt.exe.windings.json "
          f"({len(with_w)} quads, {n_excl} with EXCLUDE=TRUE"
          f"{'; missing: ' + ' '.join(missing) if missing else ''})")
    if pre_fix:
        print(f"  NOTE: {len(pre_fix)}/{len(with_w)} quads have pre-fix JSON "
              f"(no secular_context): accel_ref inferred, and trend/accel "
              f"were saved as abs() so any negative MLR value shows as "
              f"positive. Re-run them with the current lt.exe.")
    print(f"  accel_ref: {refs}")
    print(f"  area-weighted trend coefficient  {g_trend:+.5f} C/yr "
          f"({10 * g_trend:+.3f} C/decade)")
    print(f"  area-weighted accel coefficient  {g_accel:+.3e} C/yr^2 "
          f"(curvature 2*accel = {2 * g_accel:+.3e} C/yr^2)")
    print(f"  warming rate {10 * rate[0]:+.3f} C/decade at {tt[0]:.0f} -> "
          f"{10 * rate[-1]:+.3f} at {tt[-1]:.0f}; secular rise "
          f"{sec[-1] - sec[0]:+.3f} C")

    curvature_budget(quads, weights, with_w, t, span)

    plot(args, t, span, glob, recon, model, tt, sec, rate)
    print(f"\nplot: {args.out}")
    return 0


def curvature_budget(quads, weights, names, t, span):
    """Area-weighted quadratic coefficient (C/yr^2) and 1950-end rise of
    each additive model component, on the common month grid."""
    tt = t[span]
    pos = {month_key(x): i for i, x in enumerate(tt)}
    comps = {k: np.zeros(len(tt)) for k in
             ("accel term", "lowest winding", "other windings",
              "k0*F(t) (manifold)", "model total (col 2)", "data (col 3)")}
    wsum = np.zeros(len(tt))
    for n in names:
        q, w = quads[n], weights[n]
        idx = np.array([pos.get(month_key(x), -1) for x in q["dates"]])
        ok = idx >= 0
        wd = q["windings"]
        d, F = q["dates"][ok], q["manifold"][ok]
        kap = np.array(wd["k_amp_phase"])
        terms = [a * np.sin(2 * np.pi * k * F + ph) for k, a, ph in kap]
        lo = int(np.argmin(np.abs(kap[:, 0])))
        parts = {
            "accel term": wd["accel"] * (d - q["t_ref"]) ** 2,
            "lowest winding": terms[lo],
            "other windings": sum(terms) - terms[lo],
            "k0*F(t) (manifold)": wd["k0"] * F,
            "model total (col 2)": q["model"][ok],
            "data (col 3)": q["data"][ok],
        }
        for k, y in parts.items():
            comps[k][idx[ok]] += w * y
        wsum[idx[ok]] += w
    x = tt - tt.mean()
    total = None
    print("  curvature budget (area-weighted quadratic fit to each component):")
    print(f"    {'component':26s} {'quad coef C/yr^2':>17s} {'share of model':>15s}")
    rows = []
    for k, y in comps.items():
        c2 = np.polyfit(x, y / wsum, 2)[0]
        rows.append((k, c2))
        if k == "model total (col 2)":
            total = c2
    for k, c2 in rows:
        share = f"{100 * c2 / total:13.0f} %" if "total" not in k and "data" not in k else ""
        print(f"    {k:26s} {c2:+17.2e} {share:>15s}")


def plot(args, t, span, glob, recon, model, tt, sec, rate):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False,
    })
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                 gridspec_kw=dict(height_ratios=[3, 2], hspace=0.12))
    sm = args.smooth

    if glob is not None:
        a1.plot(t, running_mean(glob, sm), color=C_DATA, lw=2,
                label="Global mean SST (Kaplan, all ocean cells)")
    a1.plot(t, running_mean(recon, sm), color=C_RECON, lw=2,
            label="Stage 1: quads' data (col 3), area-prorated")
    a1.plot(t, running_mean(model, sm), color=C_MODEL, lw=2,
            label="Stage 2: quads' model (col 2), area-prorated")
    a1.set_ylabel(f"SST anomaly, °C ({sm}-month mean)")
    a1.set_title("Global SST rebuilt from the grid quads", loc="left",
                 fontsize=12, color=INK)
    a1.legend(loc="upper left", fontsize=9)

    a2.plot(tt, sec, color=C_MODEL, lw=2)
    a2.set_title("Stage 3: area-prorated secular term (trend + accel) "
                 "from the quads' MLR", loc="left", fontsize=10, color=INK2)
    a2.set_ylabel(f"Rise since {tt[0]:.0f}, °C")
    a2.set_xlabel("Year")
    a2.legend(loc="upper left", fontsize=9)
    a2.annotate(f"{10 * rate[0]:+.2f} °C/decade", (tt[0], sec[0]),
                xytext=(8, 10), textcoords="offset points", color=INK2, fontsize=9)
    a2.annotate(f"{10 * rate[-1]:+.2f} °C/decade", (tt[-1], sec[-1]),
                xytext=(-6, 6), textcoords="offset points", ha="right",
                color=INK2, fontsize=9)
    a2.set_xlim(tt[0], tt[-1])
    fig.savefig(args.out, dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
