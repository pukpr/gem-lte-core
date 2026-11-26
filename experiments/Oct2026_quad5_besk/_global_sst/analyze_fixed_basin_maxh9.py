#!/usr/bin/env python3
"""Analyze the 89-quad fixed-basin lpap, MAXH=9, NH=5, ACCEL=FALSE run (1880-2023, raw monthly)
and compare directly against the previous per-quad refined (MAXH=30, ACCEL=TRUE) 89-quad run."""
import math
import sys
from pathlib import Path

import matplotlib
import netCDF4 as nc4
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
RAW = EXP / "Oct2026_quad5_refine/_global_sst/raw"
sys.path.insert(0, str(EXP / "Sep2026"))
import global_sst_from_quads as G  # noqa: E402


def load_kaplan_from(nc_path: Path, min_year: float = 1880.0):
    ds = nc4.Dataset(nc_path)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    tv = ds.variables["time"]
    cal = nc4.num2date(tv[:], tv.units)
    dates = np.array([c.year + (c.month - 1) / 12.0 for c in cal])
    keep = dates >= min_year
    sst = ds.variables["sst"][keep]
    return np.array([G.month_key(x) for x in dates[keep]]), sst, lat, lon


def rebase_1961_1990(d):
    out = dict(d)
    for m in range(12):
        base = [
            v for k, v in d.items()
            if k % 12 == m and 1961 * 12 <= k < 1991 * 12 and np.isfinite(v)
        ]
        mu = float(np.mean(base))
        for k in d:
            if k % 12 == m:
                out[k] = d[k] - mu
    return out


def detr_lin(t, y):
    return y - np.polyval(np.polyfit(t, y, 1), t)


def detr_quad(t, y):
    tc = t - t.mean()
    return y - np.polyval(np.polyfit(tc, y, 2), tc)


def cc_rms(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    return float(np.corrcoef(a, b)[0, 1]), float(np.sqrt(np.mean((a - b) ** 2)))


def oob_significance(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok] - np.mean(x[ok]), y[ok] - np.mean(y[ok])
    n = len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    r1x = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    r1y = float(np.corrcoef(y[:-1], y[1:])[0, 1])
    n_eff = n * (1.0 - r1x * r1y) / max(1.0 + r1x * r1y, 1e-6)
    n_eff = max(3.0, min(float(n), n_eff))
    t_stat = r * math.sqrt((n_eff - 2.0) / max(1.0 - r * r, 1e-12))
    p = math.erfc(abs(t_stat) / math.sqrt(2.0))
    return r, n_eff, t_stat, p


def bp(y, lo_mo=None, hi_mo=None):
    F = np.fft.rfft(y - np.mean(y))
    freqs = np.fft.rfftfreq(len(y), d=1 / 12.0)
    keep = np.ones_like(freqs, dtype=bool)
    if lo_mo:
        keep &= freqs <= 12.0 / lo_mo
    if hi_mo:
        keep &= freqs >= 12.0 / hi_mo
    return np.fft.irfft(F * keep, n=len(y))


def main():
    quads = G.load_quads(HERE / "fixed_basin_maxh9")
    nc_1950 = G.load_kaplan(EXP / "Feb2026/k_sst/sst.mean.anom.nc")
    weights, _ = G.ocean_area_weights(sorted(quads), nc_1950)

    kk_full, sst_full, lat_full, _ = load_kaplan_from(EXP / "Feb2026/k_sst/sst.mean.anom.nc", 1880.0)
    keys = np.arange(
        min(G.month_key(q["dates"][0]) for q in quads.values()),
        max(G.month_key(q["dates"][-1]) for q in quads.values()) + 1,
    )
    mod, _ = G.prorate(quads, weights, "model", keys)

    kap = dict(zip(kk_full, G.global_mean(sst_full, lat_full)))
    had = {
        int(r[0]) * 12 + int(r[1]) - 1: float(r[2])
        for r in np.genfromtxt(RAW / "hadsst4_globe.csv", delimiter=",", skip_header=1)
    }
    noaa = {}
    for line in open(RAW / "noaa_ocean.csv"):
        x = line.strip().split(",")
        if len(x) == 2 and x[0].isdigit():
            noaa[int(x[0][:4]) * 12 + int(x[0][4:]) - 1] = float(x[1])

    S = {
        "model": rebase_1961_1990(dict(zip(keys, mod))),
        "Kaplan": rebase_1961_1990(kap),
        "HadSST": rebase_1961_1990(had),
        "NOAA": rebase_1961_1990(noaa),
    }
    obs = ("Kaplan", "HadSST", "NOAA")
    common = [
        k for k in keys
        if all(k in S[o] and np.isfinite(S[o][k]) for o in obs)
        and np.isfinite(S["model"][k])
    ]
    S["consensus"] = {k: float(np.mean([S[o][k] for o in obs])) for k in common}

    t = np.array(common) / 12.0
    V = {n: np.array([S[n][k] for k in common]) for n in S}
    ref = V["consensus"]
    m = V["model"]
    fit_m = t >= 1950.0
    oob_m = t < 1950.0

    # Per-quad std & spatial coherence
    pq_fit, pq_oob = [], []
    for q in quads.values():
        tq, mq = q["dates"], q["model"]
        f_q, o_q = tq >= 1950.0, tq < 1950.0
        pq_fit.append(np.std(detr_lin(tq[f_q], mq[f_q])))
        pq_oob.append(np.std(detr_lin(tq[o_q], mq[o_q])))
    g_fit = np.std(detr_lin(t[fit_m], m[fit_m]))
    g_oob = np.std(detr_lin(t[oob_m], m[oob_m]))
    coh_fit = (g_fit ** 2) / (np.mean(pq_fit) ** 2)
    coh_oob = (g_oob ** 2) / (np.mean(pq_oob) ** 2)

    lines = [
        f"89 quads, FIXED basin lpap, MAXH=9, NH=5, ACCEL=FALSE (raw MONTHLY, no 12-mo mean), {t[0]:.0f}-{t[-1]:.2f}",
        f"  Per-quad detrended std : 1950-2023={np.mean(pq_fit):.4f} C, 1880-1949={np.mean(pq_oob):.4f} C (ratio={np.mean(pq_oob)/np.mean(pq_fit):.3f})",
        f"  Global-mean detr std   : 1950-2023={g_fit:.4f} C, 1880-1949={g_oob:.4f} C (ratio={g_oob/g_fit:.3f})",
        f"  Spatial coherence      : 1950-2023={coh_fit:.4f}, 1880-1949={coh_oob:.4f} (drop={coh_fit/coh_oob:.2f}x)",
        "",
        "=== 1. RAW MONTHLY COMPARISON VS 3-DATASET CONSENSUS ===",
        f"{'window':16s} {'n':>5s} | {'raw CC / RMS':>16s} | {'detr(lin) CC':>13s} | {'detr(quad) CC':>14s} | {'OOB p (N_eff)':>18s}",
    ]
    windows = [
        ("1950-2023 (fit)", fit_m),
        ("1880-1949 (OOB)", oob_m),
        ("1880-1914 (OOB1)", (t >= 1880.0) & (t < 1915.0)),
        ("1915-1949 (OOB2)", (t >= 1915.0) & (t < 1950.0)),
        ("1920-1949 (OOB)", (t >= 1920.0) & (t < 1950.0)),
        ("1930-1949 (OOB)", (t >= 1930.0) & (t < 1950.0)),
        ("1880-2023 (all)", np.ones_like(t, dtype=bool)),
    ]
    for lbl, mask in windows:
        c1, r1 = cc_rms(m[mask], ref[mask])
        cd1, _ = cc_rms(detr_lin(t[mask], m[mask]), detr_lin(t[mask], ref[mask]))
        cd2, _ = cc_rms(detr_quad(t[mask], m[mask]), detr_quad(t[mask], ref[mask]))
        _, neff, _, pval = oob_significance(detr_lin(t[mask], m[mask]), detr_lin(t[mask], ref[mask]))
        lines.append(
            f"{lbl:16s} {int(mask.sum()):5d} | {c1:+.3f} / {r1:.3f} C | "
            f"{cd1:+.3f}         | {cd2:+.3f}          | p={pval:.2e} ({neff:.0f})"
        )

    lines += ["", "=== 2. FREQUENCY-BAND BREAKDOWN (DETRENDED RAW MONTHLY) ==="]
    for band, lo, hi in [("fast (<2 yr)", None, 24), ("ENSO (2-7 yr)", 24, 84), ("decadal (>7 yr)", 84, None)]:
        mf, cf = bp(detr_lin(t[fit_m], m[fit_m]), lo, hi), bp(detr_lin(t[fit_m], ref[fit_m]), lo, hi)
        mo, co = bp(detr_lin(t[oob_m], m[oob_m]), lo, hi), bp(detr_lin(t[oob_m], ref[oob_m]), lo, hi)
        lines.append(
            f"  {band:15s}: FIT CC={np.corrcoef(mf, cf)[0,1]:+.3f} (std {np.std(mf):.3f}/{np.std(cf):.3f}) | "
            f"OOB CC={np.corrcoef(mo, co)[0,1]:+.3f} (std {np.std(mo):.3f}/{np.std(co):.3f})"
        )

    lines += [
        "",
        "=== 3. LINEAR TRENDS (C/century) ===",
        f"  1950-2023: consensus {np.polyfit(t[fit_m], ref[fit_m], 1)[0]*100:+.2f} | model {np.polyfit(t[fit_m], m[fit_m], 1)[0]*100:+.2f} | Kaplan {np.polyfit(t[fit_m], V['Kaplan'][fit_m], 1)[0]*100:+.2f}",
        f"  1880-1949: consensus {np.polyfit(t[oob_m], ref[oob_m], 1)[0]*100:+.2f} | model {np.polyfit(t[oob_m], m[oob_m], 1)[0]*100:+.2f} | Kaplan {np.polyfit(t[oob_m], V['Kaplan'][oob_m], 1)[0]*100:+.2f}",
    ]

    print("\n".join(lines))
    with open(HERE / "global_sst_fixed_basin_maxh9.txt", "w") as fh:
        fh.write("# " + "\n# ".join(lines) + "\n")
        fh.write("# date model consensus Kaplan HadSST NOAA\n")
        for i, k in enumerate(common):
            fh.write(
                f"{k / 12:.4f} "
                + " ".join(f"{V[n][i]:+.4f}" for n in ("model", "consensus", "Kaplan", "HadSST", "NOAA"))
                + "\n"
            )

    # 3-panel raw monthly figure
    INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    fig, ax = plt.subplots(3, 1, figsize=(15.5, 11.5), facecolor=SURF, gridspec_kw=dict(height_ratios=[2.5, 2.5, 2.5]))
    for A in ax:
        A.set_facecolor(SURF)
        A.grid(color=GRID, lw=0.6)
        for s in ("top", "right"):
            A.spines[s].set_visible(False)

    # Panel 1: 1950-2023 raw monthly
    A = ax[0]
    A.plot(t[fit_m], ref[fit_m], color=INK, lw=1.0, alpha=0.85, label="consensus, raw monthly (Kaplan, HadSST, NOAA)")
    A.plot(t[fit_m], m[fit_m], color="#eb6834", lw=1.0, alpha=0.9, label="89-quad FIXED-basin MAXH=9 ACCEL=FALSE, raw monthly")
    c_f, r_f = cc_rms(m[fit_m], ref[fit_m])
    cd_f, _ = cc_rms(detr_lin(t[fit_m], m[fit_m]), detr_lin(t[fit_m], ref[fit_m]))
    A.set_title(
        f"1950-2023 Fit Span (Fixed Basin lpap, MAXH=9, ACCEL=FALSE, Raw Monthly): CC = {c_f:+.3f} (RMS {r_f:.3f} C, detrended {cd_f:+.3f})",
        loc="left", fontsize=11, color=INK,
    )
    A.set_ylabel("anomaly vs 1961-1990 (C)")
    A.set_xlim(1950, 2023)
    A.legend(frameon=False, ncol=2, fontsize=9, loc="upper left")

    # Panel 2: 1880-1949 OOB backward validation (detrended raw monthly)
    A = ax[1]
    t_o = t[oob_m]
    ref_od = detr_lin(t_o, ref[oob_m])
    mod_od = detr_lin(t_o, m[oob_m])
    c_od, r_od = cc_rms(mod_od, ref_od)
    _, neff_o, _, p_o = oob_significance(mod_od, ref_od)
    A.plot(t_o, ref_od, color=INK, lw=1.0, alpha=0.85, label="consensus, raw monthly (detrended on 1880-1949)")
    A.plot(t_o, mod_od, color="#eb6834", lw=1.0, alpha=0.9, label="89-quad FIXED-basin MAXH=9 OOB model, raw monthly (detrended)")
    A.axhline(0, color=MUTED, lw=0.7)
    A.set_title(
        f"1880-1949 Non-Fit OOB Backward Validation (Raw Monthly, Detrended): CC = {c_od:+.3f} (RMS {r_od:.3f} C, p = {p_o:.2e}, N_eff = {neff_o:.0f})",
        loc="left", fontsize=11, color=INK,
    )
    A.set_ylabel("detrended monthly (C)")
    A.set_xlim(1880, 1950)
    A.legend(frameon=False, ncol=2, fontsize=9, loc="upper left")

    # Panel 3: Full 1880-2023 raw monthly overview (no parabola artifact at all!)
    A = ax[2]
    A.axvspan(1880, 1950, color="#eceae3", alpha=0.7, lw=0, label="1880-1949 OOB (zero training data)")
    A.axvline(1950, color=MUTED, lw=1.0, ls=":")
    A.plot(t, ref, color=INK, lw=0.85, alpha=0.85, label="consensus, raw monthly")
    A.plot(t, m, color="#eb6834", lw=0.85, alpha=0.9, label="89-quad FIXED-basin MAXH=9 ACCEL=FALSE (unadjusted)")
    c_all, r_all = cc_rms(m, ref)
    A.set_title(
        f"Full 1880-2023 Raw Monthly Global SST (Unadjusted, No Quadratic Bend): CC = {c_all:+.3f} (RMS {r_all:.3f} C)",
        loc="left", fontsize=11, color=INK,
    )
    A.set_ylabel("anomaly vs 1961-1990 (C)")
    A.set_xlabel("Year")
    A.set_xlim(1880, 2023)
    A.legend(frameon=False, ncol=3, fontsize=8.8, loc="upper left")

    fig.tight_layout()
    fig.savefig(HERE / "global_sst_fixed_basin_maxh9.png", dpi=130, facecolor=SURF)
    print(f"plot: {HERE / 'global_sst_fixed_basin_maxh9.png'}")


if __name__ == "__main__":
    main()
