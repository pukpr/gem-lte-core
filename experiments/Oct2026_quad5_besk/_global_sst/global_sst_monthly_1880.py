#!/usr/bin/env python3
"""Raw monthly (no 12-month smoothing) global SST comparison (1950-2023) AND
70-year non-fit OOB backward extrapolation (1880-1949) from all 89 fitted quads.

Loads:
  - extrap1880/<quad>/lte_results.csv (native lt.exe F9=-1, ACCEL=TRUE, EXCLUDE=FALSE,
    TRAIN_START=1950.0, TRAIN_END=2023.0; zero pre-1950 training data)
  - Kaplan SST V2 all-ocean mean (1856-2023)
  - HadSST 4.2.0.0 GLOBE (1850-2026)
  - NOAA NCEI ERSST v5 global ocean (1850-2026)

Writes global_sst_monthly_1880.txt and global_sst_monthly_1880.png.
"""
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
    """Autocorrelation-adjusted two-sided p-value (Bretherton et al. 1999) for monthly series."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok] - np.mean(x[ok]), y[ok] - np.mean(y[ok])
    n = len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    r1x = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    r1y = float(np.corrcoef(y[:-1], y[1:])[0, 1])
    n_eff = n * (1.0 - r1x * r1y) / max(1.0 + r1x * r1y, 1e-6)
    n_eff = max(3.0, min(float(n), n_eff))
    t_stat = r * math.sqrt((n_eff - 2.0) / max(1.0 - r * r, 1e-12))
    # Normal approx for df >= 30 (exact enough for reporting)
    z = abs(t_stat)
    p = math.erfc(z / math.sqrt(2.0))
    return r, n_eff, t_stat, p


def main():
    quads = G.load_quads(HERE / "extrap1880")
    nc_1950 = G.load_kaplan(EXP / "Feb2026/k_sst/sst.mean.anom.nc")
    weights, _ = G.ocean_area_weights(sorted(quads), nc_1950)

    kk_full, sst_full, lat_full, _ = load_kaplan_from(EXP / "Feb2026/k_sst/sst.mean.anom.nc", 1880.0)
    keys = np.arange(
        min(G.month_key(q["dates"][0]) for q in quads.values()),
        max(G.month_key(q["dates"][-1]) for q in quads.values()) + 1,
    )

    # 1. Raw native lt.exe 1880-2023 model (includes accel*(t-1950)^2 backward before 1950)
    mod_raw, _ = G.prorate(quads, weights, "model", keys)

    # 2. Same model with post-1950 quadratic accel*(t-1950)^2 NOT extrapolated backward before 1950
    wsum = sum(weights[n] for n in quads)
    t_keys = keys / 12.0
    g_accel = sum(weights[n] * quads[n]["windings"]["accel"] for n in quads) / wsum
    pre_mask = t_keys < 1950.0
    parabola_pre = np.where(pre_mask, g_accel * (t_keys - 1950.0) ** 2, 0.0)
    mod_no_parab = mod_raw - parabola_pre

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
        "model_raw": rebase_1961_1990(dict(zip(keys, mod_raw))),
        "model_noparab": rebase_1961_1990(dict(zip(keys, mod_no_parab))),
        "Kaplan": rebase_1961_1990(kap),
        "HadSST": rebase_1961_1990(had),
        "NOAA": rebase_1961_1990(noaa),
    }
    obs = ("Kaplan", "HadSST", "NOAA")
    common = [
        k for k in keys
        if all(k in S[o] and np.isfinite(S[o][k]) for o in obs)
        and np.isfinite(S["model_raw"][k])
    ]
    S["consensus"] = {k: float(np.mean([S[o][k] for o in obs])) for k in common}

    t = np.array(common) / 12.0
    V = {n: np.array([S[n][k] for k in common]) for n in S}
    ref = V["consensus"]
    m_raw = V["model_raw"]
    m_nop = V["model_noparab"]

    # Descriptive linear trend gap fitted ONLY on 1950-2023
    fit_m = t >= 1950.0
    oob_m = t < 1950.0
    tc_fit = t[fit_m] - t[fit_m].mean()
    p1_gap = np.polyfit(tc_fit, ref[fit_m] - m_raw[fit_m], 1)
    m_fit_lingap = m_nop + np.polyval(p1_gap, t - t[fit_m].mean())

    lines = [
        f"89 basin-seeded quads, raw MONTHLY (no 12-mo smoothing), {t[0]:.0f}-{t[-1]:.2f} ({len(t)} months)",
        f"  Fit span : 1950-2023 ({int(fit_m.sum())} months) | Non-fit OOB span: 1880-1949 ({int(oob_m.sum())} months)",
        f"  Area-weighted quad accel = {g_accel:+.3e} C/yr^2 (parabola at 1880 = {g_accel * 70**2:+.3f} C)",
        "",
        "=== 1. RAW MONTHLY COMPARISON (NO 12-MONTH MEAN) VS 3-DATASET CONSENSUS ===",
        f"{'window':16s} {'n':>5s} | {'raw_mod CC/RMS':>16s} | {'noparab CC/RMS':>16s} | {'detr(lin) CC':>13s} | {'detr(quad) CC':>14s} | {'OOB p (N_eff)':>18s}",
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
        c1, r1 = cc_rms(m_raw[mask], ref[mask])
        c2, r2 = cc_rms(m_nop[mask], ref[mask])
        cd1, _ = cc_rms(detr_lin(t[mask], m_nop[mask]), detr_lin(t[mask], ref[mask]))
        cd2, _ = cc_rms(detr_quad(t[mask], m_raw[mask]), detr_quad(t[mask], ref[mask]))
        _, neff, _, pval = oob_significance(detr_lin(t[mask], m_nop[mask]), detr_lin(t[mask], ref[mask]))
        lines.append(
            f"{lbl:16s} {int(mask.sum()):5d} | {c1:+.3f} / {r1:.3f} C | {c2:+.3f} / {r2:.3f} C | "
            f"{cd1:+.3f}         | {cd2:+.3f}          | p={pval:.2e} ({neff:.0f})"
        )

    # Also compare against each individual dataset on 1880-1949 OOB
    lines += ["", "=== 2. 1880-1949 OOB MONTHLY CC BY INDIVIDUAL OBSERVED DATASET ==="]
    for o in ("consensus", "Kaplan", "HadSST", "NOAA"):
        c_nop, r_nop = cc_rms(m_nop[oob_m], V[o][oob_m])
        c_d1, _ = cc_rms(detr_lin(t[oob_m], m_nop[oob_m]), detr_lin(t[oob_m], V[o][oob_m]))
        c_d2, _ = cc_rms(detr_quad(t[oob_m], m_raw[oob_m]), detr_quad(t[oob_m], V[o][oob_m]))
        _, neff, _, pval = oob_significance(detr_lin(t[oob_m], m_nop[oob_m]), detr_lin(t[oob_m], V[o][oob_m]))
        lines.append(
            f"  vs {o:10s}: noparab {c_nop:+.3f} (RMS {r_nop:.3f} C) | detr(lin) {c_d1:+.3f} (p={pval:.2e}, N_eff={neff:.0f}) | detr(quad) {c_d2:+.3f}"
        )

    # Check inter-dataset agreement on 1880-1949 (noise floor of early SST observations)
    lines += ["", "=== 3. OBSERVATIONAL NOISE FLOOR ON 1880-1949 (INTER-DATASET AGREEMENT) ==="]
    for a, b in (("Kaplan", "HadSST"), ("Kaplan", "NOAA"), ("HadSST", "NOAA")):
        c_ab, r_ab = cc_rms(V[a][oob_m], V[b][oob_m])
        cd_ab, _ = cc_rms(detr_lin(t[oob_m], V[a][oob_m]), detr_lin(t[oob_m], V[b][oob_m]))
        lines.append(f"  {a:6s} vs {b:6s}: monthly CC={c_ab:+.3f} (RMS {r_ab:.3f} C) | detr(lin) CC={cd_ab:+.3f}")

    # Linear trends on 1880-1949 vs 1950-2023
    lines += [
        "",
        "=== 4. LINEAR TRENDS (C/century) ===",
        f"  1950-2023: consensus {np.polyfit(t[fit_m], ref[fit_m], 1)[0]*100:+.2f} | model {np.polyfit(t[fit_m], m_raw[fit_m], 1)[0]*100:+.2f} | Kaplan {np.polyfit(t[fit_m], V['Kaplan'][fit_m], 1)[0]*100:+.2f}",
        f"  1880-1949: consensus {np.polyfit(t[oob_m], ref[oob_m], 1)[0]*100:+.2f} | model_raw {np.polyfit(t[oob_m], m_raw[oob_m], 1)[0]*100:+.2f} | model_noparab {np.polyfit(t[oob_m], m_nop[oob_m], 1)[0]*100:+.2f} | Kaplan {np.polyfit(t[oob_m], V['Kaplan'][oob_m], 1)[0]*100:+.2f}",
    ]

    print("\n".join(lines))
    with open(HERE / "global_sst_monthly_1880.txt", "w") as fh:
        fh.write("# " + "\n# ".join(lines) + "\n")
        fh.write("# date model_raw model_noparab consensus Kaplan HadSST NOAA\n")
        for i, k in enumerate(common):
            fh.write(
                f"{k / 12:.4f} "
                + " ".join(f"{V[n][i]:+.4f}" for n in ("model_raw", "model_noparab", "consensus", "Kaplan", "HadSST", "NOAA"))
                + "\n"
            )

    # 3-panel raw monthly figure (NO 12-month smoothing)
    INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    fig, ax = plt.subplots(3, 1, figsize=(15.5, 11.5), facecolor=SURF, gridspec_kw=dict(height_ratios=[2.5, 2.5, 2.5]))
    for A in ax:
        A.set_facecolor(SURF)
        A.grid(color=GRID, lw=0.6)
        for s in ("top", "right"):
            A.spines[s].set_visible(False)

    # Panel 1: 1950-2023 raw monthly (no 12-mo mean)
    A = ax[0]
    A.plot(t[fit_m], ref[fit_m], color=INK, lw=1.0, alpha=0.85, label="consensus, raw monthly (Kaplan, HadSST, NOAA)")
    A.plot(t[fit_m], m_raw[fit_m], color="#eb6834", lw=1.0, alpha=0.9, label="89-quad model, raw monthly (unadjusted)")
    A.plot(t[fit_m], m_fit_lingap[fit_m], color="#1baf7a", lw=0.9, alpha=0.75, ls="--",
           label=f"model + descriptive linear gap ({p1_gap[0]*100:+.2f} C/cy)")
    c_f, r_f = cc_rms(m_raw[fit_m], ref[fit_m])
    c_fg, r_fg = cc_rms(m_fit_lingap[fit_m], ref[fit_m])
    A.set_title(
        f"1950-2023 Fit Span (Raw Monthly, No 12-Mo Smoothing): CC = {c_f:+.3f} (RMS {r_f:.3f} C; with linear gap {c_fg:+.3f} / {r_fg:.3f} C)",
        loc="left", fontsize=11, color=INK,
    )
    A.set_ylabel("anomaly vs 1961-1990 (C)")
    A.set_xlim(1950, 2023)
    A.legend(frameon=False, ncol=3, fontsize=8.8, loc="upper left")

    # Panel 2: 1880-1949 OOB backward validation (raw monthly, linearly detrended)
    A = ax[1]
    t_o = t[oob_m]
    ref_od = detr_lin(t_o, ref[oob_m])
    mod_od = detr_lin(t_o, m_nop[oob_m])
    c_od, r_od = cc_rms(mod_od, ref_od)
    _, neff_o, _, p_o = oob_significance(mod_od, ref_od)
    A.plot(t_o, ref_od, color=INK, lw=1.0, alpha=0.85, label="consensus, raw monthly (detrended on 1880-1949)")
    A.plot(t_o, mod_od, color="#eb6834", lw=1.0, alpha=0.9, label="89-quad OOB model, raw monthly (detrended)")
    A.axhline(0, color=MUTED, lw=0.7)
    A.set_title(
        f"1880-1949 Non-Fit OOB Backward Validation (Raw Monthly, Detrended): CC = {c_od:+.3f} (RMS {r_od:.3f} C, p = {p_o:.2e}, N_eff = {neff_o:.0f})",
        loc="left", fontsize=11, color=INK,
    )
    A.set_ylabel("detrended monthly (C)")
    A.set_xlim(1880, 1950)
    A.legend(frameon=False, ncol=2, fontsize=9, loc="upper left")

    # Panel 3: Full 1880-2023 raw monthly overview
    A = ax[2]
    A.axvspan(1880, 1950, color="#eceae3", alpha=0.7, lw=0, label="1880-1949 OOB (zero training data)")
    A.axvline(1950, color=MUTED, lw=1.0, ls=":")
    A.plot(t, ref, color=INK, lw=0.85, alpha=0.85, label="consensus, raw monthly")
    A.plot(t, m_nop, color="#eb6834", lw=0.85, alpha=0.9, label="89-quad model, raw monthly (accel=0 pre-1950)")
    A.plot(t[oob_m], m_raw[oob_m], color="#8b5cf6", lw=0.7, alpha=0.55, ls="--",
           label=f"raw parabola extrap (accel={g_accel:+.1e} C/yr^2)")
    A.set_title("Full 1880-2023 Raw Monthly Global SST: 70-Year OOB Backward Extrapolation + 73-Year Fit", loc="left", fontsize=11, color=INK)
    A.set_ylabel("anomaly vs 1961-1990 (C)")
    A.set_xlabel("Year")
    A.set_xlim(1880, 2023)
    A.legend(frameon=False, ncol=3, fontsize=8.5, loc="upper left")

    fig.tight_layout()
    fig.savefig(HERE / "global_sst_monthly_1880.png", dpi=130, facecolor=SURF)
    print(f"plot: {HERE / 'global_sst_monthly_1880.png'}")


if __name__ == "__main__":
    main()
