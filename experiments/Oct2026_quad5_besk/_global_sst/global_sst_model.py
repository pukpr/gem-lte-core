#!/usr/bin/env python3
"""Monthly global SST from the 30 basin-seeded fast-refine Kaplan quads (Oct2026_quad5_besk).

Compares:
  - model      : ocean-area-weighted mean of the 30 fitted quads (rescored unfiltered F9=-1,
                 with each quad's own fitted trend + accel; no post-hoc rescaling)
  - quad data  : ocean-area-weighted mean of the same 30 quads' observed SST
  - Kaplan     : cos(lat)-weighted mean of all valid ocean cells on the Kaplan SST V2 grid
  - HadSST     : Met Office HadSST 4.2.0.0 monthly GLOBE
  - NOAA       : NOAA NCEI NOAAGlobalTemp / ERSST v5 global ocean
  - consensus  : mean of Kaplan, HadSST, and NOAA

Any missing trend between the model and observed global SST is decomposed and treated
descriptively into:
  (1) model vs 30-quad data trend gap,
  (2) 30-quad spatial sampling vs full-ocean Kaplan trend gap, and
  (3) Kaplan vs 3-dataset consensus trend gap.

Writes global_sst_model.txt and global_sst_model.png.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
RAW = EXP / "Oct2026_quad5_refine/_global_sst/raw"
sys.path.insert(0, str(EXP / "Sep2026"))
import global_sst_from_quads as G  # noqa: E402


def rebase(d):
    """Remove each series' own 1961-1990 calendar-month means."""
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


def run12(y):
    out = np.full_like(y, np.nan)
    k = np.ones(12) / 12.0
    out[6:len(y) - 5] = np.convolve(y, k, mode="valid")
    return out


def score(a, b, m):
    ok = m & np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]), float(np.sqrt(np.mean((a[ok] - b[ok]) ** 2)))


def slope_c_per_century(t, y, m=None):
    if m is not None:
        t, y = t[m], y[m]
    return float(np.polyfit(t, y, 1)[0] * 100.0)


def main():
    quads = G.load_quads(HERE / "nofilt")
    nc = G.load_kaplan(EXP / "Feb2026/k_sst/sst.mean.anom.nc")
    weights, covered = G.ocean_area_weights(sorted(quads), nc)
    all_keys_nc, sst_nc, lat_nc, _ = nc
    total_ocean_area = float((np.cos(np.deg2rad(lat_nc))[:, None] * (~np.ma.getmaskarray(sst_nc)).mean(axis=0)).sum())
    cov_frac = sum(weights.values()) / total_ocean_area

    keys = np.arange(
        min(G.month_key(q["dates"][0]) for q in quads.values()),
        max(G.month_key(q["dates"][-1]) for q in quads.values()) + 1,
    )
    model, frac = G.prorate(quads, weights, "model", keys)
    recon, _ = G.prorate(quads, weights, "data", keys)

    kap = dict(zip(all_keys_nc, G.global_mean(sst_nc, lat_nc)))
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
        "model": rebase(dict(zip(keys, model))),
        "quad data": rebase(dict(zip(keys, recon))),
        "Kaplan": rebase(kap),
        "HadSST": rebase(had),
        "NOAA": rebase(noaa),
    }
    obs = ("Kaplan", "HadSST", "NOAA")
    common = [
        k for k in keys
        if all(k in S[o] and np.isfinite(S[o][k]) for o in obs)
        and np.isfinite(S["model"][k])
        and np.isfinite(S["quad data"][k])
    ]
    S["consensus"] = {k: float(np.mean([S[o][k] for o in obs])) for k in common}

    t = np.array(common) / 12.0
    held = (t >= 2000.0) & (t < 2005.0)
    tr = ~held
    V = {n: np.array([S[n][k] for k in common]) for n in S}
    ref = V["consensus"]
    qd = V["quad data"]
    mod = V["model"]

    # Descriptive missing-trend polynomials (fit to residual = ref - mod, and qd - mod)
    tc = t - t.mean()
    gap_cons = ref - mod
    gap_qd = qd - mod
    p1_cons = np.polyfit(tc, gap_cons, 1)
    p2_cons = np.polyfit(tc, gap_cons, 2)
    p1_qd = np.polyfit(tc, gap_qd, 1)
    mod_plus_lin = mod + np.polyval(p1_cons, tc)
    mod_plus_quad = mod + np.polyval(p2_cons, tc)
    V["model+lin_gap"] = mod_plus_lin
    V["model+quad_gap"] = mod_plus_quad

    # Area-weighted secular term from lt.exe.windings.json
    wsum = sum(weights[n] for n in quads)
    sec = sum(weights[n] * G.secular_curve(quads[n], t) for n in quads) / wsum
    sec = sec - np.mean(sec[(t >= 1961.0) & (t < 1991.0)])
    g_trend = sum(weights[n] * quads[n]["windings"]["trend"] for n in quads) / wsum
    g_accel = sum(weights[n] * quads[n]["windings"]["accel"] for n in quads) / wsum

    # Summary lines
    nq = len(quads)
    lines = [
        f"{nq} basin-seeded fast-refine quads, {t[0]:.0f}-{t[-1]:.2f} ({len(t)} months common to model, Kaplan, HadSST, NOAA)",
        f"ocean area covered by the {nq} quads: {100 * cov_frac:.1f}% of Kaplan global ocean; monthly weight present: {np.nanmean(frac):.3f}",
        "",
        "=== 1. MODEL & OBSERVATIONS VS 3-DATASET CONSENSUS (Kaplan, HadSST, NOAA) ===",
        f"{'series':16s} {'monthly CC / RMS':>18s} {'detr(lin) CC':>13s} {'detr(quad) CC':>14s} {'12-mo CC':>9s} {'detr 12-mo':>11s} {'held 2000-05':>13s}",
    ]
    for n in ("model", "model+lin_gap", "model+quad_gap", "quad data", "Kaplan", "HadSST", "NOAA"):
        a = V[n]
        c, r = score(a, ref, tr)
        cd1, _ = score(detr_lin(t, a), detr_lin(t, ref), tr)
        cd2, _ = score(detr_quad(t, a), detr_quad(t, ref), tr)
        cs, _ = score(run12(a), run12(ref), tr)
        cds, _ = score(run12(detr_lin(t, a)), run12(detr_lin(t, ref)), tr)
        ch, _ = score(a, ref, held)
        lines.append(
            f"{n:16s} {c:+.3f} / {r:.3f} C      {cd1:+.3f}         {cd2:+.3f}       {cs:+.3f}     {cds:+.3f}        {ch:+.3f}"
        )

    lines += [
        "",
        "=== 2. MODEL VS 30-QUAD DATA FIT DIRECTLY ===",
    ]
    c_mq, r_mq = score(mod, qd, tr)
    cd1_mq, _ = score(detr_lin(t, mod), detr_lin(t, qd), tr)
    cd2_mq, _ = score(detr_quad(t, mod), detr_quad(t, qd), tr)
    cs_mq, _ = score(run12(mod), run12(qd), tr)
    cds_mq, _ = score(run12(detr_lin(t, mod)), run12(detr_lin(t, qd)), tr)
    ch_mq, _ = score(mod, qd, held)
    lines.append(
        f"{'model vs qdata':16s} {c_mq:+.3f} / {r_mq:.3f} C      {cd1_mq:+.3f}         {cd2_mq:+.3f}       {cs_mq:+.3f}     {cds_mq:+.3f}        {ch_mq:+.3f}"
    )

    # Trend budget
    sl_mod = slope_c_per_century(t, mod)
    sl_qd = slope_c_per_century(t, qd)
    sl_kap = slope_c_per_century(t, V["Kaplan"])
    sl_had = slope_c_per_century(t, V["HadSST"])
    sl_noaa = slope_c_per_century(t, V["NOAA"])
    sl_cons = slope_c_per_century(t, ref)

    early = t < 1980.0
    late = t >= 1980.0
    c2_mod = float(np.polyfit(tc, mod, 2)[0])
    c2_qd = float(np.polyfit(tc, qd, 2)[0])
    c2_cons = float(np.polyfit(tc, ref, 2)[0])

    lines += [
        "",
        "=== 3. DESCRIPTIVE TREND & MISSING-TREND BUDGET (1950-2023) ===",
        f"  linear trend model     : {sl_mod:+.2f} C/century  (1950-80: {slope_c_per_century(t, mod, early):+.2f}, 1980-23: {slope_c_per_century(t, mod, late):+.2f}; quad c2={c2_mod:+.2e} C/yr^2)",
        f"  linear trend quad data : {sl_qd:+.2f} C/century  (1950-80: {slope_c_per_century(t, qd, early):+.2f}, 1980-23: {slope_c_per_century(t, qd, late):+.2f}; quad c2={c2_qd:+.2e} C/yr^2)",
        f"  linear trend Kaplan    : {sl_kap:+.2f} C/century  (1950-80: {slope_c_per_century(t, V['Kaplan'], early):+.2f}, 1980-23: {slope_c_per_century(t, V['Kaplan'], late):+.2f})",
        f"  linear trend HadSST    : {sl_had:+.2f} C/century",
        f"  linear trend NOAA      : {sl_noaa:+.2f} C/century",
        f"  linear trend consensus : {sl_cons:+.2f} C/century  (1950-80: {slope_c_per_century(t, ref, early):+.2f}, 1980-23: {slope_c_per_century(t, ref, late):+.2f}; quad c2={c2_cons:+.2e} C/yr^2)",
        "",
        f"  TOTAL missing linear trend (consensus - model): {sl_cons - sl_mod:+.2f} C/century ({(sl_cons - sl_mod) * (t[-1] - t[0]) / 100.0:+.3f} C over {t[-1] - t[0]:.0f} yr)",
        f"    (a) model vs 30-quad data gap (quad_data - model) : {sl_qd - sl_mod:+.2f} C/century",
        f"    (b) 30-quad spatial sampling  (Kaplan - quad_data): {sl_kap - sl_qd:+.2f} C/century",
        f"    (c) dataset trend offset      (consensus - Kaplan): {sl_cons - sl_kap:+.2f} C/century",
        f"  missing quadratic curvature (consensus - model)     : {c2_cons - c2_mod:+.2e} C/yr^2",
        f"  MLR secular term (windings.json): trend={g_trend * 100:+.2f} C/century, accel={g_accel:+.2e} C/yr^2, 1950-2023 rise={sec[-1] - sec[0]:+.3f} C",
        f"  observed spread (mean |dataset - consensus|)        : {np.mean([np.abs(V[o] - ref).mean() for o in obs]):.3f} C",
    ]

    print("\n".join(lines))
    with open(HERE / "global_sst_model.txt", "w") as fh:
        fh.write("# " + "\n# ".join(lines) + "\n")
        fh.write("# date model model+lin_gap consensus Kaplan HadSST NOAA quad_data\n")
        for i, k in enumerate(common):
            fh.write(
                f"{k / 12:.4f} "
                + " ".join(
                    f"{V[n][i]:+.4f}"
                    for n in ("model", "model+lin_gap", "consensus", "Kaplan", "HadSST", "NOAA", "quad data")
                )
                + "\n"
            )

    # 4-panel figure
    INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    fig, ax = plt.subplots(4, 1, figsize=(15, 13.5), facecolor=SURF, gridspec_kw=dict(height_ratios=[3.2, 2.2, 2.2, 2.0]))
    for A in ax:
        A.set_facecolor(SURF)
        A.axvspan(2000, 2005, color=GRID, alpha=0.85, lw=0)
        A.grid(color=GRID, lw=0.6)
        for s in ("top", "right"):
            A.spines[s].set_visible(False)

    # Panel 1: Full global SST
    A = ax[0]
    A.plot(t, ref, color="#b8b5ae", lw=0.7, label="consensus, monthly")
    A.plot(t, run12(ref), color=INK, lw=2.2, label="consensus (Kaplan, HadSST, NOAA), 12-mo")
    A.plot(t, run12(qd), color="#2a78d6", lw=1.6, ls="--", label=f"{nq}-quad data (area-weighted), 12-mo")
    A.plot(t, run12(mod), color="#eb6834", lw=2.1, label=f"{nq}-quad model (unadjusted), 12-mo")
    A.plot(t, run12(mod_plus_lin), color="#1baf7a", lw=1.5, ls=":", label=f"model + descriptive linear gap ({sl_cons - sl_mod:+.2f} C/cy)")
    A.set_title(
        f"Global Monthly SST: {nq} Basin-Seeded Fast-Refine Quads ({100 * cov_frac:.1f}% ocean area) vs Observed Consensus",
        loc="left", fontsize=12, color=INK,
    )
    A.set_ylabel("anomaly vs 1961-1990 (C)")
    A.legend(frameon=False, ncol=2, fontsize=9, loc="upper left")

    # Panel 2: Linearly detrended 12-mo variability
    A = ax[1]
    A.plot(t, run12(detr_lin(t, ref)), color=INK, lw=1.9, label="consensus, detrended")
    A.plot(t, run12(detr_lin(t, qd)), color="#2a78d6", lw=1.4, ls="--", label=f"{nq}-quad data, detrended")
    A.plot(t, run12(detr_lin(t, mod)), color="#eb6834", lw=1.9, label=f"{nq}-quad model, detrended")
    A.axhline(0, color=MUTED, lw=0.7)
    A.set_ylabel("detrended, 12-mo (C)")
    A.legend(frameon=False, ncol=3, fontsize=9, loc="upper left")

    # Panel 3: Residual & missing trend
    A = ax[2]
    for o, c in zip(obs, ("#2a78d6", "#1baf7a", "#e87ba4")):
        A.plot(t, run12(V[o] - ref), color=c, lw=1.1, alpha=0.85, label=f"{o} - consensus")
    A.plot(t, run12(ref - mod), color="#eb6834", lw=1.9, label="consensus - model (total gap)")
    A.plot(t, run12(qd - mod), color="#8b5cf6", lw=1.4, ls="--", label="quad_data - model")
    A.plot(t, np.polyval(p1_cons, tc), color=INK, lw=1.2, ls=":", label=f"linear gap fit ({sl_cons - sl_mod:+.2f} C/cy)")
    A.axhline(0, color=MUTED, lw=0.8)
    A.set_ylabel("residual, 12-mo (C)")
    A.legend(frameon=False, ncol=3, fontsize=8.5, loc="upper left")

    # Panel 4: Secular trend comparison
    A = ax[3]
    q_cons = np.polyval(np.polyfit(tc, ref, 2), tc)
    q_qd = np.polyval(np.polyfit(tc, qd, 2), tc)
    q_mod = np.polyval(np.polyfit(tc, mod, 2), tc)
    A.plot(t, q_cons, color=INK, lw=2.0, label=f"consensus quad trend ({sl_cons:+.2f} C/cy)")
    A.plot(t, q_qd, color="#2a78d6", lw=1.6, ls="--", label=f"{nq}-quad data quad trend ({sl_qd:+.2f} C/cy)")
    A.plot(t, q_mod, color="#eb6834", lw=2.0, label=f"{nq}-quad model quad trend ({sl_mod:+.2f} C/cy)")
    A.plot(t, sec, color="#1baf7a", lw=1.5, ls="-.", label=f"MLR secular term alone ({g_trend * 100:+.2f} C/cy)")
    A.axhline(0, color=MUTED, lw=0.7)
    A.set_ylabel("secular trend (C)")
    A.set_xlabel("Year (grey band = 2000-2005 held-out interval)")
    A.legend(frameon=False, ncol=2, fontsize=9, loc="upper left")

    for A in ax:
        A.set_xlim(t[0], t[-1])
    fig.tight_layout()
    fig.savefig(HERE / "global_sst_model.png", dpi=120, facecolor=SURF)
    print(f"plot: {HERE / 'global_sst_model.png'}")


if __name__ == "__main__":
    main()
