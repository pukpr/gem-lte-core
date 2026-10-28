#!/usr/bin/env python3
"""The CW's 6.4-yr / 6.0-yr (F - M) beat, carried over to the QBO (2026-10-26). Discussion: qbo_six_year_beat.md.

Arithmetic. The CW alias is twice the QBO alias: 2F - 26 = 2 (F - 13), so F_CW = 2q with q = 0.4219 /yr. In polar
motion (ChandlerWobble/docs/cw_suppression.md section 6.1) two 6-year signals meet:
  annual - CW = 1 - 2q = 0.1560 /yr (6.41 yr)   and   F - M = 0.1668 /yr (6.00 yr, draconic - anomalistic),
and folded by the annual (odd-harmonic) impulse, F - M sits at 1 - (F - M) = 0.8332 /yr beside the CW at 0.8438: the two
lines beat with a ~93-yr envelope (CW minima fitted at 1924-29 and 2017-22; out of sample, 1890-1972, r 0.42-0.49).
In the QBO (an odd-harmonic, sign-sensitive comb, the drive in which F - M is native) the same pair appears twice:
  band 0.84:  2q = 0.8438 (the QBO's second harmonic, = F_CW)  vs  1 - (F - M) = 0.8332
  band 0.58:  1 - q = 0.5781 (annual sideband)                vs  q + (F - M)   = 0.5887
  (and q - (F - M) = M - 13 = 0.2551, the anomalistic alias)
both separated by 2q + (F - M) - 1 = 0.0107 /yr (93 yr). The QBO record (KIT radiosondes, 1953-2026, 74 yr) cannot
resolve either pair, but a beat shows as a slow decline of the band's envelope toward a node, with a phase drift.

Method (as ChandlerWobble/docs/cw_longrecord.py): sliding 12-yr windows, 1-yr step; in each window a joint
least-squares fit of fixed lines q, 1-q, 1+q, 2-q, 2+q, 2q, annual, semiannual and a linear trend; amplitude and phase
of the 1-q and 2q bands per window (phase relative to the 1960-1985 circular mean). CW comparison: the C01 track of
cw_longrecord.json (same windows, line at f_CW).
    python3 qbo_six_year_beat.py   -> qbo_six_year_beat.png, qbo_six_year_beat.json
"""
import json
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
KIT = HERE.parent / "experiments/Oct2026_generator/QBO/kit_qbo.dat"
CWL = Path.home() / "github/pukpr/ChandlerWobble/docs/cw_longrecord.json"
CWS = Path.home() / "github/pukpr/ChandlerWobble/docs/cw_suppression.json"
LEVELS = [70, 50, 40, 30, 20, 15, 10]; ENDS = [16, 23, 30, 37, 44, 51, 58]
Y = 365.2422
q = Y / 27.212221 % 1; fFM = Y / 27.212221 - Y / 27.554550
LINES = {"q": q, "1-q": 1 - q, "1+q": 1 + q, "2-q": 2 - q, "2+q": 2 + q, "2q": 2 * q, "annual": 1.0, "semi": 2.0}
WIN = 12.0


def kit():
    rows = {}
    for ln in KIT.read_text().splitlines()[9:]:
        if len(ln) < 20 or not ln[:5].strip().isdigit():
            continue
        yy, mo = int(ln[6:8]), int(ln[8:10]); yr = (1900 if yy >= 50 else 2000) + yy
        v = {m.end() + 10: int(m.group()) for m in re.finditer(r"-?\d+", ln[10:])}; rows[(yr, mo)] = [v.get(e, np.nan) / 10 for e in ENDS]
    k = sorted(rows); return np.array([y + (m - 0.5) / 12 for y, m in k]), np.array([rows[x] for x in k])


def track(t, u):
    out = []
    for c in np.arange(np.ceil(t[0] + WIN / 2), np.floor(t[-1] - WIN / 2) + 1):
        m = (t >= c - WIN / 2) & (t < c + WIN / 2) & np.isfinite(u)
        if m.sum() < 0.9 * WIN * 12:
            continue
        tt = t[m] - c; cols = [np.ones(m.sum()), tt]
        for f in LINES.values():
            cols += [np.cos(2 * np.pi * f * t[m]), np.sin(2 * np.pi * f * t[m])]
        cf = np.linalg.lstsq(np.column_stack(cols), u[m], rcond=None)[0][2:]
        row = {"year": float(c)}
        for i, k in enumerate(LINES):
            a, b = cf[2 * i], cf[2 * i + 1]; row[k] = complex(a, -b)                 # u ~ Re[A e^{i 2 pi f t}]
        out.append(row)
    return out


def main():
    t, U = kit()
    res = {"q": q, "f_FM": fFM, "pairs": {"0.84": [2 * q, 1 - fFM], "0.58": [1 - q, q + fFM]}, "beat_period_yr": 1 / (2 * q + fFM - 1),
           "levels": {}}
    trk = {}
    for L in (20, 30, 50, 70):
        tr = track(t, U[:, LEVELS.index(L)]); trk[L] = tr
        yrs = np.array([r["year"] for r in tr]); ref = (yrs >= 1960) & (yrs <= 1985)
        lv = {}
        for band in ("q", "1-q", "2q"):
            z = np.array([r[band] for r in tr]); ph0 = np.angle(np.mean(z[ref] / np.abs(z[ref])))
            amp = np.abs(z); ph = np.degrees(np.angle(z * np.exp(-1j * ph0)))
            k = np.argmin(amp); slope = np.polyfit(yrs[yrs <= 2010], amp[yrs <= 2010], 1)[0]
            lv[band] = {"years": yrs.tolist(), "amp": amp.tolist(), "phase_deg": ph.tolist(), "min_year": float(yrs[k]),
                        "min_amp": float(amp[k]), "max_amp": float(amp.max()), "trend_to_2010_per_decade": float(10 * slope),
                        "amp_first_decade": float(amp[:10].mean()), "amp_last_decade": float(amp[-10:].mean())}
        # self-referenced: harmonic / sideband relative to the QBO line itself (removes the QBO's own phase wandering)
        zq = np.array([r["q"] for r in tr]); z2 = np.array([r["2q"] for r in tr]); z1 = np.array([r["1-q"] for r in tr])
        for name, rel in (("2q_rel", z2 * np.conj(zq) ** 2 / np.abs(zq) ** 2), ("1-q_rel", z1 * zq / np.abs(zq))):   # |rel| = |band|
            ph0 = np.angle(np.mean(rel[ref] / np.abs(rel[ref]))); amp = np.abs(rel) / np.abs(zq)          # |band| / |q| in both cases
            ph = np.degrees(np.unwrap(np.angle(rel * np.exp(-1j * ph0))))
            lv[name] = {"years": yrs.tolist(), "ratio": amp.tolist(), "phase_deg": ph.tolist(),
                        "phase_drift_deg_per_decade": float(10 * np.polyfit(yrs, ph, 1)[0]),
                        "beat_drift_expected_deg_per_decade": float((-1 if name == "2q_rel" else 1) * 10 * 360 / res["beat_period_yr"])}
        res["levels"][L] = lv
    cw = json.loads(CWL.read_text()); ctr = [r for r in cw["track"] if r["year"] >= 1950]
    cy = np.array([r["year"] for r in ctr]); ca = np.array([r["amp_mas"] for r in ctr]); cp = np.array([r["phase_deg"] for r in ctr])
    res["cw"] = {"years": cy.tolist(), "amp": ca.tolist(), "phase": cp.tolist(), "min_year_after_1950": float(cy[np.argmin(ca)]),
                 "fm_beat_minima": {k: v["F-M (0.8332)"]["beat_minima"] for k, v in json.loads(CWS.read_text())["two_line"].items()},
                 "fm_out_of_sample": cw["fm_beat_out_of_sample"]}
    # descriptive comparison: QBO band envelopes vs CW envelope on common window centres (few independent samples)
    res["corr_with_cw_amp"] = {}
    for L in (30, 50):
        for band in ("1-q", "2q"):
            qy = np.array(res["levels"][L][band]["years"]); qa = np.array(res["levels"][L][band]["amp"])
            com = np.intersect1d(qy, np.round(cy)); a1 = np.interp(com, qy, qa); a2 = np.interp(com, cy, ca)
            res["corr_with_cw_amp"][f"{L}hPa_{band}"] = {"r": float(np.corrcoef(a1, a2)[0, 1]), "n_windows": int(len(com)),
                                                         "independent_approx": int(len(com) / WIN) + 1}
    (HERE / "qbo_six_year_beat.json").write_text(json.dumps(res, indent=1, default=float))
    report(res); plot(res)


def report(res):
    print(f"q {res['q']:.4f}, F-M {res['f_FM']:.4f} /yr ({1 / res['f_FM']:.3f} yr); pairs 0.84: {res['pairs']['0.84']}, 0.58: {res['pairs']['0.58']}; beat {res['beat_period_yr']:.1f} yr")
    for L, lv in res["levels"].items():
        for band, v in lv.items():
            if band.endswith("_rel"):
                continue
            print(f"  {L} hPa {band:4s}: amp first decade {v['amp_first_decade']:.2f}, last decade {v['amp_last_decade']:.2f} m/s; "
                  f"trend to 2010 {v['trend_to_2010_per_decade']:+.2f} m/s/decade; min {v['min_amp']:.2f} at {v['min_year']:.0f} "
                  f"(max {v['max_amp']:.2f})")
    c = res["cw"]; print(f"CW (C01, 12-yr windows) min after 1950 at {c['min_year_after_1950']:.0f}; F-M two-line minima {c['fm_beat_minima']}")
    print("descriptive r(QBO band envelope, CW envelope):", {k: round(v["r"], 2) for k, v in res["corr_with_cw_amp"].items()},
          "(~5 independent windows)")
    for L in (20, 30, 50, 70):
        for nm in ("2q_rel", "1-q_rel"):
            v = res["levels"][L][nm]; y = np.array(v["years"])
            print(f"  {L} hPa {nm}: ratio " + ", ".join(f"{int(a)}:{b:.2f}" for a, b in zip(y[::6], np.array(v["ratio"])[::6]))
                  + f" | phase " + ", ".join(f"{int(a)}:{b:+.0f}" for a, b in zip(y[::6], np.array(v["phase_deg"])[::6]))
                  + f" | drift {v['phase_drift_deg_per_decade']:+.0f} deg/decade (F - M beat: {v['beat_drift_expected_deg_per_decade']:+.0f})")
    for L in (30, 50):
        v = res["levels"][L]["2q"]; y = np.array(v["years"]); p = np.array(v["phase_deg"])
        print(f"  {L} hPa 2q phase by window centre: " + ", ".join(f"{int(a)}:{b:+.0f}" for a, b in zip(y[::4], p[::4])))
        v = res["levels"][L]["1-q"]; p = np.array(v["phase_deg"])
        print(f"  {L} hPa 1-q phase by window centre: " + ", ".join(f"{int(a)}:{b:+.0f}" for a, b in zip(y[::4], p[::4])))


def plot(res):
    fig, ax = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
    c = res["cw"]; a = ax[0]; a.plot(c["years"], c["amp"], "k-", lw=1.5, label="CW amplitude (C01, 12-yr windows)")
    b = a.twinx(); b.plot(c["years"], c["phase"], color="#d94a4a", lw=1, label="CW phase"); b.set_ylabel("CW phase (deg)", color="#d94a4a")
    for k, mins in c["fm_beat_minima"].items():
        for m_ in mins:
            if 1950 < m_ < 2030: a.axvspan(m_ - 0.5, m_ + 0.5, color="#eda100", alpha=0.25)
    a.set_ylabel("CW amplitude (mas)"); a.legend(loc="lower left", fontsize=8)
    a.set_title("(a) Chandler wobble envelope; orange: minima of the CW + (F − M) two-line beat (fitted 1973–2026, several Q)", fontsize=10)
    for ax_, band, ttl, pred in ((ax[1], "2q_rel", "(b) QBO band 0.84 /yr relative to the QBO line: |2q|/|q| and phase(2q) − 2 phase(q)", -1),
                                 (ax[2], "1-q_rel", "(c) QBO band 0.58 /yr relative to the QBO line: |1−q|/|q| and phase(1−q) + phase(q)", +1)):
        b = ax_.twinx()
        for L, col in ((30, "k"), (50, "#2a78d6"), (20, "#1baf7a"), (70, "#eb6834")):
            v = res["levels"][L][band]; ax_.plot(v["years"], v["ratio"], color=col, lw=1.4, label=f"{L} hPa ratio")
            b.plot(v["years"], v["phase_deg"], color=col, lw=0.8, ls="--")
        y0 = np.array(res["levels"][30][band]["years"]); b.plot(y0, pred * 360 / res["beat_period_yr"] * (y0 - 1972), color="#eda100", lw=2.5, alpha=0.6,
                                                                label=f"F − M beat drift ({pred * 39:+.0f}°/decade)")
        for x in (2016.2, 2020.0): ax_.axvline(x, color="#d94a4a", ls=":", lw=1)
        ax_.set_ylabel("band / QBO line"); ax_.set_ylim(0, 1.5); b.set_ylabel("relative phase, unwrapped (deg, dashed)")
        ax_.legend(loc="upper left", fontsize=8, ncol=4); b.legend(loc="lower left", fontsize=8)
        ax_.set_title(ttl + "; dotted: disruptions", fontsize=10)
    ax[2].set_xlabel("window centre (12-yr windows)"); ax[2].set_xlim(1950, 2026)
    fig.tight_layout(); fig.savefig(HERE / "qbo_six_year_beat.png", dpi=105); plt.close(fig)


if __name__ == "__main__":
    main()
