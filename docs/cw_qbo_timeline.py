#!/usr/bin/env python3
"""When did the Chandler wobble collapse and the QBO go off schedule, and how does each show in its model? (2026-10-26)
Discussion: cw_qbo_timeline.md.

Both devices share the semiannual comb and the draconic month (ChandlerWobble/docs/cw_suppression.md section 5):
  CW   sign-blind drive B^2 on a same-polarity comb: alias 2 x 13.4219 - 26 = 0.8438 /yr (comb harmonic 26)
  QBO  sign-sensitive drive B on an alternating comb: alias 13.4219 - 13 = 0.4219 /yr (comb harmonic 13)
A shift dt of the comb timing rotates the CW alias forcing by 360*26*dt and the QBO's by 360*13*dt, so the QBO rotation
should be half the CW one (determined only modulo 180 deg, since the CW rotation is known modulo 360).

Inputs (nothing refitted here):
  CW   ChandlerWobble/docs/cw_suppression.json: envelope and phase (raw IERS polar motion demodulated at 0.8440 /yr,
       Gaussian sigma 1 yr) and the required forcing per epoch for several (f0, Q)
  QBO  KIT radiosonde winds (kit_qbo.dat, 70-10 hPa): demodulated at 0.4219 /yr with the same smoothing; and the
       required pulse train of the inverted QBO device (qbo_control.json, from qbo_control.py), whose alias component
       per epoch, divided by the astronomical comb's, is the QBO's "required forcing", on the same epochs as the CW.
    python3 cw_qbo_timeline.py   -> cw_qbo_timeline.png, cw_qbo_timeline.json
"""
import json
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
CWJ = Path.home() / "github/pukpr/ChandlerWobble/docs/cw_suppression.json"
KIT = HERE.parent / "experiments/Oct2026_generator/QBO/kit_qbo.dat"
Q = 365.2422 / 27.212221 % 1
EPOCHS = [1973, 1990, 2005, 2010, 2013, 2016, 2019, 2027]
LEVELS = [70, 50, 40, 30, 20, 15, 10]; ENDS = [16, 23, 30, 37, 44, 51, 58]


def kit():
    rows = {}
    for ln in KIT.read_text().splitlines()[9:]:
        if len(ln) < 20 or not ln[:5].strip().isdigit():
            continue
        yy, mo = int(ln[6:8]), int(ln[8:10]); yr = (1900 if yy >= 50 else 2000) + yy
        v = {m.end() + 10: int(m.group()) for m in re.finditer(r"-?\d+", ln[10:])}; rows[(yr, mo)] = [v.get(e, np.nan) / 10 for e in ENDS]
    k = sorted(rows); return np.array([y + (m - 0.5) / 12 for y, m in k]), np.array([rows[x] for x in k])


def demod(t, u, f, sigma=1.0):
    z = np.where(np.isfinite(u), u, 0.0) * np.exp(-2j * np.pi * f * t); w = np.isfinite(u).astype(float)
    dt = np.median(np.diff(t)); k = np.arange(-int(4 * sigma / dt), int(4 * sigma / dt) + 1) * dt
    g = np.exp(-0.5 * (k / sigma) ** 2)
    return 2 * np.convolve(z, g, "same") / np.convolve(w, g, "same")


def main():
    cw = json.loads(CWJ.read_text())
    env = cw["envelope"]; yrs = np.array([e["year"] for e in env]); camp = np.array([e["amp_mas"] for e in env])
    cph = np.array([e["phase_deg"] for e in env])
    t, U = kit()
    qd = {L: demod(t, U[:, LEVELS.index(L)], Q) for L in (70, 50, 40, 30, 20)}
    ref = {L: np.angle(np.mean(qd[L][(t >= 1973) & (t < 1990)])) for L in qd}
    ann = lambda L, y: qd[L][np.argmin(np.abs(t - (y + 0.5)))]
    out = {"qbo_alias_cpy": Q, "cw_alias_cpy": 2 * 365.2422 / 27.212221 - 26, "years": {}}
    for y in range(1975, 2026):
        out["years"][y] = {"cw_amp_mas": float(camp[yrs == y][0]) if y in yrs else None,
                           "cw_phase_deg": float(cph[yrs == y][0]) if y in yrs else None,
                           **{f"qbo{L}_amp": float(abs(ann(L, y))) for L in qd},
                           **{f"qbo{L}_phase": float(np.degrees(np.angle(ann(L, y) * np.exp(-1j * ref[L])))) for L in qd}}
    # required forcing per epoch, QBO device
    qc = json.loads((HERE / "qbo_control.json").read_text())["pulses"]
    tn, ar, aa = (np.array(qc[k]) for k in ("t", "a_req", "a_ast"))

    def comp(a, m): return np.sum(a[m] * np.exp(-2j * np.pi * Q * tn[m]))
    base = None; qf = []
    for e0, e1 in zip(EPOCHS[:-1], EPOCHS[1:]):
        m = (tn >= e0) & (tn < e1); R = comp(ar, m) / comp(aa, m); base = R if base is None else base
        rr = R / base; qf.append({"epoch": [e0, e1], "rel_amp": float(abs(rr)), "rel_phase": float(np.degrees(np.angle(rr))), "pulses": int(m.sum())})
    out["qbo_required_forcing"] = qf
    out["cw_required_forcing"] = {k: v["forcing"] for k, v in cw["cases"].items()}
    # comb-shift prediction: QBO rotation = CW rotation / 2 (mod 180)
    out["comb_shift_check"] = []
    for case, fl in out["cw_required_forcing"].items():
        rows = []
        for f, g in zip(fl, qf):
            pred = f["rel_phase"] / 2
            alt = pred + 180 if pred < 0 else pred - 180
            rows.append({"epoch": f["epoch"], "cw": f["rel_phase"], "qbo_pred": [pred, alt], "qbo_obs": g["rel_phase"],
                         "miss": float(min(abs((g["rel_phase"] - p + 180) % 360 - 180) for p in (pred, alt)))})
        out["comb_shift_check"].append({"case": case, "rows": rows})
    # disruption months: easterly (u < -5 m/s) at 40 hPa while 30 and 50 hPa are westerly, 2014-2021
    dis = [(float(tt), float(U[i, 2]), float(U[i, 3]), float(U[i, 1])) for i, tt in enumerate(t)
           if 2014 <= tt < 2021.5 and U[i, 2] < -5 and U[i, 3] > 0 and U[i, 1] > 0]
    out["disruption_months_40hPa"] = dis
    (HERE / "cw_qbo_timeline.json").write_text(json.dumps(out, indent=1, default=float))
    report(out, yrs, camp, cph, t, qd, ref)
    plot(out, yrs, camp, cph, t, U, qd, ref)


def report(out, yrs, camp, cph, t, qd, ref):
    print("year   CW amp(mas) phase | QBO30 amp phase | QBO50 amp phase | QBO70 amp phase")
    for y in range(2000, 2026):
        r = out["years"][y]
        cw_ = f"{r['cw_amp_mas']:6.0f} {r['cw_phase_deg']:+5.0f}" if r["cw_amp_mas"] is not None else "     -     -"
        print(f"{y}  {cw_}     | {r['qbo30_amp']:5.1f} {r['qbo30_phase']:+5.0f}   | {r['qbo50_amp']:5.1f} {r['qbo50_phase']:+5.0f}   | {r['qbo70_amp']:5.1f} {r['qbo70_phase']:+5.0f}")
    print("\nrequired forcing per epoch (rel. 1973-90): QBO device  | CW Q85 | CW Q32")
    c85 = out["cw_required_forcing"]["consensus 432 d, Q 85"]; c32 = out["cw_required_forcing"]["fitted f0 0.843, Q 32"]
    for g, a, b in zip(out["qbo_required_forcing"], c85, c32):
        print(f"  {g['epoch']}: QBO {g['rel_amp']:.2f} / {g['rel_phase']:+5.0f} ({g['pulses']} pulses) | CW {a['rel_amp']:.2f} / {a['rel_phase']:+5.0f} | {b['rel_amp']:.2f} / {b['rel_phase']:+5.0f}")
    print("\ncomb-shift check (QBO rotation should be CW rotation / 2, mod 180):")
    for c in out["comb_shift_check"]:
        print(f"  {c['case']}: misses " + ", ".join(f"{r['epoch'][0]}-{r['epoch'][1]} {r['miss']:.0f}" for r in c["rows"][1:]))
    print("\n40 hPa easterly while 30 and 50 hPa westerly (disruption signature), 2014-2021:")
    for tt, u40, u30, u50 in out["disruption_months_40hPa"]:
        print(f"  {tt:.2f}: u40 {u40:+.1f}, u30 {u30:+.1f}, u50 {u50:+.1f} m/s")


def plot(out, yrs, camp, cph, t, U, qd, ref):
    fig, ax = plt.subplots(4, 1, figsize=(14, 15), gridspec_kw={"height_ratios": [1, 1, 1, 1.1]})
    a = ax[0]; a.plot(yrs, camp, "k-o", ms=3, label="CW envelope (mas)"); a.set_ylabel("CW amplitude (mas)")
    b = a.twinx(); b.plot(yrs, cph, color="#d94a4a", lw=1, label="CW phase (deg)"); b.set_ylabel("CW phase (deg)", color="#d94a4a")
    a.set_title("(a) Chandler wobble at its alias 0.8438 /yr (from cw_suppression.json): collapse 2010–15, phase swing through the minimum", fontsize=10)
    for x in (2010, 2016, 2019): a.axvline(x, color="0.5", ls=":", lw=0.8)
    a.set_xlim(1975, 2026); a.legend(loc="upper right", fontsize=8)
    m = (t >= 1975) & (t <= 2026)
    a = ax[1]; b = a.twinx()
    for L, col in ((30, "k"), (50, "#2a78d6"), (70, "#1baf7a")):
        a.plot(t[m], np.abs(qd[L][m]), color=col, lw=1.2, label=f"{L} hPa amplitude")
        b.plot(t[m], np.degrees(np.angle(qd[L][m] * np.exp(-1j * ref[L]))), color=col, lw=0.8, ls="--")
    a.set_ylabel("QBO amplitude at 0.4219 /yr (m/s)"); b.set_ylabel("phase rel. 1973–90 (deg, dashed)")
    for x in (2010, 2016, 2019): a.axvline(x, color="0.5", ls=":", lw=0.8)
    a.axvspan(2015.8, 2016.6, color="#f6d5d5"); a.axvspan(2019.8, 2020.4, color="#f6d5d5")
    a.set_xlim(1975, 2026); a.legend(loc="lower left", fontsize=8)
    a.set_title("(b) QBO (KIT radiosondes) demodulated at its alias 0.4219 /yr, same smoothing; pink: the 2016 and 2019–20 disruptions", fontsize=10)
    a = ax[2]
    ep = [0.5 * (g["epoch"][0] + g["epoch"][1]) for g in out["qbo_required_forcing"]]
    c85 = out["cw_required_forcing"]["consensus 432 d, Q 85"]; c32 = out["cw_required_forcing"]["fitted f0 0.843, Q 32"]
    a.plot(ep, [g["rel_phase"] for g in out["qbo_required_forcing"]], "o-", color="k", label="QBO device: required forcing phase")
    a.plot(ep, [f["rel_phase"] / 2 for f in c85], "s--", color="#d94a4a", label="CW Q 85 phase / 2 (comb-shift prediction for the QBO)")
    a.plot(ep, [f["rel_phase"] / 2 for f in c32], "^--", color="#eb6834", label="CW Q 32 phase / 2")
    for e in EPOCHS: a.axvline(e, color="0.85", lw=0.6)
    a.axhline(0, color="0.6", lw=0.6); a.set_ylabel("phase rel. 1973–90 (deg)"); a.set_xlim(1975, 2026); a.legend(fontsize=8)
    a.set_title("(c) Required forcing per epoch: the QBO device's vs half the CW's (a comb-timing shift predicts QBO = CW/2, mod 180)", fontsize=10)
    a = ax[3]; m = (t >= 2013.5) & (t <= 2021.5)
    for i, L in enumerate((20, 30, 40, 50, 70)):
        u = U[m, LEVELS.index(L)]; a.plot(t[m], u / 15 - 2.2 * i, color="k", lw=1)
        a.fill_between(t[m], -2.2 * i, u / 15 - 2.2 * i, where=u < 0, color="#2a78d6", alpha=0.35, lw=0)
        a.fill_between(t[m], -2.2 * i, u / 15 - 2.2 * i, where=u > 0, color="#d94a4a", alpha=0.3, lw=0)
        a.text(2013.55, -2.2 * i + 0.4, f"{L} hPa", fontsize=8)
    for tt, *_ in out["disruption_months_40hPa"]:
        a.axvline(tt, color="#eda100", lw=2, alpha=0.5)
    a.set_yticks([]); a.set_xlim(2013.5, 2021.5); a.set_xlabel("year")
    a.set_title("(d) KIT winds 2013.5–2021.5 (red westerly, blue easterly); orange: 40 hPa easterly inside a 30/50 hPa westerly (disruption)", fontsize=10)
    fig.tight_layout(); fig.savefig(HERE / "cw_qbo_timeline.png", dpi=105); plt.close(fig)


if __name__ == "__main__":
    main()
