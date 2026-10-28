#!/usr/bin/env python3
"""Test (c): does the wave activity entering the QBO column carry the tidal comb lines? (2026-10-27)
Gate hypothesis: the tide sets WHEN the upward wave flux gets through (F_in = F0 [1 + eps g(M(t))]), so the envelope
of the 4-20 d waves near the tropopause should carry the comb lines: q = 0.4219 cpy (13.4219 - 13, draconic alias),
2q, and the sidebands 1 +/- q, 2 +/- q around the annual and semiannual lines.
Data: daily ERA5 u, 0N 0E 00 UTC, at 100 (when present), 70, 50, 30, 20 hPa (get_era5_point.py). Wave activity E(t) =
4-20 d band-passed u squared, averaged to months (log E is used: a gate is multiplicative).
Joint least squares for all lines + linear trend; noise level from the same fit at 400 random off-line frequencies in
0.15-3.2 cpy (at least 1.5/T from every line); SNR = amplitude / median off-line amplitude; p = fraction of off-line
amplitudes >= the line's.
CONFOUND: the QBO wind itself modulates wave activity (Doppler shift, critical-level filtering), and the observed QBO
period is close to 1/q by construction of the model, so a QBO-modulated envelope has q lines whatever forces the QBO;
seasonal (annual/semiannual) convection times a QBO modulation makes the 1 +/- q, 2 +/- q sidebands too. Control: the
same fit after removing a multiplicative season x local-state model, log E = s(month) + b(u bin, shear sign) +
b70(u70 bin) (u, shear at the level; 8 u bins). Lines that survive the control are not explained by the local or
70 hPa QBO state. The clean place for the test is 100 hPa, below most of the QBO (ERA5 100 hPa is downloading).
    python3 envelope_lines.py   -> envelope_lines.png, envelope_lines.json"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import butter, filtfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
Q = 13.4219 - 13
LINES = {"q": Q, "2q": 2 * Q, "1-q": 1 - Q, "1": 1.0, "1+q": 1 + Q, "2-q": 2 - Q, "2": 2.0, "2+q": 2 + Q}
LV = [L for L in (100, 70, 50, 30, 20) if (HERE / f"era_qbo_u{L}_point0E_daily.dat").exists()]
rng = np.random.default_rng(1)


def monthly(td, x):
    ym = np.floor(td * 12 + 1e-9).astype(int); k = np.unique(ym)
    return (k + 0.5) / 12.0, np.array([x[ym == j].mean() for j in k])


def design(t, freqs):
    cols = [np.ones_like(t), t - t.mean()]
    for f in freqs: cols += [np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t)]
    return np.array(cols).T


def line_fit(t, y):
    f = list(LINES.values()); c = np.linalg.lstsq(design(t, f), y, rcond=None)[0][2:].reshape(-1, 2); amp = np.hypot(*c.T)
    T = t[-1] - t[0]; off = []
    while len(off) < 400:
        g = rng.uniform(0.15, 3.2)
        if min(abs(g - v) for v in f) > 1.5 / T:
            cg = np.linalg.lstsq(design(t, f + [g]), y, rcond=None)[0][-2:]; off.append(np.hypot(*cg))
    off = np.array(off)
    return {n: {"amp": float(a), "snr": float(a / np.median(off)), "p": float(np.mean(off >= a))} for n, a in zip(LINES, amp)}


def spectrum(t, y, fr):
    y = y - np.polyval(np.polyfit(t, y, 1), t); w = np.hanning(len(t))
    return np.array([abs(np.sum(w * y * np.exp(-2j * np.pi * f * t))) * 2 / w.sum() for f in fr])


out, spec = {}, {}
b, a = butter(3, [1 / 20, 1 / 4], btype="band", fs=1.0)
U = {L: np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat") for L in LV}
td = U[LV[0]][:, 0]
fr = np.linspace(0.1, 3.3, 1200)
for L in LV:
    d = U[L]; tL, uL = d[:, 0], d[:, 1]
    tm, E = monthly(tL, filtfilt(b, a, uL - uL.mean()) ** 2); logE = np.log(E)
    _, um = monthly(tL, uL)
    # local shear from the neighbouring ERA levels (log-p), monthly
    i = LV.index(L); hi = LV[i + 1] if i + 1 < len(LV) else L; lo = LV[i - 1] if i > 0 else L   # LV runs 100 -> 20 hPa
    mon = lambda x: np.interp(tm, *monthly(U[x][:, 0], U[x][:, 1]))
    shear = np.sign(mon(hi) - mon(lo))                                  # + = westerly shear (u increases upward)
    u70 = np.interp(tm, *monthly(U[70][:, 0], U[70][:, 1]))
    bins = lambda x: np.digitize(x, np.quantile(x, np.linspace(0, 1, 9)[1:-1]))
    month = (np.floor(tm * 12) % 12).astype(int)
    X = np.column_stack([np.eye(12)[month], np.eye(8)[bins(um)] * (shear[:, None] > 0), np.eye(8)[bins(um)] * (shear[:, None] <= 0), np.eye(8)[bins(u70)][:, 1:]])
    resid = logE - X @ np.linalg.lstsq(X, logE, rcond=None)[0]
    out[L] = {"raw": line_fit(tm, logE), "after_state_control": line_fit(tm, resid), "wind_itself": line_fit(tm, um),
              "state_model_r2": float(1 - resid.var() / logE.var()), "months": int(len(tm)), "span": [float(tm[0]), float(tm[-1])]}
    spec[L] = (spectrum(tm, logE, fr), spectrum(tm, resid, fr))
json.dump({"q": Q, "lines": LINES, "levels": out}, open(HERE / "envelope_lines.json", "w"), indent=1)

fig, ax = plt.subplots(len(LV), 1, figsize=(14, 3.0 * len(LV)), facecolor=SURF, sharex=True)
ax = np.atleast_1d(ax)
for i, L in enumerate(LV):
    s0, s1 = spec[L]; a_ = ax[i]
    a_.plot(fr, s0, color=MUTED, lw=1, label="log wave activity (4–20 d), raw")
    a_.plot(fr, s1, color="#d94a4a", lw=1.3, label="after removing season × QBO-state model")
    top = max(s0.max(), s1.max())
    for n, f in LINES.items():
        a_.axvline(f, color="#2a78d6", lw=0.8, ls=":"); a_.text(f, top * 1.02, n, ha="center", fontsize=8, color="#2a78d6")
    o = out[L]
    txt = "  ".join(f"{n}: {o['raw'][n]['snr']:.1f}→{o['after_state_control'][n]['snr']:.1f}" for n in LINES)
    a_.set_title(f"{L} hPa — line SNR raw→controlled: {txt}  (state model R² {o['state_model_r2']:.2f})", loc="left", fontsize=8.5, color=INK)
    a_.set_ylabel("amplitude", color=MUTED); a_.set_ylim(0, top * 1.12); a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
    if i == 0: a_.legend(fontsize=8.5, frameon=False, loc="upper right")
ax[-1].set_xlabel("frequency (cycles per year)", color=MUTED)
fig.tight_layout(); fig.savefig(HERE / "envelope_lines.png", dpi=105, facecolor=SURF)

print(f"q = {Q:.4f} cpy; SNR = amplitude / median off-line amplitude (p = off-line fraction >= line)")
for L, o in out.items():
    print(f"== {L} hPa ({o['months']} months, {o['span'][0]:.1f}-{o['span'][1]:.1f}); season x state model R2 {o['state_model_r2']:.2f}")
    for kind in ("wind_itself", "raw", "after_state_control"):
        print(f"  {kind:20s} " + " ".join(f"{n}:{o[kind][n]['snr']:4.1f}(p{o[kind][n]['p']:.2f})" for n in LINES))
print("wrote envelope_lines.png")
