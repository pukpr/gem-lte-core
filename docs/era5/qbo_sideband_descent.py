"""Why the QBO 2 +/- q sideband mix changes between 30 and 70 hPa: calendar-locked timing gate + descending response.
(2026-10-26; needs qbo_sideband_rules.json from qbo_sideband_rules.py.)
The sidebands at 2 +/- q come from the semiannual parts of the QBO's seasonal amplitude modulation (AM, alpha_2) and
timing modulation (PM, phi_2). If the timing gate is fixed in calendar time for the whole column and the AM descends
with the QBO (delay tau below 10 hPa), their relative phase rotates by -720 deg x tau(yr): over the ~8 months between
30 and 70 hPa, more than a full swap, while the 2.37-yr carrier itself shifts by only ~0.28 cycle.
Tests: (1) descent delay tau(z) from lag correlation with 10 hPa (independent of the seasonal analysis) -> predicted
winding -720 tau vs the measured AM-PM winding; (2) forward model: one common timing gate (mean phi_2 over levels) and
the 10 hPa AM delayed by tau (B) or with measured depth per level (A) -> predicted lower/upper sideband ratio.
usage: python3 qbo_sideband_descent.py [path to kit_qbo.dat]"""
import json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "../../experiments/Oct2026_generator/QBO/kit_qbo.dat")
LEVELS = [70, 50, 40, 30, 20, 15, 10]; ENDS = [16, 23, 30, 37, 44, 51, 58]
rows = {}
for ln in open(SRC).read().splitlines()[9:]:
    if len(ln) < 20 or not ln[:5].strip().isdigit():
        continue
    yy, mm = int(ln[6:8]), int(ln[8:10]); yr = (1900 if yy >= 50 else 2000) + yy
    v = {m.end() + 10: int(m.group()) for m in re.finditer(r"-?\d+", ln[10:])}; rows[(yr, mm)] = [v.get(e, np.nan) / 10 for e in ENDS]
U = np.array([rows[k] for k in sorted(rows)]); ok = np.all(np.isfinite(U), axis=1)
R = json.load(open(os.path.join(HERE, "qbo_sideband_rules.json")))["levels"]
LV = LEVELS[::-1]

tau = {}
ref = U[ok, LEVELS.index(10)] - U[ok, LEVELS.index(10)].mean()
for j, L in enumerate(LEVELS):
    x = U[ok, j] - U[ok, j].mean(); cc = [np.corrcoef(ref[:len(ref) - l], x[l:])[0, 1] for l in range(25)]; tau[L] = int(np.argmax(cc)) / 12.0


def mods(L):
    a = np.array(R[str(L)]["amp_by_month"]) * np.exp(1j * np.radians(R[str(L)]["phase_by_month_deg"])); r = a / a.mean() - 1; m = np.arange(12)
    return np.mean(r.real * np.exp(-4j * np.pi * m / 12)), np.mean(r.imag * np.exp(-4j * np.pi * m / 12))


M = {L: mods(L) for L in LV}
w = np.degrees(np.unwrap([np.angle(M[L][0] / M[L][1]) for L in LV])); w -= w[0]
gate = np.mean([M[L][1] for L in LV]); a10 = M[10][0]
ratio = lambda al, ph: abs(al - 1j * ph) / abs(al + 1j * ph)                     # = lower (2-q) / upper (2+q)
out = {}
print(" hPa  delay   winding pred / meas (sd)    lower/upper: observed   model A   model B")
for i, L in enumerate(LV):
    rot = np.exp(-1j * 2 * np.pi * 2 * tau[L]); obs = ratio(*M[L])
    A = ratio(abs(M[L][0]) * np.exp(1j * np.angle(a10)) * rot, gate); B = ratio(a10 * rot, gate)
    out[L] = {"tau_months": tau[L] * 12, "winding_pred_deg": -720 * tau[L], "winding_meas_deg": float(w[i]),
              "winding_sd_deg": R[str(L)]["AM_minus_PM_phase2_sd_deg"], "lower_over_upper_obs": obs, "model_A": A, "model_B": B}
    print(f"{L:4d}  {tau[L] * 12:4.1f} mo   {-720 * tau[L]:6.0f} / {w[i]:5.0f} ({R[str(L)]['AM_minus_PM_phase2_sd_deg']:.0f})"
          f"           {obs:6.2f}     {A:6.2f}    {B:6.2f}")
p = np.array([-720 * tau[L] for L in LV]); print(f"winding: r = {np.corrcoef(p, w)[0, 1]:.2f}, slope meas/pred = {np.dot(w, p) / np.dot(p, p):.2f}")
json.dump(out, open(os.path.join(HERE, "qbo_sideband_descent.json"), "w"), indent=1)
