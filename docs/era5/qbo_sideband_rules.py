"""Sideband selection rules of the QBO by altitude: why 2+q survives at one level and 2-q at another (2026-10-26).
Data   kit_qbo.dat (FU Berlin / KIT monthly radiosonde winds: Canton 1953-67, Gan 1967-75, Singapore 1976-),
       levels 70, 50, 40, 30, 20, 15, 10 hPa, 0.1 m/s; parsed by fixed value end-columns (flags end two columns later).
Lines  q = 13.4219 mod 1 = 0.4219 cpy (draconic alias). Sidebands k +/- q for k = 1, 2.
Method For each calendar month fit u = Re[a e^{i 2 pi q t}] + const over all years (cos, sin regressors): a(month) is the
       QBO's complex amplitude as a function of season. With B_j = <a(month) e^{+i 2 pi j month/12}>, a term B_j
       e^{-i 2 pi j t} gives the line |q - j|: B_{+2} <-> 2-q (lower), B_{-2} <-> 2+q (upper), B_{+1} <-> 1-q,
       B_{-1} <-> 1+q (checked against a joint least-squares fit of the five lines).
       Writing a = abar (1 + alpha + i phi) (alpha: seasonal amplitude modulation, phi: seasonal timing modulation),
       |B_{+j}| ~ |alpha_j + i phi_j| and |B_{-j}| ~ |alpha_j - i phi_j|: one sideband cancels when the amplitude and
       timing modulations have equal depth and are 90 deg apart; the SIGN of that 90 deg picks which sideband survives
       (empirically: AM - PM phase near +90 deg -> upper 2+q survives, near -90 deg -> lower 2-q survives).
       Uncertainties: block bootstrap over years (whole years resampled), 300 draws.
Output qbo_sideband_rules.png, qbo_sideband_rules.json
usage: python3 qbo_sideband_rules.py [path to kit_qbo.dat]"""
import json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "../../experiments/Oct2026_generator/QBO/kit_qbo.dat")
LEVELS = [70, 50, 40, 30, 20, 15, 10]; ENDS = [16, 23, 30, 37, 44, 51, 58]
q = (365.242199 / 27.212221) % 1

rows = {}
for ln in open(SRC).read().splitlines()[9:]:
    if len(ln) < 20 or not ln[:5].strip().isdigit():
        continue
    yy, mm = int(ln[6:8]), int(ln[8:10]); yr = (1900 if yy >= 50 else 2000) + yy
    vals = {m.end() + 10: int(m.group()) for m in re.finditer(r"-?\d+", ln[10:])}
    rows[(yr, mm)] = [vals.get(e, np.nan) / 10.0 for e in ENDS]
keys = sorted(rows); t = np.array([y + (m - 0.5) / 12 for y, m in keys]); U = np.array([rows[k] for k in keys]); month = np.array([m for _, m in keys])
print(f"{len(t)} months {t[0]:.2f}-{t[-1]:.2f}; q = {q:.4f} cpy ({1 / q:.3f} yr)")
chk = np.loadtxt(os.path.join(os.path.dirname(SRC), "kit_qbo_u30.dat"))
common = np.intersect1d(np.round(chk[:, 0], 3), np.round(t - 0.5 / 12, 3))
a30 = dict(zip(np.round(t - 0.5 / 12, 3), U[:, 3])); b30 = dict(zip(np.round(chk[:, 0], 3), chk[:, 1]))
print(f"parse check vs kit_qbo_u30.dat: max |diff| over {len(common)} months = {max(abs(a30[c] - b30[c]) for c in common if np.isfinite(a30[c])):.2f} m/s")


mm = np.arange(12)


def per_month(u, tt=None, mo=None):
    tt_all = t if tt is None else tt; mo_all = month if mo is None else mo
    a = np.empty(12, complex)
    for m in range(1, 13):
        s = (mo_all == m) & np.isfinite(u); tt = tt_all[s]
        X = np.column_stack([np.cos(2 * np.pi * q * tt), np.sin(2 * np.pi * q * tt), np.ones(s.sum()), tt - tt.mean()])
        c = np.linalg.lstsq(X, u[s], rcond=None)[0]; a[m - 1] = c[0] - 1j * c[1]        # u = Re[a e^{i 2 pi q t}]
    return a


def joint(u):
    s = np.isfinite(u); tt = t[s]; fr = {"q": q, "1-q": 1 - q, "1+q": 1 + q, "2-q": 2 - q, "2+q": 2 + q}
    cols = [np.ones(s.sum()), tt - tt.mean()] + [f(2 * np.pi * k * tt) for k in (1, 2, 3, 4) for f in (np.cos, np.sin)]
    for v in fr.values(): cols += [np.cos(2 * np.pi * v * tt), np.sin(2 * np.pi * v * tt)]
    c = np.linalg.lstsq(np.column_stack(cols), u[s], rcond=None)[0][10:]
    return {k: float(np.hypot(c[2 * i], c[2 * i + 1])) for i, k in enumerate(fr)}


def stats(u, tt_mask=None):
    a = per_month(u); ab = a.mean(); r = a / ab - 1
    al2 = np.mean(r.real * np.exp(-2j * np.pi * 2 * mm / 12)); ph2 = np.mean(r.imag * np.exp(-2j * np.pi * 2 * mm / 12))
    up = abs(np.mean(a * np.exp(-2j * np.pi * 2 * mm / 12))); lo = abs(np.mean(a * np.exp(2j * np.pi * 2 * mm / 12)))
    return up / lo, float(np.degrees(np.angle(al2 / ph2)))


rng = np.random.default_rng(1); years = np.floor(t).astype(int); uy = np.unique(years)
res = {"q_cpy": q, "levels": {}}
fig, ax = plt.subplots(1, 4, figsize=(20, 4.8))
print(f"{'hPa':>4s} {'|abar|':>7s} | joint LS amp (m/s): {'1-q':>6s} {'1+q':>6s} {'2-q':>6s} {'2+q':>6s} | |B-2|/|B+2| | "
      f"semiannual AM depth, PM depth (rad), AM-PM phase (deg) +/- bootstrap sd")
for j, L in enumerate(LEVELS):
    u = U[:, j]
    if np.isfinite(u).sum() < 300:
        continue
    a = per_month(u); ab = a.mean(); r = a / ab - 1                      # r = alpha + i phi (small-modulation form)
    alpha, phi = r.real, r.imag
    Bj = {k: np.mean(a * np.exp(2j * np.pi * k * mm / 12)) for k in (-2, -1, 1, 2)}
    boot = []
    for _ in range(300):
        ys = rng.choice(uy, len(uy)); idx = np.concatenate([np.where(years == y)[0] for y in ys])
        # resampled years keep their own times so the q phase is preserved
        aa = per_month(u[idx], t[idx], month[idx]); r_ = aa / aa.mean() - 1
        a2 = np.mean(r_.real * np.exp(-2j * np.pi * 2 * mm / 12)); p2 = np.mean(r_.imag * np.exp(-2j * np.pi * 2 * mm / 12))
        boot.append((abs(np.mean(aa * np.exp(-2j * np.pi * 2 * mm / 12))) / abs(np.mean(aa * np.exp(2j * np.pi * 2 * mm / 12))),
                     np.angle(a2 / p2)))
    boot = np.array(boot); lr = np.log10(boot[:, 0].real); ph_b = boot[:, 1].real
    circ_sd = float(np.degrees(np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * ph_b)))))))
    al2 = np.mean(alpha * np.exp(-2j * np.pi * 2 * mm / 12)); ph2 = np.mean(phi * np.exp(-2j * np.pi * 2 * mm / 12))
    al1 = np.mean(alpha * np.exp(-2j * np.pi * mm / 12)); ph1 = np.mean(phi * np.exp(-2j * np.pi * mm / 12))
    J = joint(u)
    rel2 = float(np.degrees(np.angle(al2 / ph2))); rel1 = float(np.degrees(np.angle(al1 / ph1)))
    up_lo = abs(Bj[-2]) / abs(Bj[2]); lo95, hi95 = np.percentile(lr, [2.5, 97.5])
    res["levels"][L] = {"abar": float(abs(ab)), "joint": J, "B": {str(k): float(abs(v)) for k, v in Bj.items()},
                        "upper_over_lower": float(up_lo), "upper_over_lower_95": [float(10 ** lo95), float(10 ** hi95)],
                        "AM_minus_PM_phase2_sd_deg": circ_sd,
                        "AM2": float(2 * abs(al2)), "PM2_rad": float(2 * abs(ph2)), "AM_minus_PM_phase2_deg": rel2,
                        "AM1": float(2 * abs(al1)), "PM1_rad": float(2 * abs(ph1)), "AM_minus_PM_phase1_deg": rel1,
                        "amp_by_month": np.abs(a).tolist(), "phase_by_month_deg": np.degrees(np.angle(a / ab)).tolist()}
    print(f"{L:4d} {abs(ab):7.2f} | {'':19s}{J['1-q']:6.2f} {J['1+q']:6.2f} {J['2-q']:6.2f} {J['2+q']:6.2f} | "
          f"up/lo {up_lo:5.2f} [{10 ** lo95:.2f}-{10 ** hi95:.2f}] | AM2 {2 * abs(al2):.3f}  PM2 {2 * abs(ph2):.3f}  rel {rel2:+5.0f} +/- {circ_sd:3.0f}  "
          f"(annual: AM1 {2 * abs(al1):.3f} PM1 {2 * abs(ph1):.3f} rel {rel1:+5.0f})")
    col = plt.cm.viridis(j / (len(LEVELS) - 1))
    ax[0].plot(mm + 1, np.abs(a) / abs(ab), "o-", color=col, label=f"{L} hPa")
    ax[1].plot(mm + 1, np.degrees(np.angle(a / ab)), "o-", color=col, label=f"{L} hPa")
    ax[2].bar(j - 0.2, J["2-q"], 0.4, color="C0", label="2 - q (1.578 cpy)" if j == 0 else None)
    ax[2].bar(j + 0.2, J["2+q"], 0.4, color="C3", label="2 + q (2.422 cpy)" if j == 0 else None)
ax[0].set_xlabel("calendar month"); ax[0].set_ylabel("|a(month)| / |mean|"); ax[0].set_title("QBO amplitude by season (AM)", fontsize=10)
ax[1].set_xlabel("calendar month"); ax[1].set_ylabel("deg"); ax[1].set_title("QBO phase by season (timing / PM)", fontsize=10); ax[1].legend(fontsize=7)
lv = [L for L in LEVELS if L in res["levels"]]; ax[2].set_xticks(range(len(lv))); ax[2].set_xticklabels([f"{L}" for L in lv])
ax[2].set_xlabel("hPa"); ax[2].set_ylabel("line amplitude (m/s)"); ax[2].legend(fontsize=8); ax[2].set_title("2 - q vs 2 + q sideband by level (joint LS)", fontsize=10)
Z = {70: 18.5, 50: 20.6, 40: 22.0, 30: 23.9, 20: 26.5, 15: 28.5, 10: 31.0}                     # approx. altitude, km
ph = np.unwrap(np.radians([res["levels"][L]["AM_minus_PM_phase2_deg"] for L in lv[::-1]]))   # from 10 hPa downward
zz = [Z[L] for L in lv[::-1]]; sd = [res["levels"][L]["AM_minus_PM_phase2_sd_deg"] for L in lv[::-1]]
ax[3].errorbar(np.degrees(ph), zz, xerr=sd, fmt="o-", color="k", capsize=3)
for x in (-630, -450, -270, -90, 90):          # +90 (mod 360): upper 2+q survives (red); -90: lower 2-q (blue)
    ax[3].axvline(x, color="C3" if (x - 90) % 360 == 0 else "C0", lw=0.8, ls=":")
for L, x, z in zip(lv[::-1], np.degrees(ph), zz): ax[3].annotate(f"{L} hPa", (x, z), xytext=(5, -3), textcoords="offset points", fontsize=7)
ax[3].set_xlabel("AM - PM semiannual phase, unwrapped from 10 hPa (deg)"); ax[3].set_ylabel("altitude (km, approx.)")
ax[3].set_title("Winding of the AM-PM phase with depth\n(dotted: quadrature; red = upper 2+q only, blue = lower 2-q only)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "qbo_sideband_rules.png"), dpi=120)
json.dump(res, open(os.path.join(HERE, "qbo_sideband_rules.json"), "w"), indent=1); print("wrote qbo_sideband_rules.png/.json")
