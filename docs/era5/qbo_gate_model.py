"""Does an equinox sample-and-hold of the draconic lunar factor produce the QBO's column-wide seasonal timing gate?
(2026-10-26)  Observed target (qbo_sideband_rules.json, KIT radiosondes 1953-2026): the QBO timing (phase at
q = 0.4219 cpy, per calendar month) is the same at all levels, most retarded around April, with a semiannual depth of
~0.1 rad; this calendar-locked gate plus the descending response explains the 2 +/- q sideband mix by level
(qbo_sideband_descent.py).

Forcing   beta(t) = sin F, F = Moon's mean argument of latitude (Meeus: 93.2720950 + 483202.0175233 T deg), i.e. the
          sign and size of the Moon's ecliptic latitude with the real draconic phase.
Comb      sample beta at the equinox instants (Meeus mean-equinox formulae, so the 186.4 / 178.8 d spacing is real) and
          hold until the next sample of that comb:
            annual-Mar, annual-Sep (one equinox a year), semi-unsigned (both, same sign),
            semi-signed (+ at the March / northward solar crossing, - at September / southward)
Response  none (pure sample-and-hold) or a first-order low-pass (60 d), then monthly means, 1953-2026
Analysis  identical to qbo_sideband_rules.py: per calendar month u = Re[a e^{i 2 pi q t}] + const + trend; AM = |a|/|abar|,
          PM = arg(a/abar); semiannual depth and phase; month of most retarded timing; correlation of the 12-month
          PM curve with the observed column-mean PM curve; joint-LS sideband amplitudes 1 +/- q, 2 +/- q.
Output    qbo_gate_model.png, qbo_gate_model.json
usage: python3 qbo_gate_model.py"""
import json, os
import datetime as dt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
q = (365.242199 / 27.212221) % 1
OBS = json.load(open(os.path.join(HERE, "qbo_sideband_rules.json")))["levels"]
LV = [str(L) for L in (10, 15, 20, 30, 40, 50, 70)]
mm = np.arange(12)


def jd(d):                                                # datetime -> Julian day
    return d.toordinal() + 1721424.5 + (d.hour + d.minute / 60) / 24


def equinoxes(y):                                         # Meeus ch. 27, mean equinoxes (JDE, years 1000-3000)
    Y = (y - 2000) / 1000
    mar = 2451623.80984 + 365242.37404 * Y + 0.05169 * Y ** 2 - 0.00411 * Y ** 3 - 0.00057 * Y ** 4
    sep = 2451810.21715 + 365242.01767 * Y - 0.11575 * Y ** 2 + 0.00337 * Y ** 3 + 0.00078 * Y ** 4
    return mar, sep


d0 = dt.datetime(1953, 1, 1); days = np.arange((dt.datetime(2026, 10, 1) - d0).days)
JD = jd(d0) + days; T = (JD - 2451545.0) / 36525
beta = np.sin(np.radians(93.2720950 + 483202.0175233 * T))
EQ = [(e, k) for y in range(1952, 2028) for k, e in enumerate(equinoxes(y))]          # k = 0 March, 1 September
spacing = np.mean([equinoxes(y)[1] - equinoxes(y)[0] for y in range(1953, 2027)])
print(f"q = {q:.4f} cpy; mean March->September equinox spacing {spacing:.2f} d")

COMBS = {"annual-Mar": lambda k: 1 if k == 0 else None, "annual-Sep": lambda k: 1 if k == 1 else None,
         "semi-unsigned": lambda k: 1, "semi-signed": lambda k: 1 if k == 0 else -1}


def held(comb):
    x = np.full(len(JD), np.nan); ev = [(e, comb(k)) for e, k in EQ if comb(k) is not None]
    for (e, s), (e2, _) in zip(ev[:-1], ev[1:]):
        val = s * np.sin(np.radians(93.2720950 + 483202.0175233 * (e - 2451545.0) / 36525))
        x[(JD >= e) & (JD < e2)] = val
    return x


def lowpass(x, tau=60.0):
    y = np.empty_like(x); y[0] = 0.0; a = np.exp(-1 / tau)
    for i in range(1, len(x)): y[i] = a * y[i - 1] + (1 - a) * x[i]
    return y


def monthly(x):
    out, tt, mo = [], [], []
    for y in range(1953, 2027):
        for m in range(1, 13):
            a = (dt.datetime(y, m, 1) - d0).days; b = (dt.datetime(y + (m == 12), m % 12 + 1, 1) - d0).days
            if b > len(x): break
            out.append(np.nanmean(x[a:b])); tt.append(y + (m - 0.5) / 12); mo.append(m)
    return np.array(out), np.array(tt), np.array(mo)


def analyse(u, t, month):
    a = np.empty(12, complex)
    for m in range(1, 13):
        s = (month == m) & np.isfinite(u); tt = t[s]
        X = np.column_stack([np.cos(2 * np.pi * q * tt), np.sin(2 * np.pi * q * tt), np.ones(s.sum()), tt - tt.mean()])
        c = np.linalg.lstsq(X, u[s], rcond=None)[0]; a[m - 1] = c[0] - 1j * c[1]
    ab = a.mean(); r = a / ab - 1
    ph2 = np.mean(r.imag * np.exp(-4j * np.pi * mm / 12)); ph1 = np.mean(r.imag * np.exp(-2j * np.pi * mm / 12))
    s = np.isfinite(u); tt = t[s]; fr = [q, 1 - q, 1 + q, 2 - q, 2 + q]
    cols = [np.ones(s.sum()), tt - tt.mean()] + [f(2 * np.pi * k * tt) for k in (1, 2, 3, 4) for f in (np.cos, np.sin)]
    for v in fr: cols += [np.cos(2 * np.pi * v * tt), np.sin(2 * np.pi * v * tt)]
    c = np.linalg.lstsq(np.column_stack(cols), u[s], rcond=None)[0][10:]; amp = [np.hypot(c[2 * i], c[2 * i + 1]) for i in range(5)]
    return {"pm_deg": np.degrees(np.angle(a / ab)), "am": np.abs(a) / abs(ab), "PM2": 2 * abs(ph2), "PM2_phase": float(np.degrees(np.angle(ph2))),
            "PM1": 2 * abs(ph1), "PM1_phase": float(np.degrees(np.angle(ph1))),
            "lines": dict(zip(["q", "1-q", "1+q", "2-q", "2+q"], [a_ / amp[0] for a_ in amp]))}


# observed column-mean gate (each level's PM curve, weighted equally)
pm_obs = np.degrees(np.angle(np.mean([np.exp(1j * np.radians(OBS[L]["phase_by_month_deg"])) for L in LV], axis=0)))
ph2_obs = np.mean(np.radians(pm_obs) * np.exp(-4j * np.pi * mm / 12)); ph1_obs = np.mean(np.radians(pm_obs) * np.exp(-2j * np.pi * mm / 12))
print(f"OBSERVED column gate: most retarded month {int(np.argmin(pm_obs)) + 1}, PM2 {2 * abs(ph2_obs):.3f} rad @ {np.degrees(np.angle(ph2_obs)):+.0f} deg, "
      f"PM1 {2 * abs(ph1_obs):.3f} rad @ {np.degrees(np.angle(ph1_obs)):+.0f} deg")
res = {"observed": {"pm_deg": pm_obs.tolist(), "PM2": float(2 * abs(ph2_obs)), "PM2_phase": float(np.degrees(np.angle(ph2_obs))),
                    "PM1": float(2 * abs(ph1_obs)), "PM1_phase": float(np.degrees(np.angle(ph1_obs)))}, "models": {}}
print(f"{'model':26s} {'retarded':>8s} {'PM2 rad':>8s} {'ph2':>5s} {'PM1 rad':>8s} {'ph1':>5s} {'r(PM curve)':>11s} | lines/q: 1-q 1+q 2-q 2+q")
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].plot(mm + 1, pm_obs, "k-", lw=3, label="observed (KIT, column mean)")
for i, (name, comb) in enumerate(COMBS.items()):
    x = held(comb)
    for resp, xr in (("S&H", x), ("S&H + 60d low-pass", lowpass(np.nan_to_num(x)))):
        u, t, mo = monthly(xr); A = analyse(u, t, mo); r = float(np.corrcoef(A["pm_deg"], pm_obs)[0, 1]); key = f"{name} | {resp}"
        res["models"][key] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in A.items()}; res["models"][key]["r_pm"] = r
        L_ = A["lines"]
        print(f"{key:26s} {int(np.argmin(A['pm_deg'])) + 1:8d} {A['PM2']:8.3f} {A['PM2_phase']:+5.0f} {A['PM1']:8.3f} {A['PM1_phase']:+5.0f} {r:11.2f} | "
              f"{L_['1-q']:.2f} {L_['1+q']:.2f} {L_['2-q']:.2f} {L_['2+q']:.2f}")
        if resp == "S&H":
            ax[0].plot(mm + 1, A["pm_deg"], "o-", color=f"C{i}", lw=1, label=f"{name} (r = {r:.2f})")
            ax[1].plot(mm + 1, A["am"], "o-", color=f"C{i}", lw=1, label=name)
        if name == "semi-signed" and resp == "S&H":
            ax[2].plot(t[:180], u[:180], color="C3", lw=1, label="semi-signed S&H (monthly)")
ax[0].set_xlabel("calendar month"); ax[0].set_ylabel("timing (deg)"); ax[0].set_title("Seasonal timing gate (PM): models vs observed", fontsize=10); ax[0].legend(fontsize=7)
ax[1].set_xlabel("calendar month"); ax[1].set_ylabel("|a(month)| / mean"); ax[1].set_title("Seasonal amplitude (AM), models", fontsize=10); ax[1].legend(fontsize=7)
ax[2].set_xlabel("year"); ax[2].set_title("Forcing: draconic latitude sampled at equinoxes, held", fontsize=10); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "qbo_gate_model.png"), dpi=120)
json.dump(res, open(os.path.join(HERE, "qbo_gate_model.json"), "w"), indent=1); print("wrote qbo_gate_model.png/.json")

# ---- mixed model: comb asymmetry rho = b/a (September / March impulse; -1 signed, 0 March-only, +1 unsigned) and gated
#      fraction g (rest = the ungated q line of the same signal), 60 d response; fitted to the observed 12-month PM curve
print("\nmixed model fit (rho, g), response 60 d:")
best = None
for rho in np.linspace(-1, 1, 21):
    x = held(lambda k, rho=round(rho, 3): 1 if k == 0 else (rho if rho != 0 else None))
    u, t, mo = monthly(lowpass(np.nan_to_num(x), 60)); s = np.isfinite(u)
    X = np.column_stack([np.cos(2 * np.pi * q * t[s]), np.sin(2 * np.pi * q * t[s])]); cont = X @ np.linalg.lstsq(X, u[s], rcond=None)[0]
    for gg in np.linspace(0, 1, 21):
        mix = u.copy(); mix[s] = (1 - gg) * cont + gg * u[s]; B = analyse(mix, t, mo)
        err = float(np.sqrt(np.mean(((B["pm_deg"] - pm_obs + 180) % 360 - 180) ** 2)))
        if best is None or err < best[0]: best = (err, float(rho), float(gg), B)
err, rho, gg, B = best; null = float(np.sqrt(np.mean(pm_obs ** 2)))
print(f"  best rho {rho:+.1f}, g {gg:.2f}: rms {err:.1f} deg (no gate {null:.1f}); retarded month {int(np.argmin(B['pm_deg'])) + 1}; "
      f"PM2 {B['PM2']:.3f} (obs {2 * abs(ph2_obs):.3f}), PM1 {B['PM1']:.3f} (obs {2 * abs(ph1_obs):.3f})")
res["mixed_fit_60d"] = {"rho": rho, "g": gg, "rms_deg": err, "null_rms_deg": null, "pm_deg": B["pm_deg"].tolist(), "PM2": B["PM2"], "PM1": B["PM1"]}
ax[0].plot(mm + 1, B["pm_deg"], "k--", lw=2, label=f"mixed fit: rho {rho:+.1f}, g {gg:.2f} (rms {err:.1f} deg)"); ax[0].legend(fontsize=7)

# ---- control: the gate is a property of the comb timing, not of the lunar phase
print("control, draconic phase shifted (annual-Mar, 60 d):")
base = None; res["control_r"] = []
for off in (0, 77, 151, 233, 300):
    x = np.full(len(JD), np.nan); ev = [e for e, k in EQ if k == 0]
    for e, e2 in zip(ev[:-1], ev[1:]):
        x[(JD >= e) & (JD < e2)] = np.sin(np.radians(93.2720950 + off + 483202.0175233 * (e - 2451545.0) / 36525))
    u, t, mo = monthly(lowpass(np.nan_to_num(x), 60)); A = analyse(u, t, mo); base = A["pm_deg"] if base is None else base
    r = float(np.corrcoef(A["pm_deg"], base)[0, 1]); res["control_r"].append(r)
    print(f"  offset {off:3d} deg: retarded month {int(np.argmin(A['pm_deg'])) + 1}, PM1 {A['PM1']:.3f}, PM2 {A['PM2']:.3f}, r {r:.3f}")
fig.savefig(os.path.join(HERE, "qbo_gate_model.png"), dpi=120)
json.dump(res, open(os.path.join(HERE, "qbo_gate_model.json"), "w"), indent=1)
