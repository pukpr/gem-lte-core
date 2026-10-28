#!/usr/bin/env python3
"""The QBO read as a control problem (2026-10-26). Discussion: qbo_control.md.
Supplementary view after the Chandler-wobble counter-torque analysis (github.com/pukpr/ChandlerWobble, docs/
cw_counter_torque.py): run the plant backwards to find the input that produced the observed behaviour, budget it
against the lunar forcing, and see what a feedback controller acting at the same instants would have to do.

Plant   the QBO tokamak device of qbo_tokamak.py, linear-edge variant (COMB|linear, parameters from qbo_tokamak.json):
        pulses a_n at a semiannual comb (trigger day of year + half a year) -> leaky integrator (tau_i) -> edge ->
        advection-diffusion penetration (w, eta, depth d, damping tau_d) to 30 / 50 / 70 hPa. Linear, so each level's
        wind is u_L(t) = c_L sum_n a_n k_L(t - t_n) + const.
1. Inverse plant in the comb basis: ridge-regularised least squares for the pulse amplitudes a_req,n (one per half
   year, 1940-2026) that reproduce the observed winds at all three levels together (KIT radiosondes, standardised);
   ridge parameter by generalised cross-validation.
2. Compare with the astronomical pulses a_ast,n = sin(beta_moon) sin(dec_sun) (solar distance)^3 at the same
   instants (Meeus elements; same functions as qbo_tokamak.py): pulse-by-pulse correlation by era, and the two pulse
   trains' components at the 2.37-yr alias (0.4219 cpy in the half-year-sampled train).
3. Residual input r_n = a_req,n - kappa a_ast,n: what the "controller" (nature) had to add beyond the lunar comb, in
   particular across the 2015-16 and 2019-20 disruptions.
4. Budget: the zonal-mean zonal component of the degree-2 lunar tidal force is identically zero (the force is the
   gradient of a potential with no zonal-mean longitude dependence), whereas the QBO needs several m/s per month at
   its transitions; so any lunar role is a gate on wave-driven momentum, not a momentum source.
5. Control of the device from 1990: feedback u_n = -g x (integrator state just before the pulse) at the comb instants,
   against feedforward u_n = -f a_ast,n (the controller knows the ephemeris). Unlike the high-Q Chandler wobble, the
   QBO device is kick-dominated (half its state decays each half year), so feedback barely acts and only feedforward
   cancels it (exact cancellation at f = 1 is an identity for a linear plant).
    python3 qbo_control.py   -> qbo_control.png, qbo_control.json
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import fftconvolve

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
LEGS = HERE.parent / "experiments/Oct2026_generator/legs"
FILES = {30: (LEGS / "qbo30_point0E/kit_qbo_u30.dat", 1.0), 50: (LEGS / "qbo50/qbo50.dat", 100.0), 70: (LEGS / "qbo70_point0E/kit_qbo_u70.dat", 1.0)}
DZ = {30: 0.0, 50: 3.5, 70: 5.9}
DT = 1 / 73.0
TG = np.arange(1900.0, 2027.0, DT)
ALIAS = 365.2422 / 27.212221 % 1                         # 0.4219 cpy
rad = np.radians
P = json.load(open(HERE / "qbo_tokamak.json"))["COMB|linear"]["x"]
DOY, _, W, ETA, D0, TAU_D, TAU_I = P


def astro(t):                                            # as qbo_tokamak.astro (Meeus mean elements)
    T = (t - 2000.0) * 365.25 / 36525.0
    D, M, Mp, F = (rad(v) for v in (297.8501921 + 445267.1114034 * T, 357.5291092 + 35999.0502909 * T,
                                    134.9633964 + 477198.8675055 * T, 93.2720950 + 483202.0175233 * T))
    beta = rad(5.128) * np.sin(F) + rad(0.281) * np.sin(Mp + F) + rad(0.278) * np.sin(Mp - F) + rad(0.173) * np.sin(2 * D - F)
    lams = rad(280.466 + 36000.770 * T) + rad(1.915) * np.sin(M) + rad(0.020) * np.sin(2 * M)
    return np.sin(beta), np.sin(rad(23.439)) * np.sin(lams), (1 + 0.0167 * np.cos(M)) ** 3


def kernel(L):
    """response at level L to a unit pulse at t = 0: leaky integrator then first-passage penetration"""
    s = np.arange(0, int(15 / DT)) * DT; integ = np.exp(-s / TAU_I)
    z = D0 + DZ[L]; ss = s[1:]
    g = np.concatenate([[0.0], z / np.sqrt(4 * np.pi * ETA * ss ** 3) * np.exp(-(z - W * ss) ** 2 / (4 * ETA * ss) - ss / TAU_D) * DT])
    return fftconvolve(integ, g)[: len(s)]


def main():
    tn = np.array([y + DOY / 365.2422 + h for y in range(1940, 2027) for h in (0.0, 0.5)]); tn = tn[tn < 2026.8]
    sb, sds, dist = astro(tn); a_ast = sb * sds * dist
    K = {L: kernel(L) for L in FILES}
    obs = {}
    for L, (f, sc) in FILES.items():
        a = np.loadtxt(f); a = a[np.isfinite(a[:, 1]) & (a[:, 0] >= 1953)]; obs[L] = (a[:, 0], a[:, 1] * sc)
    sd = {L: obs[L][1].std() for L in obs}

    def design(L):                                        # columns: response at the observation times to each pulse
        t, _ = obs[L]; X = np.zeros((len(t), len(tn)))
        for j, t0 in enumerate(tn):
            lag = t - t0; m = (lag >= 0) & (lag < 15)
            X[m, j] = np.interp(lag[m], np.arange(len(K[L])) * DT, K[L])
        return X
    Xs = {L: design(L) for L in obs}
    # alternating LS: per-level gain c_L (c_30 = 1) and pulses a; ridge by GCV on the stacked system
    c = {L: 1.0 for L in obs}
    for it in range(4):
        A = np.vstack([c[L] * Xs[L] for L in obs]); y = np.concatenate([(obs[L][1] - obs[L][1].mean()) / sd[L] for L in obs])
        A = A - A.mean(0)
        U, S, Vt = np.linalg.svd(A, full_matrices=False); uy = U.T @ y
        best = None
        for lam in np.logspace(-4, 2, 61) * S[0] ** 2:
            fil = S ** 2 / (S ** 2 + lam); res = y - U @ (fil * uy); gcv = len(y) * np.sum(res ** 2) / (len(y) - fil.sum()) ** 2
            if best is None or gcv < best[0]: best = (gcv, lam, fil)
        a_req = Vt.T @ (best[2] / S * uy)
        for L in obs:
            p = Xs[L] @ a_req; p -= p.mean(); yy = (obs[L][1] - obs[L][1].mean()) / sd[L]
            c[L] = float(np.dot(p, yy) / np.dot(p, p)) if L != 30 else 1.0
    fitr = {L: float(np.corrcoef(Xs[L] @ a_req, obs[L][1])[0, 1]) for L in obs}
    dof = float(best[2].sum())
    print(f"inverse plant: {len(tn)} pulses, effective dof {dof:.0f}; reconstruction r 30/50/70 = "
          + " / ".join(f"{fitr[L]:.2f}" for L in obs) + f"; gains c = {c}")

    # 2. compare with astronomical pulses
    out = {"plant": dict(zip(["doy", "theta", "w", "eta", "d", "tau_d", "tau_i"], P)), "pulses": len(tn), "dof": dof,
           "reconstruction_r": fitr, "eras": {}}
    for a0, a1 in ((1955, 2006), (2006, 2026.8), (1955, 2026.8)):
        m = (tn >= a0) & (tn < a1); r = float(np.corrcoef(a_req[m], a_ast[m])[0, 1])
        # surrogate: shift the astronomical train by a random number of pulses (keeps its spectrum)
        rs = [np.corrcoef(a_req[m], np.roll(a_ast, k)[m])[0, 1] for k in range(3, len(tn) - 3)]
        pct = float(np.mean(np.abs(rs) < abs(r)) * 100)
        out["eras"][f"{a0}-{a1:.0f}"] = {"r": r, "shift_percentile": pct}
        print(f"pulse-by-pulse r(a_req, a_ast) {a0}-{a1:.0f}: {r:+.2f}  (|r| percentile vs shifted trains {pct:.0f})")

    def comp(a, m, f=ALIAS):
        return np.sum(a[m] * np.exp(-2j * np.pi * f * tn[m])) / m.sum()
    for a0, a1 in ((1955, 2006), (2006, 2026.8)):
        m = (tn >= a0) & (tn < a1); R = comp(a_req, m) / comp(a_ast, m)
        out["eras"][f"{a0}-{a1:.0f}"]["alias_ratio"] = [float(abs(R)), float(np.degrees(np.angle(R)))]
        print(f"  alias component a_req / a_ast {a0}-{a1:.0f}: |ratio| {abs(R) / abs(comp(a_req, (tn >= 1955) & (tn < 2006)) / comp(a_ast, (tn >= 1955) & (tn < 2006))):.2f} (rel. to 1955-2005), phase {np.degrees(np.angle(R)):+.0f} deg")
    # 3. residual input
    m0 = (tn >= 1955) & (tn < 2006); kap = float(np.dot(a_req[m0], a_ast[m0]) / np.dot(a_ast[m0], a_ast[m0]))
    r_n = a_req - kap * a_ast; s0 = r_n[m0].std()
    big = [(float(t), float(v / s0)) for t, v in zip(tn, r_n) if abs(v) > 2 * s0 and t >= 1955]
    out.update(kappa=kap, residual_sd_1955_2005=float(s0), large_residual_pulses=big,
               pulses={"t": tn.tolist(), "a_req": a_req.tolist(), "a_ast": a_ast.tolist()})
    print(f"residual input r_n = a_req - {kap:.2f} a_ast; pulses beyond 2 sd (1955-): " + ", ".join(f"{t:.2f} ({v:+.1f})" for t, v in big))
    rms_era = {e: float(np.sqrt(np.mean(r_n[(tn >= a) & (tn < b)] ** 2)) / s0) for e, (a, b) in {"1955-2005": (1955, 2006), "2006-2014": (2006, 2015), "2015-2021": (2015, 2021.5), "2021-2026": (2021.5, 2026.8)}.items()}
    out["residual_rms_rel"] = rms_era; print("residual rms relative to 1955-2005:", {k: round(v, 2) for k, v in rms_era.items()})

    # 4. budget
    t30, u30 = obs[30]; dudt = np.abs(np.diff(u30)) / np.diff(t30) / 12        # m/s per month
    GM, dm, a_e = 6.674e-11 * 7.346e22, 3.844e8, 6.371e6
    tid = 1.5 * GM * a_e / dm ** 3; per_month = tid * 30.44 * 86400
    out["budget"] = {"qbo_du_dt_p90_m_s_per_month": float(np.percentile(dudt, 90)), "qbo_du_dt_median": float(np.median(dudt)),
                     "lunar_tidal_horizontal_accel_max_m_s2": tid, "same_as_m_s_per_month_if_unidirectional": per_month,
                     "zonal_mean_zonal_component": 0.0}
    print(f"budget: QBO |du/dt| at 30 hPa median {np.median(dudt):.1f}, 90th pct {np.percentile(dudt, 90):.1f} m/s per month; "
          f"lunar tidal horizontal acceleration max {tid:.1e} m/s2 = {per_month:.1f} m/s per month if it never reversed; "
          f"zonal-mean zonal component exactly 0")

    # 5. control of the device (astronomical pulses, linear plant), as in cw_counter_torque.pulsed_control: feedback on
    #    the measured state at the comb instants, against feedforward that knows the coming lunar pulse (the ephemeris)
    idx = np.minimum(np.round((tn - TG[0]) / DT).astype(int), len(TG) - 1); kick = dict(zip(idx, a_ast)); rr = np.exp(-DT / TAU_I)

    def run(gain, ff=0.0, on=1990.0):
        acc = 0.0; Dv = np.empty(len(TG))
        for i in range(len(TG)):
            acc *= rr
            if i in kick:
                u = (-gain * acc - ff * kick[i]) if TG[i] >= on else 0.0
                acc += kick[i] + u
            Dv[i] = acc
        return Dv
    zk = np.arange(1, int(15 / DT)) * DT
    gpen = np.concatenate([[0.0], D0 / np.sqrt(4 * np.pi * ETA * zk ** 3) * np.exp(-(D0 - W * zk) ** 2 / (4 * ETA * zk) - zk / TAU_D) * DT])
    pen = lambda Dv: fftconvolve(Dv, gpen)[: len(TG)]
    m = TG > 2000; u_free = pen(run(0.0)); s_free = np.std(u_free[m])
    gains = [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
    fb = {g: float(np.std(pen(run(g))[m]) / s_free) for g in gains}
    ffr = {f: float(np.std(pen(run(0.0, f))[m]) / s_free) for f in (0.25, 0.5, 0.75, 1.0)}
    u_fb, u_ff = pen(run(0.8)), pen(run(0.0, 0.75))
    out["control"] = {"state_surviving_half_year": float(np.exp(-0.5 / TAU_I)), "feedback_residual_by_gain": fb,
                      "feedforward_residual_by_fraction": ffr}
    print(f"control: state surviving a half year {np.exp(-0.5 / TAU_I):.2f}; feedback residual by gain {fb}; feedforward residual {ffr}")
    late = (u_fb, u_ff)
    (HERE / "qbo_control.json").write_text(json.dumps(out, indent=1, default=float))
    plot(tn, a_req, a_ast, kap, r_n, s0, obs, Xs, c, sd, u_free, late, out)


def plot(tn, a_req, a_ast, kap, r_n, s0, obs, Xs, c, sd, u_free, late, out):
    fig = plt.figure(figsize=(15, 13)); gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 0.9, 1])
    ax = fig.add_subplot(gs[0, :]); t, y = obs[30]
    ax.plot(t, (y - y.mean()) / sd[30], color="k", lw=1.1, label="observed 30 hPa (sd units)")
    p = Xs[30] @ a_req; ax.plot(t, p - p.mean(), color="#d94a4a", lw=1, label=f"plant driven by the required pulses (r {out['reconstruction_r'][30]:.2f})")
    ax.set_xlim(1953, 2027); ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.set_title("(a) Running the plant backwards: the semiannual pulse train that reproduces the observed QBO (30/50/70 hPa jointly)", fontsize=10)
    ax = fig.add_subplot(gs[1, :]); sc = 1 / np.std(a_req[tn > 1955])
    ax.vlines(tn, 0, a_req * sc, color="#d94a4a", lw=2, label="required pulses a_req")
    ax.plot(tn + 0.06, kap * a_ast * sc, "o", ms=3, color="#2a78d6", label=f"astronomical pulses x {kap:.2f} (sin β · sin δ_sun · dist³)")
    ax.axvspan(2015.8, 2016.6, color="0.85"); ax.axvspan(2019.8, 2020.4, color="0.85")
    e = out["eras"]
    ax.set_title(f"(b) Required vs astronomical pulse amplitudes: r = {e['1955-2006']['r']:+.2f} (1955–2005), {e['2006-2027']['r']:+.2f} (2006–2026)", fontsize=10)
    ax.set_xlim(1953, 2027); ax.legend(fontsize=8, loc="upper left", ncol=2); ax.set_ylabel("pulse (norm.)")
    ax = fig.add_subplot(gs[2, :])
    ax.vlines(tn, 0, r_n / s0, color="0.3", lw=2); ax.axhline(2, color="C3", ls=":", lw=0.8); ax.axhline(-2, color="C3", ls=":", lw=0.8)
    ax.axvspan(2015.8, 2016.6, color="0.85"); ax.axvspan(2019.8, 2020.4, color="0.85"); ax.set_xlim(1953, 2027)
    ax.set_ylabel("residual (sd of 1955–2005)"); ax.set_xlabel("year")
    rr = out["residual_rms_rel"]
    ax.set_title(f"(c) Extra input beyond the lunar comb, r_n = a_req − κ a_ast (grey: 2016 and 2020 disruptions). rms vs 1955–2005: "
                 f"2006–14 {rr['2006-2014']:.2f}, 2015–21 {rr['2015-2021']:.2f}, 2021–26 {rr['2021-2026']:.2f}", fontsize=10)
    u_fb, u_ff = late; C = out["control"]
    ax = fig.add_subplot(gs[3, 0]); m = TG > 1975; sc = 1 / np.std(u_free[TG > 1965])
    ax.plot(TG[m], u_free[m] * sc, color="0.6", lw=0.8, label="device, lunar comb only")
    ax.plot(TG[m], u_fb[m] * sc, color="#d94a4a", lw=1, label=f"+ feedback at the comb instants (gain 0.8): {100 * C['feedback_residual_by_gain'][0.8]:.0f}% left")
    ax.plot(TG[m], u_ff[m] * sc, color="#2a78d6", lw=1, label=f"+ feedforward, 75% of the coming lunar pulse: {100 * C['feedforward_residual_by_fraction'][0.75]:.0f}% left")
    ax.axvline(1990, color="k", ls=":"); ax.legend(fontsize=7.5, loc="upper right"); ax.set_xlabel("year")
    ax.set_title("(d) Control from 1990: feedback barely acts, feedforward (ephemeris) does", fontsize=10)
    ax = fig.add_subplot(gs[3, 1]); g = sorted(C["feedback_residual_by_gain"]); f = sorted(C["feedforward_residual_by_fraction"])
    ax.plot(g, [C["feedback_residual_by_gain"][k] for k in g], "o-", color="#d94a4a", label="feedback gain (state just before the pulse)")
    ax.plot(f, [C["feedforward_residual_by_fraction"][k] for k in f], "s-", color="#2a78d6", label="feedforward fraction of the lunar pulse")
    ax.set_xlabel("gain / fraction"); ax.set_ylabel("residual QBO (rms, after 2000)"); ax.set_ylim(0, 1.15); ax.legend(fontsize=8)
    ax.set_title(f"(e) Why: only {100 * C['state_surviving_half_year']:.0f}% of the state survives a half year, so the signal is\n"
                 "the newest kick, which feedback at the same instant cannot see (Chandler: high-Q, feedback works)", fontsize=9.5)
    fig.tight_layout(); fig.savefig(HERE / "qbo_control.png", dpi=105); plt.close(fig)


if __name__ == "__main__":
    main()
