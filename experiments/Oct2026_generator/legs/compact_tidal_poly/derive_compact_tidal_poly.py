#!/usr/bin/env python3
"""Derive & benchmark compact, Meeus-free gravitational polynomial representations of dLOD.

Benchmarks on Feb2026/dlod3.dat (20,829 daily samples, 1962-2019) against:
  - 52-line empirical lpap (104 free params): CC ~ 0.9717
  - Compact Model A (Kepler-Hill closed-form gravitational potential, 4-8 physical params, ZERO Meeus tables):
      Phi(t) = u(t)^3 * (3*sin^2(delta(t)) - 1) + w_S * (3*sin^2(delta_S(t)) - 1)
      with u(t) = 1 + e_r cos(M') + ev_r cos(2D - M') + var_r cos(2D)
           lam  = L + e_l sin(M') + ev_l sin(2D - M') + var_l sin(2D)
           sin(delta) = sin(eps)*sin(lam) + cos(eps)*sin(i)*sin(F + dlam)
  - Compact Model B (Pure Polynomial in 4 fundamental orbital angles L, M', F, D — LINEAR in coefficients!):
      Degree-2/3 polynomial basis built directly from (cos 2L, sin 2L), (cos M', sin M'),
      (cos Omega, sin Omega), (cos 2D, sin 2D), (cos 2Ls, sin 2Ls) — greedy forward selection
      showing exact CC vs parameter count from 2 up to 16 terms!
  - Climate transfer test: fit the same compact polynomial to R5_atl4_cold_H10_l03's fitted tidal sum!
"""
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

ROOT = Path("/home/paul/gemini/gem-lte-core/experiments")
DLOD_PATH = ROOT / "Feb2026/dlod3.dat"
R5_DIR = ROOT / "Oct2026_generator/legs/baltic_3pt_overnight/repro_0693/R5_atl4_cold_H10_l03"
OUT_DIR = ROOT / "Oct2026_generator/legs/compact_tidal_poly"
OUT_DIR.mkdir(parents=True, exist_ok=True)

dlod = np.loadtxt(DLOD_PATH)
t_yr, y_dlod = dlod[:, 0], dlod[:, 1]

# Exact linear fundamental orbital arguments (Julian centuries T from J2000.0 = JD 2451545.0)
# Anchor matching gem-zonal.adb: 1990.5 <-> JD 2448074.761995, tropical year = 365.2422484 d
D2R = math.pi / 180.0
YR_D = 365.2422484

def fundamental_angles(t_decimal_yr, tau_days=0.0):
    """Return ONLY the 4 linear fundamental lunar/solar angles (radians):
    L (mean tropical longitude), Mp (mean anomaly), F (argument of latitude), D (mean elongation),
    plus Ls = L - D (Sun mean longitude) and Omega = L - F (ascending node).
    Purely linear in time (standard IAU/Chapront mean rates) — ZERO Meeus perturbation tables!
    """
    jd = 2448074.761995 + (t_decimal_yr - 1990.5) * YR_D - tau_days
    T = (jd - 2451545.0) / 36525.0
    # Linear mean angles in degrees (epoch J2000.0 + linear rate only!)
    L  = (218.3164477 + 481267.88123421 * T) * D2R
    Mp = (134.9633964 + 477198.8675055  * T) * D2R
    F  = ( 93.2720950 + 483202.0175233  * T) * D2R
    D  = (297.8501921 + 445267.1114034  * T) * D2R
    Ls = L - D
    Om = L - F
    return L, Mp, F, D, Ls, Om

def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])

# ── MODEL A: Compact Closed-Form Kepler-Hill Gravitational Potential (ZERO Meeus!) ──
# 4-term minimal (Kepler only: e_r, e_l, tau, w_S) vs 8-term Kepler+Hill (add evection & variation)
def eval_model_a(params, t_yr, mode="kepler_hill"):
    eps = 23.439291 * D2R
    inc = 5.145396 * D2R
    if mode == "kepler_only":
        # 4 physical params: tau, e_r, e_l, w_S
        tau, e_r, e_l, w_S = params
        ev_r = var_r = ev_l = var_l = 0.0
    elif mode == "kepler_evect":
        # 6 physical params: tau, e_r, e_l, ev_r, ev_l, w_S
        tau, e_r, e_l, ev_r, ev_l, w_S = params
        var_r = var_l = 0.0
    else:  # kepler_hill (8 physical params)
        tau, e_r, e_l, ev_r, ev_l, var_r, var_l, w_S = params

    L, Mp, F, D, Ls, _ = fundamental_angles(t_yr, tau_days=tau)
    # Inverse distance a/r (3-term Kepler+Hill, cubed in gravitational tidal potential!)
    u = 1.0 + e_r * np.cos(Mp) + ev_r * np.cos(2.0 * D - Mp) + var_r * np.cos(2.0 * D)
    # Ecliptic longitude perturbation dlam (3-term Kepler+Hill)
    dlam = e_l * np.sin(Mp) + ev_l * np.sin(2.0 * D - Mp) + var_l * np.sin(2.0 * D)
    lam = L + dlam
    # Exact spherical trig for sin(declination) with inclined lunar orbit
    sin_b = math.sin(inc) * np.sin(F + dlam)
    cos_b = np.sqrt(np.maximum(0.0, 1.0 - sin_b * sin_b))
    sin_dec = math.sin(eps) * cos_b * np.sin(lam) + math.cos(eps) * sin_b
    # Zonal P2^0 lunar + solar gravitational potential
    phi_m = (u ** 3) * (3.0 * sin_dec * sin_dec - 1.0)
    sin_dec_s = math.sin(eps) * np.sin(Ls)
    phi_s = 3.0 * sin_dec_s * sin_dec_s - 1.0
    raw = phi_m + w_S * phi_s
    # dlod3.dat has the 18.61-yr nodal drift (cos Omega) removed in the reference calibration:
    Om = L - F
    X_nod = np.column_stack([np.ones_like(raw), np.cos(Om), np.sin(Om)])
    b_n, *_ = np.linalg.lstsq(X_nod, raw, rcond=None)
    return raw - X_nod @ b_n

print("=" * 105)
print("MODEL A: Compact Closed-Form Kepler-Hill Gravitational Potential (ZERO Meeus Tables)")
print("=" * 105)

# Scan tau_days over [-14, +14] (full Mm cycle, checking BOTH signs of correlation) to find true epoch+lag!
taus = np.linspace(-14.0, 14.0, 561)
best_tau, best_c0 = 0.0, -1.0
for tau_try in taus:
    c_try = abs(corr(eval_model_a([tau_try, 0.0549, 0.1098, 0.0100, 0.0222, 0.0080, 0.0115, 0.08], t_yr, "kepler_hill"), y_dlod))
    if c_try > best_c0:
        best_c0, best_tau = c_try, float(tau_try)
phi_phys0 = eval_model_a([best_tau, 0.0549, 0.1098, 0.0, 0.0, 0.0, 0.0, 0.08], t_yr, "kepler_hill")
phi_phys1 = eval_model_a([best_tau, 0.0549, 0.1098, 0.0100, 0.0222, 0.0080, 0.0115, 0.08], t_yr, "kepler_hill")
print(f"Best lag alignment tau = {best_tau:+.3f} d (epoch + tidal lag)")
print(f"Pure textbook Kepler (0 fitted orbital params, e=0.0549):                CC = {abs(corr(phi_phys0, y_dlod)):+.5f}")
print(f"Pure textbook Kepler+Evect+Var (0 fitted orbital params, 3 terms):       CC = {abs(corr(phi_phys1, y_dlod)):+.5f}")

model_a_res = {}
for mode, p0, npar in [
    ("kepler_only",  [best_tau, 0.0549, 0.1098, 0.46], 4),
    ("kepler_evect", [best_tau, 0.0549, 0.1098, 0.010, 0.022, 0.46], 6),
    ("kepler_hill",  [best_tau, 0.0549, 0.1098, 0.010, 0.022, 0.008, 0.011, 0.46], 8),
]:
    res = minimize(lambda p: -abs(corr(eval_model_a(p, t_yr, mode), y_dlod)), p0, method="Nelder-Mead",
                   options={"maxiter": 600, "xatol": 1e-4, "fatol": 1e-6})
    y_m = eval_model_a(res.x, t_yr, mode)
    c = corr(y_m, y_dlod)
    # Sign-align
    if c < 0: y_m = -y_m; c = -c
    model_a_res[mode] = {"npar": npar, "cc": c, "params": [float(x) for x in res.x]}
    print(f"Fitted {mode:14s} ({npar} physical params):                          CC = {c:+.5f}  params={[round(x,4) for x in res.x]}")

# ── MODEL B: Pure Polynomial Expansion in Fundamental Orbital Coordinates (LINEAR in coeffs!) ──
print("\n" + "=" * 105)
print("MODEL B: Compact Gravitational Polynomial Basis in (2L, M', Omega, 2D, 2Ls) — Linear in Coeffs!")
print("=" * 105)

# Evaluate at the single physical lag tau = 0.92 d (can also optimize tau as 1 scalar!)
tau_opt = model_a_res["kepler_hill"]["params"][0]
L, Mp, F, D, Ls, Om = fundamental_angles(t_yr, tau_days=tau_opt)

# Elementary gravitational coordinates (ONLY 5 fundamental cosines/sines):
c2L, s2L = np.cos(2 * L), np.sin(2 * L)
cM,  sM  = np.cos(Mp),    np.sin(Mp)
cOm, sOm = np.cos(Om),    np.sin(Om)
c2D, s2D = np.cos(2 * D), np.sin(2 * D)
c2S, s2S = np.cos(2 * Ls), np.sin(2 * Ls)
cF,  sF  = np.cos(2 * F), np.sin(2 * F)

# Candidate low-degree polynomial monomials directly from expanding (1 + e cM + ...)^3 * (3 sin^2(delta) - 1):
# Note: cos(2L)*cos(M') splits into Mt (2L+M') and 2L-M'; using (c, s) pairs or direct gravitational products!
CANDS = [
    ("c2L [Mf: 13.66d]",                c2L),
    ("s2L [Mf quad]",                   s2L),
    ("cM [Mm: 27.55d]",                 cM),
    ("sM [Mm quad]",                    sM),
    ("c2L*cOm + s2L*sOm [Mf_nodal]",    c2L * cOm + s2L * sOm),   # = cos(2L - Om) = cos(L+F), 13.633d!
    ("c2L*cM - s2L*sM [Mt: 9.13d]",     c2L * cM - s2L * sM),     # = cos(2L + M'), 9.133d!
    ("c2L*cM + s2L*sM [2L-M': 27.09d]", c2L * cM + s2L * sM),     # = cos(2L - M')
    ("c2D*cM + s2D*sM [MSm: 31.81d]",   c2D * cM + s2D * sM),     # = cos(2D - M'), evection!
    ("c2D [MSf: 14.77d]",               c2D),                     # = cos(2D), variation!
    ("c2S [Ssa: 182.6d]",               c2S),                     # = cos(2Ls), semiannual!
    ("s2S [Ssa quad]",                  s2S),
    ("cOm [18.61y Nodal]",              cOm),
    ("cM^2 [2Mm: 13.78d]",              2 * cM * cM - 1.0),
    ("(c2L*cM-s2L*sM)*cOm [Mt_nodal]",  (c2L * cM - s2L * sM) * cOm + (s2L * cM + c2L * sM) * sOm),
    ("(c2L*c2D-s2L*s2D)*cM [Mtm: 9.56d]", (c2L * c2D - s2L * s2D) * cM + (s2L * c2D + c2L * s2D) * sM),
    ("c2F [2D draconic: 13.61d]",       cF),
    ("s2L*cOm - c2L*sOm [Mf_nod quad]", s2L * cOm - c2L * sOm),
    ("s2L*cM + c2L*sM [Mt quad]",       s2L * cM + c2L * sM),
    ("c2L*(c2M) [4T/2M: 6.85d]",        c2L * (2 * cM * cM - 1) - s2L * (2 * sM * cM)),
    ("c2L*c2D+s2L*s2D [2L-2D: 182.6]",  c2L * c2D + s2L * s2D),
]

# Greedy forward selection to show exact parsimony curve (K = 1 .. 16 terms)
selected = []
remaining = list(range(len(CANDS)))
poly_curve = []
for step in range(1, 17):
    best_cc, best_j, best_b = -1.0, None, None
    for j in remaining:
        cols = [np.ones_like(t_yr)] + [CANDS[i][1] for i in selected + [j]]
        X = np.column_stack(cols)
        b, *_ = np.linalg.lstsq(X, y_dlod, rcond=None)
        c = corr(X @ b, y_dlod)
        if c > best_cc:
            best_cc, best_j, best_b = c, j, b
    selected.append(best_j)
    remaining.remove(best_j)
    poly_curve.append((step, CANDS[best_j][0], best_cc, best_b))
    print(f"  K={step:2d} terms | + {CANDS[best_j][0]:36s} -> dLOD CC = {best_cc:.5f}")

# ── PURE MONOMIAL POLYNOMIAL TEST (strictly products of 4 real coordinates x_k without trig identities!) ──
print("\n" + "=" * 105)
print("MODEL C: Pure Degree-2 Multivariate Polynomial in 6 Real Orbital Coordinates (u_L, v_L, u_M, v_M, u_Om, u_2D)")
print("=" * 105)
# Just 3 complex phasors: z_L = e^{i2L}, z_M = e^{iM'}, z_Om = e^{iOm}, z_2D = e^{i2D}
# Test compact polynomial: P(t) = a1*c2L + a2*cM + a3*c2L*cOm + a4*s2L*sOm + a5*c2L*cM + a6*s2L*sM + a7*c2D*cM + a8*s2D*sM + a9*c2D + a10*c2S
X_pure = np.column_stack([
    np.ones_like(t_yr),
    c2L, s2L,
    cM, sM,
    c2L * cOm, s2L * sOm,
    c2L * cM, s2L * sM,
    c2D * cM, s2D * sM,
    c2D, c2S, s2S,
])
b_pure, *_ = np.linalg.lstsq(X_pure, y_dlod, rcond=None)
cc_pure = corr(X_pure @ b_pure, y_dlod)
print(f"13-term pure degree-2 polynomial in (c2L,s2L,cM,sM,cOm,sOm,c2D,s2D,c2S,s2S): dLOD CC = {cc_pure:.5f}")

# ── CLIMATE TRANSFER TEST: How well does our compact polynomial explain R5_atl4_cold_H10_l03's fitted tidal sum? ──
print("\n" + "=" * 105)
print("CLIMATE TRANSFER TEST: Compact Polynomial vs R5_atl4_cold_H10_l03 (NAO+Baltic+AMO+TNA) Fitted Tidal Sum")
print("=" * 105)
r5_cmp = np.loadtxt(R5_DIR / "dlod_compare.csv", delimiter=",")
t_r5, y_r5_ref, y_r5_fit = r5_cmp[:, 0], r5_cmp[:, 1], r5_cmp[:, 2]
L5, Mp5, F5, D5, Ls5, Om5 = fundamental_angles(t_r5, tau_days=tau_opt)
# Build top-12 gravitational polynomial basis on t_r5
def build_poly12(t_arr, tau):
    L_, Mp_, F_, D_, Ls_, Om_ = fundamental_angles(t_arr, tau_days=tau)
    c2L_, s2L_ = np.cos(2 * L_), np.sin(2 * L_)
    cM_,  sM_  = np.cos(Mp_),    np.sin(Mp_)
    cOm_, sOm_ = np.cos(Om_),    np.sin(Om_)
    c2D_, s2D_ = np.cos(2 * D_), np.sin(2 * D_)
    c2S_, s2S_ = np.cos(2 * Ls_), np.sin(2 * Ls_)
    c2F_, s2F_ = np.cos(2 * F_), np.sin(2 * F_)
    return np.column_stack([
        np.ones_like(t_arr),
        c2L_, s2L_,
        cM_, sM_,
        c2L_ * cOm_ + s2L_ * sOm_,   # Mf nodal
        c2L_ * cM_ - s2L_ * sM_,     # Mt
        c2D_ * cM_ + s2D_ * sM_,     # MSm evection
        c2D_,                         # MSf variation
        c2S_, s2S_,                   # Ssa (both phases — captures 90 deg Atlantic shift!)
        c2F_, s2F_,                   # 2D half-draconic
        (c2L_ * cM_ - s2L_ * sM_) * cOm_ + (s2L_ * cM_ + c2L_ * sM_) * sOm_,  # Mt nodal
    ])

X_r5 = build_poly12(t_r5, tau_opt)
b_r5_ref, *_ = np.linalg.lstsq(X_r5, y_r5_ref, rcond=None)
b_r5_fit, *_ = np.linalg.lstsq(X_r5, y_r5_fit, rcond=None)
print(f"13-coeff compact gravitational polynomial -> dLOD ref: CC = {corr(X_r5 @ b_r5_ref, y_r5_ref):.5f}")
print(f"13-coeff compact gravitational polynomial -> R5 fit:   CC = {corr(X_r5 @ b_r5_fit, y_r5_fit):.5f}")

# ── PLOT PARSIMONY CURVE ─────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2))
ks = [x[0] for x in poly_curve]
ccs = [x[2] for x in poly_curve]
ax1.plot(ks, ccs, "o-", color="#1f77b4", lw=2.2, label="Model B: Linear gravitational polynomial (K coeffs)")
ax1.scatter([4, 6, 8], [model_a_res[m]["cc"] for m in ("kepler_only","kepler_evect","kepler_hill")],
            color="#d95f02", marker="s", s=85, zorder=5, label="Model A: Closed-form Kepler-Hill u³(3sin²δ−1) (4/6/8 params)")
ax1.axhline(0.9717, color="#2ca02c", ls="--", lw=1.3, label="52-line empirical lpap (104 params): CC = 0.9717")
ax1.set_xlabel("Number of Free Parameters (K)")
ax1.set_ylabel("Daily dLOD Correlation CC (1962–2019, N=20,829)")
ax1.set_title("Parsimony Curve: Compact Gravitational Basis vs 104-Param lpap")
ax1.set_ylim(0.84, 0.985)
ax1.grid(True, alpha=0.3)
ax1.legend(loc="lower right", fontsize=8.2)

# Right panel: 1-year daily overlay (2000.0-2001.0) of measured dlod3.dat vs 8-param Model A vs 10-term Model B
sel = (t_yr >= 2000.0) & (t_yr <= 2001.0)
y_a8 = eval_model_a(model_a_res["kepler_hill"]["params"], t_yr, "kepler_hill")
# Scale-match y_a8 to y_dlod
b_sc = np.polyfit(y_a8, y_dlod, 1)
y_a8_s = np.polyval(b_sc, y_a8)
X10 = np.column_stack([np.ones_like(t_yr)] + [CANDS[i][1] for i in selected[:10]])
y_b10 = X10 @ poly_curve[9][3]
ax2.plot(t_yr[sel], y_dlod[sel], color="k", lw=1.1, alpha=0.8, label="Measured dLOD (dlod3.dat)")
ax2.plot(t_yr[sel], y_a8_s[sel], color="#d95f02", lw=1.2, label=f"Model A: 8-param Kepler-Hill (CC={model_a_res['kepler_hill']['cc']:.4f})")
ax2.plot(t_yr[sel], y_b10[sel], color="#1f77b4", lw=1.1, ls="--", label=f"Model B: 10-term Polynomial (CC={poly_curve[9][2]:.4f})")
ax2.set_xlabel("Year (daily samples, 2000)")
ax2.set_ylabel("dLOD (ms)")
ax2.set_title("2000 Daily dLOD: Zero-Meeus Compact Gravitational Models")
ax2.grid(True, alpha=0.3)
ax2.legend(loc="upper right", fontsize=8)

fig.suptitle("Compact Meeus-Free Gravitational Representation of dLOD (Replacing 52-Line / 104-Param lpap)", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT_DIR / "compact_tidal_parsimony.png", dpi=120)
plt.close(fig)
print("Saved", OUT_DIR / "compact_tidal_parsimony.png")
