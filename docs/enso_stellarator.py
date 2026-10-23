#!/usr/bin/env python3
"""Stellarator-like (winding-selective, wavenumber > 0) emulation of the equatorial Pacific (2026-10-22).

Field model on the common manifold F(t) (AMO/PDO/NINO4/IODE four-leg run, base k = 0.2017):
    zeta(x, t) = sum_h  [ sum_n B_hn Y_n(x) ] * (a_h cos 2 pi h k F + b_h sin 2 pi h k F)
with basin standing modes Y_n(x) = cos(n pi x), x = (lon - 120E) / 160 deg across the basin (n = 0, 1, 2).
Structures (fit 1950.0-2005.9, score 2006.0-2022.5):
  A  shared pattern     one spatial pattern P(x) for every winding (rank 1: the axisymmetric analogue)
  B1 single mode        each winding couples to ONE basin mode n_h (chosen on training data)
  B  stellarator        each winding has its own mix of n = 0, 1, 2 (3 spatial dof per winding)
  C  free               each winding has a free profile at every location (ceiling)
Windings: h = 1..9 x base and the slow winding (0.581 x base). Data: 8 equatorial Kaplan quads 130E-90W
(Oct2026_dlod), quadratic detrend, standardized. Scores: held-out correlation per location, and pattern correlation
of the held-out field at each time (mean over times).
    python3 enso_stellarator.py  -> enso_stellarator.png, stdout
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EX = HERE.parent / "experiments"
RUN = EX / "Oct2026_generator/legs/amo_pdo_nino4_iode"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
QUADS = [("kN000_E130", 130), ("kN000_E150", 150), ("kN000_E170", 170), ("kN000_W170", 190), ("kN000_W150", 210),
         ("kN000_W130", 230), ("kN000_W110", 250), ("kN000_W090", 270)]
p = json.load(open(RUN / "lt.exe.p")); KB = abs(p["ltep"][-1]); KS = abs(p["ltep"][0])
man = np.loadtxt(RUN / "lte_results.1.csv", delimiter=","); TM, FM = man[:, 0], man[:, 3]
t = np.arange(1950.0, 2022.5, 1 / 12) + 1 / 24
Y = []
for q, _ in QUADS:
    a = np.loadtxt(EX / "Oct2026_dlod" / q / f"{q}.dat"); y = np.interp(t, a[:, 0], a[:, 1])
    y = y - np.polyval(np.polyfit(t, y, 2), t); Y.append(y / y.std())
Y = np.array(Y).T                                                     # time x location
x = (np.array([lon for _, lon in QUADS]) - 120.0) / 160.0
MODES = np.array([np.cos(n * np.pi * x) for n in range(3)]).T          # location x n
F = np.interp(t, TM, FM)
KS_ = [h * KB for h in range(1, 10)] + [KS]; NAMES = [f"{h}x" for h in range(1, 10)] + ["slow"]
T = np.stack([np.stack([np.cos(2 * np.pi * k * F), np.sin(2 * np.pi * k * F)], 1) for k in KS_], 1)   # time x H x 2
tr = t < 2006.0; te = ~tr
H, L = len(KS_), Y.shape[1]


def fit_free(Ytr, Ttr):
    X = Ttr.reshape(len(Ttr), -1); X1 = np.column_stack([X, np.ones(len(X))])
    C = np.linalg.lstsq(X1, Ytr, rcond=None)[0]; return C                               # (2H+1) x L


def predict_free(C, Tt):
    X = Tt.reshape(len(Tt), -1); return np.column_stack([X, np.ones(len(X))]) @ C


def fit_spatial(S_list, Ytr, Ttr):
    """each winding h: coefficient profile in span(S_list[h]) (location x d_h); solve jointly by least squares"""
    cols = []
    for h in range(H):
        for j in range(2):
            for d in range(S_list[h].shape[1]):
                cols.append((h, j, d))
    A = np.zeros((Ytr.size, len(cols) + L))
    for c, (h, j, d) in enumerate(cols):
        A[:, c] = np.outer(Ttr[:, h, j], S_list[h][:, d]).ravel()
    for l in range(L):                                                  # per-location mean
        e = np.zeros((len(Ttr), L)); e[:, l] = 1; A[:, len(cols) + l] = e.ravel()
    w = np.linalg.lstsq(A, Ytr.ravel(), rcond=None)[0]
    return cols, w


def predict_spatial(S_list, cols, w, Tt):
    out = np.zeros((len(Tt), L))
    for c, (h, j, d) in enumerate(cols):
        out += w[c] * np.outer(Tt[:, h, j], S_list[h][:, d])
    return out + w[len(cols):]


def scores(P_, Yt):
    loc = np.array([np.corrcoef(P_[:, l], Yt[:, l])[0, 1] for l in range(L)])
    pat = np.nanmean([np.corrcoef(P_[i] - P_[i].mean(), Yt[i] - Yt[i].mean())[0, 1] for i in range(len(Yt))])
    return loc, pat


RES = {}
# C free
C = fit_free(Y[tr], T[tr]); RES["C free (16 per winding)"] = (predict_free(C, T), H * 2 * L)
# A shared pattern: rank-1 of the free coefficients, then refit amplitudes with the pattern fixed
U, s, Vt = np.linalg.svd(C[:-1].T, full_matrices=False); Pshared = U[:, :1]
cols, w = fit_spatial([Pshared] * H, Y[tr], T[tr]); RES["A shared pattern (rank 1)"] = (predict_spatial([Pshared] * H, cols, w, T), L + 2 * H)
# B stellarator: n = 0, 1, 2 mix per winding
SB = [MODES] * H; cols, w = fit_spatial(SB, Y[tr], T[tr]); RES["B stellarator (n=0,1,2 per winding)"] = (predict_spatial(SB, cols, w, T), 6 * H)
# B1 single mode per winding, chosen on training fit
best = []
for h in range(H):
    errs = []
    for n in range(3):
        S = [MODES] * H; S = [S[i] if i != h else MODES[:, n:n + 1] for i in range(H)]
        cl, ww = fit_spatial(S, Y[tr], T[tr]); errs.append(np.mean((predict_spatial(S, cl, ww, T[tr]) - Y[tr]) ** 2))
    best.append(int(np.argmin(errs)))
S1 = [MODES[:, n:n + 1] for n in best]; cols, w = fit_spatial(S1, Y[tr], T[tr])
RES["B1 one basin mode per winding"] = (predict_spatial(S1, cols, w, T), 2 * H)

print(f"basin modes cos(n pi x), x = (lon - 120E)/160; nodes: n=1 at 200E, n=2 at 160E and 240E")
print(f"{'structure':40s} {'params':>6s} | held-out r by location (130E ... 90W)                      | mean r | pattern r")
for k, (P_, npar) in RES.items():
    loc, pat = scores(P_[te], Y[te]); loctr, pattr = scores(P_[tr], Y[tr])
    print(f"{k:40s} {npar:6d} | " + " ".join(f"{v:+.2f}" for v in loc) + f" | {loc.mean():+.3f} | {pat:+.3f}   (train mean r {loctr.mean():+.3f}, pattern {pattr:+.3f})")

# selectivity: energy share of n = 0, 1, 2 in each winding's B (stellarator) profile, from the free fit projected on modes
print("\nwinding -> basin-mode selectivity (share of each winding's free-fit profile in n = 0 / 1 / 2; residual = not in n <= 2)")
Q, _ = np.linalg.qr(MODES)
for h in range(H):
    prof = C[2 * h:2 * h + 2].T                                   # location x 2 (cos, sin)
    tot = (prof ** 2).sum(); proj = MODES @ np.linalg.lstsq(MODES, prof, rcond=None)[0]
    sh = [((np.outer(MODES[:, n], np.linalg.lstsq(MODES[:, n:n + 1], prof, rcond=None)[0][0])) ** 2).sum() / tot for n in range(3)]
    amp = np.sqrt(tot / L)
    print(f"   {NAMES[h]:5s} amp {amp:.3f}: n0 {sh[0]:.2f}  n1 {sh[1]:.2f}  n2 {sh[2]:.2f}   in span(n<=2) {(proj ** 2).sum() / tot:.2f}   chosen single mode n = {best[h]}")

fig, axs = plt.subplots(1, 2, figsize=(15, 4.8), facecolor=SURF)
ax = axs[0]; lon = [l for _, l in QUADS]
for (k, (P_, _)), c in zip(RES.items(), ("#52514e", "#eda100", "#2a78d6", "#1baf7a")):
    loc, _ = scores(P_[te], Y[te]); ax.plot(lon, loc, marker="o", color=c, lw=1.5, label=k)
ax.set_xticks(lon); ax.set_xticklabels(["130E", "150E", "170E", "170W", "150W", "130W", "110W", "90W"])
ax.set_ylabel("held-out r (2006-2022)", color=MUTED); ax.legend(frameon=False, fontsize=8.5)
ax.set_title("held-out skill along the equator", loc="left", fontsize=10.5, color=INK)
ax = axs[1]; xx = np.linspace(0, 1, 100)
for n, c in zip(range(3), ("#2a78d6", "#eb6834", "#1baf7a")):
    ax.plot(120 + 160 * xx, np.cos(n * np.pi * xx), color=c, lw=1.8, label=f"basin mode n = {n}")
ax.axvline(160, color=MUTED, ls=":", lw=1); ax.text(161, -0.95, "observed node ~160E", fontsize=8.5, color=MUTED)
ax.set_title("the spatial modes the windings can select", loc="left", fontsize=10.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
for a in axs:
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "enso_stellarator.png", dpi=95, facecolor=SURF)
print("plot: enso_stellarator.png")
