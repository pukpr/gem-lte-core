#!/usr/bin/env python3
"""How the high harmonics act in each leg: each leg's model split into the harmonics up to 30 and above 30 (same
windings, each leg's own amplitudes), the part of the data each explains, train and held out 2000-2005, and the
agreement of the two legs' amplitudes and phases (sign-free: the 5-day series is the INVERTED SOI)."""
import json, numpy as np
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); p = json.load(open("lt.exe.p")); kb = abs(p["ltep"][-1]); W = [json.load(open(f"lt.exe.windings.{i}.json"))["k_amp_phase"] for i in (1, 2)]
hi = [j for j, x in enumerate(W[0]) if j >= 2 and round(abs(x[0]) / kb) > 30]; lo = [j for j in range(len(W[0])) if j not in hi]
Z = [np.array([x[1] * np.exp(1j * x[2]) for x in w]) for w in W]; ag = lambda m: float(np.real(np.vdot(Z[0][m], Z[1][m])) / np.sqrt(np.vdot(Z[0][m], Z[0][m]).real * np.vdot(Z[1][m], Z[1][m]).real))
print(f"high harmonics: {sorted(round(abs(W[0][j][0]) / kb) for j in hi)}; leg agreement of amplitudes and phases (magnitude matters; sign flips with the inverted series): all {ag(list(range(len(W[0])))):+.2f}, up to 30 {ag(lo):+.2f}, above 30 {ag(hi):+.2f}")
for i, nm, N in ((1, "CRU SOI monthly", 12), (2, "SOI 5-day", 73)):
    a = np.loadtxt(f"lte_results.{i}.csv", delimiter=","); t, m, y, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]; g = (t >= 2000) & (t < 2005); w = W[i - 1]
    part = lambda js: sum(w[j][1] * np.sin(2 * np.pi * w[j][0] * F + w[j][2]) for j in js); H = part(hi); Lw = part(lo)
    r = y - np.polyval(np.polyfit(Lw[~g], y[~g], 1), Lw); print(f"{nm:16s} model {cc(m[~g], y[~g]):+.3f} train / {cc(m[g], y[g]):+.3f} held out | high-harmonic part (pre-lag) against data less the low part: train {cc(H[~g], r[~g]):+.3f}, held out {cc(H[g], r[g]):+.3f}; its share of the model's winding variance {100 * H.var() / (H + Lw).var():.0f}%")
