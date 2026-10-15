#!/usr/bin/env python3
"""Darwin and Tahiti scores with the mean seasonal cycle taken out of both model and data (it dominates the raw
pressure and would otherwise inflate the correlations), plus the agreement of the two legs' responses."""
import json, numpy as np
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); Z = []
for i, nm in ((1, "Darwin"), (2, "Tahiti")):
    a = np.loadtxt(f"lte_results.{i}.csv", delimiter=","); t, m, y = a[:, 0], a[:, 1], a[:, 2]; g = (t >= 2000) & (t < 2005); mo = np.round((t % 1) * 12).astype(int) % 12
    ds = lambda v: v - np.array([v[(mo == k) & ~g].mean() for k in range(12)])[mo]; dm, dy = ds(m), ds(y)
    print(f"{nm}: train {cc(m[~g], y[~g]):+.3f}, held out 2000-05 {cc(m[g], y[g]):+.3f}; before 1880 {cc(m[t < 1880], y[t < 1880]):+.3f}, 1880 on (outside 2000-05) {cc(m[(t >= 1880) & ~g], y[(t >= 1880) & ~g]):+.3f}")
    w = json.load(open(f"lt.exe.windings.{i}.json"))["k_amp_phase"]; Z.append(np.array([x[1] * np.exp(1j * x[2]) for x in w]))
print("response agreement Darwin vs Tahiti %+.2f (opposite sign expected for a see-saw)" % float(np.real(np.vdot(Z[0], Z[1])) / np.sqrt(np.vdot(Z[0], Z[0]).real * np.vdot(Z[1], Z[1]).real)))
a1, a2 = np.loadtxt("lte_results.1.csv", delimiter=","), np.loadtxt("lte_results.2.csv", delimiter=","); k1 = {round(t * 12): f for t, f in a1[:, [0, 3]]}; c = np.array([(k1[round(t * 12)], f) for t, f in a2[:, [0, 3]] if round(t * 12) in k1]); print(f"the two legs' manifolds on common months CC {cc(c[:, 0], c[:, 1]):+.5f}")
