#!/usr/bin/env python3
"""Do the two legs use the same response? Compares the regressed amplitude and phase of every winding between
lt.exe.windings.1.json and lt.exe.windings.2.json (shared windings, per-leg regression).
Agreement = normalized complex correlation of the two amplitude vectors (1 = same response up to scale, 0 = unrelated,
-1 = opposite). Each leg's amplitudes are first divided by its data's standard deviation. usage: leg_agreement.py [dir]"""
import json, sys, numpy as np
d = sys.argv[1] if len(sys.argv) > 1 else "."; W = [json.load(open(f"{d}/lt.exe.windings.{i}.json")) for i in (1, 2)]; base = min((abs(k) for k, _, _ in W[0]["k_amp_phase"][1:2]), default=1.0)
sd = [np.loadtxt(f"{d}/lte_results.{i}.csv", delimiter=",")[:, 2].std() for i in (1, 2)]
Z = [np.array([a * np.exp(1j * p) for _, a, p in w["k_amp_phase"]]) / s for w, s in zip(W, sd)]; K = np.array([k for k, _, _ in W[0]["k_amp_phase"]]); o = np.argsort(np.abs(K))
agree = lambda m: float(np.real(np.vdot(Z[0][m], Z[1][m])) / np.sqrt(np.vdot(Z[0][m], Z[0][m]).real * np.vdot(Z[1][m], Z[1][m]).real))
print(f"{'winding':>9s} {'harmonic':>8s} | {'leg 1 amp':>9s} {'phase':>6s} | {'leg 2 amp':>9s} {'phase':>6s} | {'phase difference':>16s}")
for i in o:
    dp = np.degrees(np.angle(Z[1][i] / Z[0][i])); h = K[i] / K[1]; print(f"{K[i]:9.4f} {(format(round(h)) if abs(h - round(h)) < 1e-6 and round(h) >= 1 else 'low'):>8s} | {abs(Z[0][i]):9.3f} {np.degrees(np.angle(Z[0][i])):+6.0f} | {abs(Z[1][i]):9.3f} {np.degrees(np.angle(Z[1][i])):+6.0f} | {dp:+13.0f} deg")
all_ = np.ones(len(K), bool); lo = np.abs(K) <= 30 * abs(K[1]); print(f"\nagreement of the two legs' responses: all windings {agree(all_):+.2f}; harmonics up to 30 {agree(lo):+.2f}; above 30 {agree(~lo):+.2f}")
print(f"level {W[0]['level'] / sd[0]:+.3f} / {W[1]['level'] / sd[1]:+.3f}; slope in manifold (k0) {W[0]['k0'] / sd[0]:+.4f} / {W[1]['k0'] / sd[1]:+.4f}")
