#!/usr/bin/env python3
"""With the integration now starting in 1860, find the starting value (init at 1860) whose manifold from 1880 on best
matches the seed's own manifold (integrated from 1880), so the seed keeps its calibrated meaning on the overlap."""
import json, os, subprocess, numpy as np
ref = np.loadtxt("first_attempt_mean_only/lte_results.1.csv", delimiter=","); kr = {round(t * 12): f for t, f in ref[:, [0, 3]]}; P = json.load(open("seed.p"))
def score(v):
    p = dict(P); p["init"] = float(v); json.dump(p, open("lt.exe.p", "w"), indent=2); subprocess.run("TEST_ONLY=TRUE ./lte_run.sh 60 2 > scan.out 2>&1", shell=True, stdin=subprocess.DEVNULL)
    a = np.loadtxt("lte_results.1.csv", delimiter=","); c = np.array([(f, kr[round(t * 12)]) for t, f in a[:, [0, 3]] if round(t * 12) in kr]); return np.sqrt(np.mean((c[:, 0] - c[:, 1]) ** 2)), np.corrcoef(c[:, 0], c[:, 1])[0, 1]
grid = np.linspace(-8, 8, 17); R = [(v, *score(v)) for v in grid]; best = min(R, key=lambda r: r[1]); print("coarse:", " ".join(f"{v:+.0f}:{r:.2f}" for v, r, c in R))
for step in (0.5, 0.1, 0.02):
    vs = best[0] + step * np.arange(-4, 5); R = [(v, *score(v)) for v in vs]; best = min(R + [best], key=lambda r: r[1])
print(f"best init at 1860: {best[0]:+.3f}; manifold from 1880 on against the seed's: rms difference {best[1]:.4f}, CC {best[2]:.5f}"); p = dict(P); p["init"] = float(best[0]); json.dump(p, open("seed_1860.p", "w"), indent=2); json.dump(p, open("lt.exe.p", "w"), indent=2)
