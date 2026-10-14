#!/usr/bin/env python3
"""The legs.png figure (5-day panel and monthly panel, model against data, standardized) for a common-response run
(LEG_COMMON): five legs as in soi_stationarity/common_real. Each stretch is drawn from the leg whose interval it is:
1881-1936 leg 2, 1936-1992 leg 3, 1992-2005 leg 1 (leader), 2005-2015 leg 4, 2015-2027 leg 5. In every stretch outside
1992-2005 the harmonics above 30 carry the amplitudes and phases regressed on 1992-2005; only the slow terms are fitted
there. The monthly panel shows calendar-month means of the 5-day-grid model against the monthly data.
usage: legs_common_plot.py RUN_DIR OUT.png [title note]"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d, out = sys.argv[1], sys.argv[2]; note = sys.argv[3] if len(sys.argv) > 3 else ""; cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v, r: (v - r.mean()) / r.std()
L = {i: np.loadtxt(f"{d}/lte_results.{i}.csv", delimiter=",") for i in range(1, 6)}; SPAN = [(2, 1881, 1936), (3, 1936, 1992), (1, 1992, 2005), (4, 2005, 2015), (5, 2015, 2027)]
def stitch(legs):
    rows = [L[i][(L[i][:, 0] >= a0) & (L[i][:, 0] < a1)] for i, a0, a1 in SPAN if i in legs]; return np.vstack(rows)
five = stitch((1, 4, 5)); t5, m5, y5 = five[:, 0], five[:, 1], five[:, 2]; allr = stitch((1, 2, 3, 4, 5)); mi = np.floor((allr[:, 0] + 2.5 / 365.25 - 1880) * 12 + 1e-9).astype(int); u = np.unique(mi); tm = 1880 + u / 12.0; mm = np.array([allr[mi == k, 1].mean() for k in u]); ym = np.array([allr[mi == k, 2].mean() for k in u])
INK, MUTED, SURF, GRID, BL, OR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"; fig, ax = plt.subplots(2, 1, figsize=(17, 8.6), facecolor=SURF, gridspec_kw=dict(hspace=0.3))
for A, t, m, y, nm in ((ax[0], t5, m5, y5, "5-day SOI"), (ax[1], tm, mm, ym, "monthly means")):
    w = (t >= 1992) & (t < 2005); A.set_facecolor(SURF); A.axvspan(1992, 2005, color=GRID, alpha=0.7, lw=0); A.plot(t, z(y, y), color=BL, lw=1.0, label="data (standardized)"); A.plot(t, z(m, m), color=OR, lw=1.0, label="model (standardized)")
    A.set_title(f"{nm}: inside 1992-2005 {cc(m[w], y[w]):+.3f}, outside {cc(m[~w], y[~w]):+.3f}" + (f", after 2005 {cc(m[t >= 2005], y[t >= 2005]):+.3f}, before 1992 {cc(m[t < 1992], y[t < 1992]):+.3f}" if (t < 1992).sum() > 50 else f", after 2005 {cc(m[t >= 2005], y[t >= 2005]):+.3f}"), loc="left", fontsize=11, color=INK)
    A.legend(frameon=False, fontsize=9, ncol=2, loc="upper right"); A.set_xlim(t[0], t[-1]); A.tick_params(colors=MUTED)
    for a0 in (1936, 2015): A.axvline(a0, color=MUTED, lw=0.6, ls=":") if t[0] < a0 < t[-1] else None
    for s in ("top", "right"): A.spines[s].set_visible(False)
    for s in ("left", "bottom"): A.spines[s].set_color(GRID)
fig.suptitle("Common response: harmonics above 30 fixed from 1992-2005 (shaded), slow terms fitted per stretch (dotted lines)" + (f"  --  {note}" if note else ""), x=0.125, ha="left", fontsize=11.5, color=INK, y=0.955)
fig.savefig(out, dpi=100, bbox_inches="tight", facecolor=SURF); print("plot:", out)
print(f"5-day: inside {cc(m5[(t5 >= 1992) & (t5 < 2005)], y5[(t5 >= 1992) & (t5 < 2005)]):+.3f}, after 2005 {cc(m5[t5 >= 2005], y5[t5 >= 2005]):+.3f} | monthly means: before 1992 {cc(mm[tm < 1992], ym[tm < 1992]):+.3f}, inside {cc(mm[(tm >= 1992) & (tm < 2005)], ym[(tm >= 1992) & (tm < 2005)]):+.3f}, after {cc(mm[tm >= 2005], ym[tm >= 2005]):+.3f}")
