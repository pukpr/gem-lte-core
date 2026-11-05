"""Cross-terms of the LTE modulation in Fourier space (Jacobi-Anger expansion).

The LTE response is y(t) = k0*F(t) + sum_w A_w sin(2 pi k_w F(t) + phi_w), with F the forcing manifold.
If F(t) = c + sum_j a_j cos(2 pi f_j t + th_j) (a few aliased tidal lines), then exactly

  exp(i 2 pi k F) = exp(i 2 pi k c) * prod_j sum_n i^n J_n(2 pi k a_j) exp(i n (2 pi f_j t + th_j))

so the response is a set of lines at every integer combination sum_j n_j f_j, with amplitude
prod_j J_{n_j}(2 pi k a_j). The original lines survive as the n = (0..,+-1,..0) terms ("vestiges"),
weighted by J_1 * prod J_0; everything else is cross-terms. This script

  1. extracts the line content of the ENSO benchmark manifold F (lte_results.csv col 4) and labels
     each line with the nearest annual-impulse alias of a tidal period in lt.exe.p,
  2. expands the slow windings analytically to a given order and checks the line sum against the
     directly computed response,
  3. reports which cross-terms build each interannual peak, and compares the predicted spectrum
     with the direct response and the NINO4 data.

  python3 lte_crossterms.py [--run DIR] [--lines 14] [--order 4]
"""
import argparse, itertools, json
import numpy as np
from scipy.special import jv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import peak_estimate as pe

YEAR = 365.2422


def extract_lines(t, x, n, fmax=1.0, fmin=0.03):
    """Greedy sinusoid extraction with joint least-squares re-fit at each step."""
    T = t[-1] - t[0]
    grid = np.arange(fmin, fmax, 1 / (16 * T))
    fs, r = [], x - x.mean()
    for _ in range(n):
        C, S = np.cos(2 * np.pi * np.outer(t, grid)), np.sin(2 * np.pi * np.outer(t, grid))
        p = (r @ C) ** 2 + (r @ S) ** 2
        f0 = grid[np.argmax(p)]
        fine = np.linspace(f0 - 1 / (8 * T), f0 + 1 / (8 * T), 401)
        pf = [(r @ np.cos(2 * np.pi * f * t)) ** 2 + (r @ np.sin(2 * np.pi * f * t)) ** 2 for f in fine]
        fs.append(fine[int(np.argmax(pf))])
        B = np.column_stack([np.ones_like(t)] + [g(2 * np.pi * f * t) for f in fs for g in (np.cos, np.sin)])
        coef = np.linalg.lstsq(B, x, rcond=None)[0]
        r = x - B @ coef
    c = coef[0]
    ab = coef[1:].reshape(-1, 2)
    a = np.hypot(ab[:, 0], ab[:, 1]); th = np.arctan2(-ab[:, 1], ab[:, 0])   # a cos(wt + th)
    return c, np.array(fs), a, th, 1 - r.var() / x.var()


TIDES = {"Mf": 13.660791, "Mf'": 13.633416, "Mm": 27.554550, "Trop": 27.321582, "Drac": 27.212221,
         "Drac/2": 13.606111, "Mt": 9.132931, "Msf": 14.765294, "Trop/3": 9.107194}
BASE = {k: (YEAR / P) % 1 for k, P in TIDES.items()}     # fractional cpy: what an annual impulse keeps
BASE["18.6y"] = YEAR / 6798.38


def alias_labels(fs, lpap=None, tol=0.004):
    """Label each frequency as the simplest n1*b1 + n2*b2 + m (m integer: annual comb) of the alias bases."""
    keys = list(BASE)
    cands = []
    for k in keys:
        for n in (1, 2):
            for m in range(-3, 4):
                cands.append((abs(n * BASE[k] + m), (f"{n if n > 1 else ''}{k}" + (f"{m:+d}" if m else "")), n))
    for i, k1 in enumerate(keys):
        for k2 in keys[i + 1:]:
            for s in (1, -1):
                for m in range(-3, 4):
                    cands.append((abs(BASE[k1] + s * BASE[k2] + m), f"{k1}{'+' if s > 0 else '-'}{k2}" + (f"{m:+d}" if m else ""), 2))
    out = []
    for f in fs:
        ok = [c for c in cands if abs(c[0] - f) < tol]
        if ok:
            g, lab, _ = min(ok, key=lambda c: (c[2], abs(c[0] - f))); out.append((lab, g))
        else:
            out.append(("unlabelled", np.nan))
    return out


def expand(c, fs, a, th, k, A, phi, order):
    """Lines of A sin(2 pi k F + phi): returns {freq: complex amplitude} (y = sum Re(Z e^{i 2pi f t}))."""
    m = len(fs); lines = {}
    z = 2 * np.pi * k * a
    base = A * np.exp(1j * (2 * np.pi * k * c + phi))
    # orders: vectors n with sum|n_j| <= order
    def vecs():
        for tot in range(order + 1):
            for idx in itertools.combinations_with_replacement(range(m), tot):
                cnt = np.bincount(idx, minlength=m)
                nz = np.nonzero(cnt)[0]
                for signs in itertools.product((1, -1), repeat=len(nz)):
                    n = np.zeros(m, int); n[nz] = cnt[nz] * np.array(signs); yield n
    for n in vecs():
        amp = np.prod([jv(nj, zj) * (1j ** nj) * np.exp(1j * nj * tj) for nj, zj, tj in zip(n, z, th)])
        f = float(n @ fs)
        # sin(X) = Im(e^{iX}) = (e^{iX} - e^{-iX})/(2i): line at f with -i/2 * base*amp and its conjugate at -f
        Zf = -0.5j * base * amp
        key = round(f, 6)
        lines[key] = lines.get(key, 0) + Zf
        lines[-key] = lines.get(-key, 0) + np.conj(Zf)
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="../../experiments/Oct2026_generator/legs/nino4_slow_then_harm/benchmark_enso_2026-11-01")
    ap.add_argument("--lines", type=int, default=30)
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--out", default="lte_crossterms.png")
    a = ap.parse_args()
    d = np.genfromtxt(f"{a.run}/lte_results.csv", delimiter=",")
    t, data, F = d[:, 0], d[:, 2], d[:, 3]
    w = json.load(open(f"{a.run}/lt.exe.windings.json")); p = json.load(open(f"{a.run}/lt.exe.p"))
    k0 = w["k0"]; wind = w["k_amp_phase"]
    base_k = min(kk for kk, _, _ in wind)
    slow = [x for x in wind if x[0] < 3 * base_k]      # 1x and 2x: expanded analytically
    fast = [x for x in wind if x[0] >= 3 * base_k]     # 15x: modulation index ~19, left numerical

    c, fs, amps, th, r2 = extract_lines(t, F, a.lines)
    labs = alias_labels(fs, p["lpap"])
    print(f"manifold F: {a.lines} lines explain {100 * r2:.1f}% of its variance")
    order = np.argsort(-amps)
    for i in order:
        lab, g = labs[i]
        print(f"  {fs[i]:.4f} cpy  {1 / fs[i]:6.2f} yr  a={amps[i]:6.2f}  index 2pi k a = {2 * np.pi * base_k * amps[i]:.2f}"
              f"  <- {lab}" + (f" ({1 / g:.3f} yr)" if g == g else ""))

    # direct responses
    Fl = c + sum(ai * np.cos(2 * np.pi * fi * t + ti) for fi, ai, ti in zip(fs, amps, th))
    y_slow = k0 * F + sum(A * np.sin(2 * np.pi * kk * F + ph) for kk, A, ph in slow)
    y_fast = sum(A * np.sin(2 * np.pi * kk * F + ph) for kk, A, ph in fast)
    # analytic line sum (slow windings, on the line model of F)
    lines = {}
    for i in range(len(fs)):
        for s in (1, -1):
            key = round(s * fs[i], 6); lines[key] = lines.get(key, 0) + 0.5 * k0 * amps[i] * np.exp(1j * s * th[i])
    for kk, A, ph in slow:
        for f, Z in expand(c, fs, amps, th, kk, A, ph, a.order).items():
            lines[f] = lines.get(f, 0) + Z
    y_ana = sum((Z * np.exp(2j * np.pi * f * t)).real for f, Z in lines.items())
    y_slow_l = k0 * Fl + sum(A * np.sin(2 * np.pi * kk * Fl + ph) for kk, A, ph in slow)
    cc = lambda u, v: np.corrcoef(u - u.mean(), v - v.mean())[0, 1]
    print(f"\nanalytic expansion to order {a.order}: {len(lines) // 2} lines; r with direct slow response on the line-model F = {cc(y_ana, y_slow_l):.4f};"
          f" on the true F = {cc(y_ana, y_slow):.3f}")
    print(f"fast (15x) part: {100 * y_fast.var() / (y_slow + y_fast).var():.0f}% of response variance, left numerical")

    # which cross-terms build each band peak
    pos = sorted(((f, abs(Z)) for f, Z in lines.items() if 1 / 8 <= f <= 1 / 2), key=lambda x: -x[1])
    def describe(f):
        best = None
        for n in itertools.product(range(-3, 4), repeat=min(4, len(fs))):
            pass
        return ""
    # brute-force label: find the n-vector (order <= a.order) producing each strong line, over the top 6 F lines
    top = order[:8]
    combos = {}
    for tot in range(1, a.order + 1):
        for idx in itertools.combinations_with_replacement(top, tot):
            cnt = {}
            for j in idx: cnt[j] = cnt.get(j, 0) + 1
            for signs in itertools.product((1, -1), repeat=len(cnt)):
                n = {j: s * cnt[j] for j, s in zip(cnt, signs)}
                f = round(sum(nj * fs[j] for j, nj in n.items()), 6)
                if f > 0 and f not in combos: combos[f] = n
    fmt_n = lambda n: " ".join(f"{'+' if v > 0 else '-'}{abs(v) if abs(v) > 1 else ''}[{1 / fs[j]:.2f}]" for j, v in sorted(n.items()))
    print("\nstrongest response lines in 2-8 yr (amplitude, period, combination of manifold lines [yr]):")
    for f, A in pos[:16]:
        n = combos.get(round(f, 6)); lab = alias_labels([f])[0][0]
        print(f"  {A:.3f}  {1 / f:5.2f} yr  {lab:14s} {fmt_n(n) if n else ''}")

    # spectra: predicted lines (smoothed) vs direct response vs data
    fig, ax = plt.subplots(figsize=(11, 4.6))
    fgrid = np.linspace(1 / 12, 1, 4000)
    T = t[-1] - t[0]; sig = 1 / T
    S = np.zeros_like(fgrid)
    for f, Z in lines.items():
        if f > 0: S += abs(Z) ** 2 * np.exp(-0.5 * ((fgrid - f) / sig) ** 2)
    for lab, x, col, lw in [("direct response, slow windings", y_slow, "magenta", 0.5), # 1.6
                            ("direct response, all windings", y_slow + y_fast, "red", 1.4),
                            ("NINO4 data", data, "blue", 1.2)]:
        fm, Pk = pe.multitaper(pe.prep(x)); P = Pk.mean(0); m = (fm > 1 / 12) & (fm < 1)
        ax.plot(1 / fm[m], P[m] / P[m].max(), color=col, lw=lw, label=lab)
    m = fgrid > 1 / 12
    ax.plot(1 / fgrid[m], S[m] / S[m].max(), color="#0a7c84", lw=0.2, ls="--", label=f"Jacobi-Anger lines, order {a.order} (smoothed)")
    for i in order[:6]:
        ax.axvline(1 / fs[i], color="green", lw=0.8, alpha=0.6)
    ax.plot([], [], color="green", lw=0.8, label="manifold lines (top 6)")
    ax.set_xlim(1.2, 12); ax.set_xlabel("period (yr)"); ax.set_ylabel("normalised power"); ax.grid(alpha=.3)
    ax.legend(fontsize=8); ax.set_title("LTE modulation of the ENSO benchmark manifold: vestiges and cross-terms")
    fig.tight_layout(); fig.savefig(a.out, dpi=130); print(f"\nplot: {a.out}")


if __name__ == "__main__":
    main()
