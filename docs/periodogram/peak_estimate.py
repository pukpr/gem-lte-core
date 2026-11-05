"""High-resolution estimates of the interannual spectral peaks of a Nino index.

pg.py uses Welch with nperseg=256 months, so its frequency bins are 1/21.3 cpy apart and
its "peaks" fall on bin centres (64/18 = 3.56 yr, 64/12 = 5.33 yr). This script locates the
peaks without that binning, by several methods, and compares them with tidal aliases:

  periodogram  full-record, Hann taper, zero-padded 16x
  multitaper   DPSS NW=3 (5 tapers), zero-padded 16x; jackknife spread over tapers
  burg         maximum-entropy (Burg AR) spectrum, over a range of AR orders
  ssa-esprit   singular spectrum analysis: frequencies of the leading oscillatory
               eigen-pairs via ESPRIT (rotational invariance of the signal subspace)

  python3 peak_estimate.py                      # both records, plot peak_estimate.png
  python3 peak_estimate.py --cpc file.txt       # CPC ersst5.nino.mth.91-20.ascii (1950-)
"""
import argparse
import numpy as np
from scipy.signal import windows, find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

YEAR = 365.2422
ALIASES = {  # annual-impulse aliases (yr): |f - round(f)| with f = YEAR / P
    "Mf 13.661 d": 13.660791, "Mm 27.555 d": 27.554550, "Drac/3 9.071 d": 9.070740,
}
alias_yr = {k: 1 / abs(YEAR / P - round(YEAR / P)) for k, P in ALIASES.items()}
BAND = (1 / 8.0, 1 / 2.0)  # cpy


def load_dat(path, t0=None, t1=None):
    d = np.loadtxt(path)
    t, x = d[:, 0], d[:, 1]
    m = np.ones_like(t, bool)
    if t0: m &= t >= t0
    if t1: m &= t < t1
    return t[m], x[m]


def load_cpc(path):
    rows = [l.split() for l in open(path) if l[:1].isdigit()]
    t = np.array([int(r[0]) + (int(r[1]) - 1) / 12 for r in rows])
    x = np.array([float(r[7]) for r in rows])  # NINO4 ANOM
    return t, x


def prep(x):
    t = np.arange(len(x))
    return x - np.polyval(np.polyfit(t, x, 2), t)  # quadratic detrend


def band_peaks(f, p, n=2):
    m = (f >= BAND[0]) & (f <= BAND[1])
    fb, pb = f[m], p[m]
    idx, _ = find_peaks(pb)
    idx = idx[np.argsort(pb[idx])[::-1]][:n]
    out = []
    for i in idx:  # parabolic refinement on log power
        if 0 < i < len(pb) - 1:
            y0, y1, y2 = np.log(pb[i - 1:i + 2]); d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            out.append(fb[i] + d * (fb[1] - fb[0]))
        else:
            out.append(fb[i])
    return sorted(out)


def periodogram(x, pad=16):
    w = windows.hann(len(x)); N = pad * len(x)
    P = np.abs(np.fft.rfft(x * w, N)) ** 2
    return np.fft.rfftfreq(N, 1 / 12), P


def multitaper(x, NW=3, pad=16):
    K = 2 * NW - 1; tap = windows.dpss(len(x), NW, K); N = pad * len(x)
    Pk = np.array([np.abs(np.fft.rfft(x * w, N)) ** 2 for w in tap])
    return np.fft.rfftfreq(N, 1 / 12), Pk


def burg(x, order):
    x = np.asarray(x, float); N = len(x)
    a = np.array([1.0]); E = x @ x / N
    f, b = x[1:].copy(), x[:-1].copy()
    for _ in range(order):
        k = -2 * (f @ b) / (f @ f + b @ b)
        a = np.concatenate([a, [0]]); a = a + k * a[::-1]
        E *= 1 - k * k
        f, b = f[1:] + k * b[1:], b[:-1] + k * f[:-1]
    return a, E


def burg_psd(x, order, nf=20000):
    a, E = burg(x, order)
    f = np.linspace(0, 6, nf)
    z = np.exp(-2j * np.pi * np.outer(f / 12, np.arange(len(a))))
    return f, E / np.abs(z @ a) ** 2


def ssa_esprit(x, M=240, r=12):
    """Frequencies (cpy) of the leading r SSA components, via ESPRIT on the signal subspace."""
    N = len(x); K = N - M + 1
    H = np.array([x[i:i + M] for i in range(K)]).T
    U, s, _ = np.linalg.svd(H, full_matrices=False)
    Us = U[:, :r]
    Phi = np.linalg.lstsq(Us[:-1], Us[1:], rcond=None)[0]
    lam = np.linalg.eigvals(Phi)
    fr = np.abs(np.angle(lam)) / (2 * np.pi) * 12
    damp = np.abs(lam)
    keep = fr > 1e-6
    return np.unique(np.round(fr[keep], 6)), s[:r] ** 2 / (s ** 2).sum()


def analyse(name, t, x):
    x = prep(x); T = len(x) / 12
    print(f"\n== {name}: {t[0]:.1f}-{t[-1] + 1/12:.1f}, {T:.1f} yr, Rayleigh 1/T = {1/T:.4f} cpy "
          f"(+-{(1/T) * 3.8**2 / 2:.2f} yr at 3.8 yr)")
    res = {}
    f, P = periodogram(x); res["periodogram"] = (f, P, band_peaks(f, P))
    f, Pk = multitaper(x); Pm = Pk.mean(0)
    jk = []
    for k in range(len(Pk)):
        pk = band_peaks(f, np.delete(Pk, k, 0).mean(0), 1); jk += pk
    res["multitaper"] = (f, Pm, band_peaks(f, Pm))
    res["mt_jack"] = np.std(jk) * np.sqrt(len(Pk) - 1)
    orders = [24, 36, 48, 60, 84, 120]
    bp = {}
    for o in orders:
        fb, Pb = burg_psd(x, o); bp[o] = (fb, Pb, band_peaks(fb, Pb))
    res["burg"] = bp
    fr, var = ssa_esprit(x)
    res["ssa"] = sorted(f for f in fr if BAND[0] <= f <= BAND[1])
    fmt = lambda fs: ", ".join(f"{1/f:.2f} yr" for f in fs)
    print(f"  periodogram (Hann, 16x pad):  {fmt(res['periodogram'][2])}")
    print(f"  multitaper NW=3:              {fmt(res['multitaper'][2])}  (jackknife sd of top peak {res['mt_jack']:.4f} cpy)")
    for o in orders:
        print(f"  Burg MEM order {o:3d} mo:       {fmt(bp[o][2])}")
    print(f"  SSA (M=240 mo) + ESPRIT:      {fmt(res['ssa'])}")
    return res


def windows_scan(t, x, width=60, step=10):
    out = []
    for a in np.arange(t[0], t[-1] - width + 1e-6, step):
        m = (t >= a) & (t < a + width)
        f, Pk = multitaper(prep(x[m])); pk = band_peaks(f, Pk.mean(0), 2)
        out.append((a, pk))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", default="../../experiments/Feb2026/nino4/nino4.dat")
    ap.add_argument("--cpc", default=None, help="CPC ersst5.nino.mth.91-20.ascii (as pg.py downloads)")
    ap.add_argument("--out", default="peak_estimate.png")
    a = ap.parse_args()
    print("annual-impulse aliases: " + ", ".join(f"{k} -> {v:.3f} yr" for k, v in alias_yr.items()))
    sets = [("NINO4 ERSST 1880-2024", *load_dat(a.dat))]
    if a.cpc:
        sets.append(("NINO4 CPC ERSSTv5 1950-", *load_cpc(a.cpc)))
    allres = [(n, t, x, analyse(n, t, x)) for n, t, x in sets]

    n, t, x, _ = allres[0]
    print(f"\n== 60-yr windows, multitaper peaks ({n})")
    for a0, pk in windows_scan(t, x):
        print(f"  {a0:.0f}-{a0 + 60:.0f}: " + ", ".join(f"{1/f:.2f}" for f in pk))

    fig, axes = plt.subplots(len(allres), 1, figsize=(11, 4.2 * len(allres)), squeeze=False)
    for ax, (n, t, x, r) in zip(axes[:, 0], allres):
        f, P, _ = r["periodogram"]; m = (f > 1 / 12) & (f < 1)
        ax.plot(1 / f[m], P[m] / P[m].max(), color="#9aa5b8", lw=1, label="periodogram (Hann, padded)")
        f, P, _ = r["multitaper"]; m = (f > 1 / 12) & (f < 1)
        ax.plot(1 / f[m], P[m] / P[m].max(), color="#1b4f72", lw=2, label="multitaper NW=3")
        fb, Pb, _ = r["burg"][60]; m = (fb > 1 / 12) & (fb < 1)
        ax.plot(1 / fb[m], Pb[m] / Pb[m].max(), color="#b5640a", lw=1.6, label="Burg MEM, order 60 mo")
        for f in r["ssa"]:
            ax.axvline(1 / f, color="#0a7c84", lw=1, ls=":")
        ax.axvline(alias_yr["Mf 13.661 d"], color="#c0392b", lw=1.6, label=f"Mf annual alias {alias_yr['Mf 13.661 d']:.3f} yr")
        ax.plot([], [], color="#0a7c84", ls=":", label="SSA-ESPRIT frequencies")
        ax.set_xlim(1.2, 10); ax.set_xlabel("period (yr)"); ax.set_ylabel("normalised power")
        ax.set_title(n); ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(a.out, dpi=130)
    print(f"\nplot: {a.out}")


if __name__ == "__main__":
    main()
