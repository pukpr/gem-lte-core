#!/usr/bin/env python3
"""winding_rank.py — automated ranked winding numbers for a series, given a
manifold. Companion to winding_scalogram.py (shares its transform).

Unlike the scalogram (visual), this ranks candidate winding numbers
objectively, with the null model upgraded from white noise to AR(1)
persistence surrogates (rho=0.97, matched length) — the red-noise floor
that white-noise normalization was shown to underestimate (see
WINDING_SCALOGRAM_FEASIBILITY.md).

Ridge qualification (all three required):
  peak  : max log2(power / AR1-floor) >= --min-snr  (default 2.0 bits)
  sharp : ridge FWHM in M            <= --max-fwhm  (default 0.12)
  steady: fraction of time windows whose local ridge (within +/-0.06 of M,
          above --cont-bits=1 bit) sits ON M (+/-0.03) >= --min-cont (0.70)

Input is either an existing analysis subdir (--index DIR: ranks col 3 Data
against col 4 Forcing of its lte_results.csv) or a bare two-column .dat
series (--dat FILE) remapped onto the manifold of a reference subdir
(--manifold-from DIR, default nino4). In both cases only the *data* panel
is ranked, never the Model panel (no circularity).

    ./winding_rank.py --index nino4
    ./winding_rank.py --dat /tmp/soi.dat --manifold-from nino4
    ./winding_rank.py --index pdo_iaaft_detuned_surrogate   # expect: none pass
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import standardize          # noqa: E402
from winding_scalogram import winding_transform    # noqa: E402

ROOT = Path(__file__).resolve().parent


def read_results_csv(path: Path):
    """(year, data, forcing) from an lte_results.csv."""
    data = np.loadtxt(path, delimiter=",", usecols=(0, 2, 3))
    year, obs, forcing = data.T
    return year, obs, forcing


def read_dat(path: Path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def ar1_floor(t, forcing, m_grid, t0_grid, sigma, n_reps=24, rho=0.97,
              seed=3):
    """Per-M power floor from AR(1) persistence surrogates (see docstring)."""
    rng = np.random.default_rng(seed)
    n = len(t)
    floor = np.zeros(len(m_grid))
    s = np.sqrt(1 - rho * rho)
    for _ in range(n_reps):
        e = rng.standard_normal(n)
        z = np.zeros(n)
        for i in range(1, n):
            z[i] = rho * z[i - 1] + s * e[i]
        Gn, _, _ = winding_transform(t, standardize(z), forcing, m_grid,
                                     t0_grid, sigma)
        floor += np.mean(np.abs(Gn) ** 2, axis=1)
    return floor / n_reps


def find_ridges(m_grid, prof, lp, min_snr, max_fwhm):
    """Candidate local maxima of the time-mean |G| profile, greedily
    de-duplicated within +/-max_fwhm, passing the peak and sharpness tests.
    Yields (M, peak_bits, fwhm)."""
    dm = m_grid[1] - m_grid[0]
    order = np.argsort(-prof)
    taken = []
    for i in order:
        if prof[i] < 0.2 * prof.max():
            break
        if any(abs(m_grid[i] - t) < max_fwhm for t in taken):
            continue
        pk = lp[i].max()
        taken.append(m_grid[i])
        if pk < min_snr:
            continue
        half = prof[i] / 2
        lo = i
        while lo > 0 and prof[lo - 1] >= half and abs(m_grid[lo - 1] -
                                                      m_grid[i]) <= max_fwhm:
            lo -= 1
        hi = i
        while (hi < len(prof) - 1 and prof[hi + 1] >= half and
               abs(m_grid[hi + 1] - m_grid[i]) <= max_fwhm):
            hi += 1
        yield float(m_grid[i]), float(pk), float((hi - lo) * dm)


def continuity(lp, m_grid, M, cont_bits=1.0, halfband=0.06):
    band = (m_grid >= M - halfband) & (m_grid <= M + halfband)
    if band.sum() < 2:  # near grid edge: fall back to nearest cells
        i0 = int(np.argmin(np.abs(m_grid - M)))
        band = np.zeros(len(m_grid), bool)
        band[max(i0 - 3, 0):min(i0 + 4, len(m_grid))] = True
    idxs = np.where(band)[0]
    on_m = 0
    for j in range(lp.shape[1]):
        col = lp[idxs, j]
        k = int(np.argmax(col))
        if abs(m_grid[idxs[k]] - M) <= 0.03 and col[k] > cont_bits:
            on_m += 1
    return on_m / lp.shape[1]


def rank_series(t, x, forcing, m_max=5.0, dm=0.01, sigma=15.0, t0_step=5.0,
                min_snr=2.0, max_fwhm=0.12, min_cont=0.70, n_reps=24):
    m_grid = np.arange(0.0, m_max + dm / 2, dm)
    t0_grid = np.arange(t[0] + sigma / 2, t[-1] - sigma / 2 + 1e-9, t0_step)
    if len(t0_grid) < 3:
        raise SystemExit(f"record too short for sigma={sigma}yr windows")
    G, _, _ = winding_transform(t, standardize(x), forcing, m_grid, t0_grid,
                                sigma)
    floor = ar1_floor(t, forcing, m_grid, t0_grid, sigma, n_reps=n_reps)
    lp = np.log2(np.maximum(np.abs(G) ** 2 / floor[:, None], 1e-9))
    prof = np.abs(G).mean(axis=1)
    out = []
    for M, pk, fwhm in find_ridges(m_grid, prof, lp, min_snr, max_fwhm):
        cont = continuity(lp, m_grid, M)
        out.append({"M": round(M, 4), "peak_bits": round(pk, 2),
                    "fwhm_M": round(fwhm, 3), "continuity": round(cont, 2),
                    "pass": bool(fwhm <= max_fwhm and cont >= min_cont)})
    out.sort(key=lambda r: -r["peak_bits"])
    return out


def load_case(args):
    """Returns (label, t, x, forcing)."""
    if args.index:
        t, x, forcing = read_results_csv(
            ROOT / args.index / "lte_results.csv")
        if args.manifold_from:
            mf_dir = Path(args.manifold_from)
            if not mf_dir.is_absolute():
                mf_dir = ROOT / mf_dir
            mt, _, mforcing = read_results_csv(mf_dir / "lte_results.csv")
            forcing = np.interp(t, mt, mforcing)
        return args.index, t, x, forcing
    dat = Path(args.dat)
    t, x = read_dat(dat)
    mf_dir = Path(args.manifold_from or "nino4")
    if not mf_dir.is_absolute():
        mf_dir = ROOT / mf_dir
    mt, _, mforcing = read_results_csv(mf_dir / "lte_results.csv")
    m = (t >= mt[0]) & (t <= mt[-1])
    t, x = t[m], x[m]
    forcing = np.interp(t, mt, mforcing)
    return dat.stem, t, x, forcing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=
                                 argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--index", help="analysis subdir with lte_results.csv")
    g.add_argument("--dat", help="bare two-column (year value) .dat file")
    ap.add_argument("--manifold-from", default=None,
                    help="subdir whose lte_results.csv col 4 supplies the "
                         "manifold (default: own col 4, or nino4 for --dat)")
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--sigma", type=float, default=15.0)
    ap.add_argument("--t0-step", type=float, default=5.0)
    ap.add_argument("--min-snr", type=float, default=2.0,
                    help="min peak log2 power / AR1 floor (bits)")
    ap.add_argument("--max-fwhm", type=float, default=0.12)
    ap.add_argument("--min-cont", type=float, default=0.70)
    ap.add_argument("--json", action="store_true", help="JSON output")
    args = ap.parse_args()

    label, t, x, forcing = load_case(args)
    rows = rank_series(t, x, forcing, args.m_max, args.dm, args.sigma,
                       args.t0_step, args.min_snr, args.max_fwhm,
                       args.min_cont)
    n_win = int((t[-1] - t[0] - args.sigma) / args.t0_step) + 1
    if args.json:
        print(json.dumps({"label": label, "n_windows": n_win,
                          "ridges": rows}, indent=1))
        return
    print(f"{label}: {len(rows)} candidate ridge(s) over {n_win} time "
          f"windows — null = AR1 rho=0.97; pass needs peak>={args.min_snr} "
          f"bits, fwhm<={args.max_fwhm}, cont>={args.min_cont}")
    print(f"{'M':>8} {'peak_bits':>10} {'fwhm_M':>7} {'cont':>5}  pass")
    for r in rows:
        print(f"{r['M']:>8.4f} {r['peak_bits']:>10.2f} {r['fwhm_M']:>7.3f} "
              f"{r['continuity']:>5.2f}  {'YES' if r['pass'] else 'no'}")
    passed = [r["M"] for r in rows if r["pass"]]
    if passed:
        print(f"\nRANKED winding numbers (by peak power): {passed}")
    else:
        print("\nRANKED winding numbers: NONE — no persistent sharp ridge "
              "above the red-noise null")


if __name__ == "__main__":
    main()
