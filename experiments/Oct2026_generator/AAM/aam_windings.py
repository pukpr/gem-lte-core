#!/usr/bin/env python3
"""Winding spectra of NINO4 and AAM series on one shared manifold F(t).

    python3 aam_windings.py RUN_DIR [--series name=file.dat ...] [--kmax 0.5] [--surr 200]

RUN_DIR is a 2-leg run (leg 1 NINO4) whose lte_results.1.csv column 4 is the
shared manifold F(t). For each series X (default: the run's own leg data,
column 3, plus any --series .dat files sampled on the same months) and
winding k, over the months where every series has data:

    G_X(k) = sum_t X(t) exp(-i 2 pi k F(t)) / N        (X mean- and trend-removed)

|G_X(k)| is the winding spectrum (the projection onto the model's own basis
sin(2 pi k F + phase)). Each spectrum is tested against phase-randomized
surrogates of X (same time-domain power spectrum, alignment with F broken):
z(k) = (|G| - mean_surr) / sd_surr, and the fraction of k bins above the
surrogates' 95th percentile (5% expected by chance).

Pairs: winding coherence between NINO4 and each other series,
C(k) = |sum_w G_N G_X*| / sqrt(sum |G_N|^2 sum |G_X|^2) over a +-dk window,
and the phase of G_X relative to G_N at the run's fitted windings
(lt.exe.windings.*.json): an AAM that is ENSO-driven shares NINO4's phases.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def load_leg(run, k):
    a = np.loadtxt(run / f"lte_results.{k}.csv", delimiter=",")
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def month_key(t):
    return np.rint(np.asarray(t) * 12).astype(int)


def prep(x, t):
    x = x - np.polyval(np.polyfit(t, x, 1), t)
    return x / x.std()


def spectrum(x, f, ks):
    return np.exp(-2j * np.pi * np.outer(ks, f)) @ x / len(x)


def surrogates(x, n, rng):
    X = np.fft.rfft(x)
    ph = rng.uniform(0, 2 * np.pi, (n, len(X)))
    ph[:, 0] = 0.0
    if len(x) % 2 == 0:
        ph[:, -1] = 0.0
    return np.fft.irfft(np.abs(X) * np.exp(1j * (np.angle(X) + ph)), len(x))


def coherence(gn, gx, half):
    c = np.full(len(gn), np.nan)
    for i in range(half, len(gn) - half):
        s = slice(i - half, i + half + 1)
        c[i] = abs(np.sum(gn[s] * np.conj(gx[s]))) / np.sqrt(np.sum(abs(gn[s]) ** 2) * np.sum(abs(gx[s]) ** 2))
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--series", action="append", default=[], help="name=path.dat, sampled on the run's months")
    ap.add_argument("--kmax", type=float, default=0.5)
    ap.add_argument("--dk", type=float, default=0.001)
    ap.add_argument("--surr", type=int, default=200)
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--exclude", nargs=2, type=float, default=None, help="drop a held-out window")
    ap.add_argument("--manifold", default=None,
                    help="lte_results csv whose column 4 is used as F(t) instead of the run's own "
                         "(e.g. a NINO4-only fit, so that the AAM series is out of sample for F)")
    a = ap.parse_args()
    run = Path(a.run)
    rng = np.random.default_rng(1)

    t1, _, n4, F = load_leg(run, 1)
    t2, _, x2, _ = load_leg(run, 2)
    label2 = run.name.split("_")[1].upper()
    series = {"NINO4": dict(zip(month_key(t1), n4)), label2: dict(zip(month_key(t2), x2))}
    for s in a.series:
        name, path = s.split("=", 1)
        d = np.loadtxt(path)
        series[name] = dict(zip(month_key(d[:, 0]), d[:, 1]))
    if a.manifold:
        m = np.loadtxt(a.manifold, delimiter=",")
        t1, F = m[:, 0], m[:, 3]
    fmap = dict(zip(month_key(t1), F))
    keys = sorted(set.intersection(*(set(k for k, v in s.items() if v != 0) for s in series.values())) & set(fmap))
    t = np.array(keys) / 12.0
    sel = np.ones(len(t), bool)
    if a.start is not None:
        sel &= t >= a.start
    if a.end is not None:
        sel &= t <= a.end
    if a.exclude:
        sel &= ~((t >= a.exclude[0]) & (t <= a.exclude[1]))
    keys = [k for k, s in zip(keys, sel) if s]
    t = t[sel]
    f = np.array([fmap[k] for k in keys])
    ks = np.arange(a.dk, a.kmax + a.dk / 2, a.dk)
    span = np.ptp(f)
    print(f"{run.name}: {len(t)} months {t[0]:.2f}-{t[-1]:.2f}; F range {span:.2f} "
          f"(winding resolution ~{1 / span:.3f}); k 0-{a.kmax}")

    W = [json.load(open(run / f"lt.exe.windings.{i}.json")) for i in (1, 2)]
    fit_k = [abs(k) for k, _, _ in W[0]["k_amp_phase"]]

    G, Z = {}, {}
    print(f"\n{'series':8s} {'frac>95%':>8s}  {'max z':>6s} at k   z at fitted windings " +
          " ".join(f"{k:.4f}" for k in fit_k))
    for name, s in series.items():
        x = prep(np.array([s[k] for k in keys]), t)
        g = spectrum(x, f, ks)
        sg = np.abs(np.array([spectrum(y, f, ks) for y in surrogates(x, a.surr, rng)]))
        z = (abs(g) - sg.mean(0)) / sg.std(0)
        p95 = np.percentile(sg, 95, axis=0)
        G[name], Z[name] = g, z
        zf = [z[np.argmin(abs(ks - k))] for k in fit_k]
        print(f"{name:8s} {np.mean(abs(g) > p95):8.3f}  {z.max():6.2f} at {ks[np.argmax(z)]:.3f}   " +
              " ".join(f"{v:+6.2f}" for v in zf))

    half = max(1, int(round(2.0 / span / a.dk)))   # +-2 resolution cells
    print(f"\nwinding coherence with NINO4 (window +-{half * a.dk:.3f} in k) and phase at fitted windings:")
    for name in list(series)[1:]:
        c = coherence(G["NINO4"], G[name], half)
        x = prep(np.array([series[name][k] for k in keys]), t)
        c0 = [np.nanmean(coherence(G["NINO4"], spectrum(y, f, ks), half)) for y in surrogates(x, 50, rng)]
        ph = [np.degrees(np.angle(G[name][i] / G["NINO4"][i])) for i in (np.argmin(abs(ks - k)) for k in fit_k)]
        both = (Z["NINO4"] > 2) & (Z[name] > 2)
        only = (Z[name] > 2) & (Z["NINO4"] < 1)
        print(f"  {name:8s} mean coherence {np.nanmean(c):.3f} (surrogate null {np.mean(c0):.3f} +- {np.std(c0):.3f}); "
              "phase vs NINO4 at fitted k: " +
              " ".join(f"{p:+5.0f}" for p in ph) +
              f"; bins z>2 in both {both.mean():.3f}, {name} only {only.mean():.3f}")
        if only.any():
            idx = np.where(only)[0]
            groups = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            peaks = sorted(((Z[name][g].max(), ks[g][np.argmax(Z[name][g])]) for g in groups), reverse=True)[:6]
            print("    strongest windings in " + name + " but not NINO4: " +
                  ", ".join(f"k={k:.3f} (z {zz:.1f})" for zz, k in peaks))

    print("\nfitted windings (lt.exe.windings.*.json), amplitude (phase deg):")
    for j, (k, _, _) in enumerate(W[0]["k_amp_phase"]):
        a1, p1 = W[0]["k_amp_phase"][j][1:]
        a2, p2 = W[1]["k_amp_phase"][j][1:]
        print(f"  k={k:+.4f}: NINO4 {a1:.3f} ({np.degrees(p1):+4.0f})  {label2} {a2:.3f} ({np.degrees(p2):+4.0f})  "
              f"ratio {a2 / a1:.3f}, phase diff {(np.degrees(p2 - p1) + 180) % 360 - 180:+4.0f}")


if __name__ == "__main__":
    main()
