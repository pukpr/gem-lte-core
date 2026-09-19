#!/usr/bin/env python3
"""ws.py — generalized winding-scalogram fitting tool for any climate
index, synthesizing the two sibling projects developed this session:

  gem-lte-core (this repo): the production-matching LTE forward pipeline
    (tide_sum -> impulse_delta -> iir -> bessel, validated against the
    Ada source), genuine train/holdout cross-validation (not just
    in-sample r), F9 filtering, IR decompensation including the
    UNCOMPENSATED basis-folding mode (found necessary for tpi), and the
    trend/seasonal linear decomposition technique (iode/tna) for
    separating nuisance from genuine signal.

  ~/eval/winding_scalogram: a locally-weighted (Gaussian-kernel)
    regression that lets amp/phase drift smoothly in time -- much richer
    DESCRIPTIVE fits than one global regression -- paired with a
    mandatory honesty discipline: IAAFT phase-randomized surrogates run
    through the identical design, reporting "honest r" = real r minus
    surrogate mean (their own README: a negative-control surrogate
    outscored real PDO on the naive design -- in-sample r alone is not
    evidence of locking).

Self-contained: all LTE primitives are inlined below (not imported from
either repo), so this script runs from any directory given only an
index's own subdirectory (<index>/<index>.dat, <index>/lt.exe.p,
optionally <index>/lt.exe.resp) somewhere findable -- see
`find_index_root`.

Usage
-----
    ./ws.py --index nino4                        # both modes, one index
    ./ws.py --index nino4 --mode global           # production-matching only
    ./ws.py --index nino4 --mode local            # local kernel fit only
    ./ws.py --index brestexcl --root /some/path
    ./ws.py --all                                 # sweep every index found
    ./ws.py --index pna --plot                    # also save a PNG overlay

Modes
-----
  global  Fits regression_factors/lte_response (the exact production
          design) with a genuine held-out window (default: 2000-2010
          excluded from the fit, scored only there). Honors each
          index's own UNCOMPENSATED flag (folds the IR lag-12 operator
          into the design basis and regresses directly against raw
          data) vs the standard DR-pre-emphasis mode. Reports
          train_r/val_r -- the number that should match a real
          production CC.

  local   Kernel-weighted (Gaussian, width --sigma years) regression
          re-fit at every time point, on the index's OWN (possibly
          irregular/gapped) date grid -- no data is fabricated across
          real gaps, unlike re-gridding onto a uniform month axis.
          Reports r_raw, r_var (quadratic-detrended), dCC (first-
          difference correlation), and -- mandatory, not optional --
          "honest_r" = r_raw minus the mean r obtained by running the
          IDENTICAL fit on --n-surrogates IAAFT phase-randomized clones
          of the real data. A small or negative honest_r means the fit
          quality is descriptive texture, not evidence of locking.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import sys
from pathlib import Path

import numpy as np

# ===========================================================================
# Section 1 -- portable LTE forward-model primitives
# (ported and validated against src/gem-lte-primitives*.adb this session;
# see gem-lte-core/experiments/Feb2026/lte_forward.py for the fuller,
# --verify-checked original this is drawn from)
# ===========================================================================

YEAR_IN_DAYS = 365.2422484
DRACONIC = 27.212220815
TROPICAL = 27.321661554
ANOMALISTIC = 27.554549886
N_PERIOD = 1.0 / (1.0 / DRACONIC - 1.0 / TROPICAL)
P_PERIOD = 1.0 / (1.0 / TROPICAL - 1.0 / ANOMALISTIC)
JERK_REFERENCE_PERIOD = 13.66083077

# Doodson_Args: (s, h, p, N, Year_Multiplier) -- 42 constituents. Periods
# MUST be recomputed from these + the current year length, NOT read from
# lt.exe.p's own stored lpap[:,0] column -- found by direct comparison
# that one term's stored/recomputed periods disagree by a factor of ~2
# (a sign-crossing near-zero-denominator case), which silently wrecks the
# whole forcing sum if lpap[:,0] is trusted directly (as a naive port of
# winding_scalogram's marquee.py does -- fine for its own local-kernel-
# regression use since that's flexible enough to mask it, NOT fine for
# global-mode production matching).
DOODSON_ARGS = [
    (1, 0, 0, 0, 1.0), (1, 0, 0, 1, 1.0), (0, 0, 2, 2, 1.0), (2, 0, 0, 1, 1.0),
    (2, 0, 0, 0, 1.0), (2, 0, 0, 2, 1.0), (1, 0, -1, 0, 1.0), (2, 0, -2, 0, 1.0),
    (0, 0, 0, 1, 1.0), (0, 0, 2, 0, 1.0), (1, -2, 1, 0, 1.0), (0, 0, 2, 1, 1.0),
    (1, 0, -1, 1, 1.0), (1, 0, -1, -1, 1.0), (0, 0, 1, 1, 1.0), (1, 0, 1, 1, 1.0),
    (0, 0, 1, -1, 1.0), (0, 0, -1, 0, 1.0), (0, 0, -2, 1, 1.0), (3, 0, -1, 0, 1.0),
    (3, 0, -1, 1, 1.0), (3, 0, -1, 2, 1.0), (0, 0, 0, 2, 1.0), (0, 0, 1, 2, 1.0),
    (3, -2, 1, 0, 1.0), (3, 0, -3, 0, 1.0), (3, -2, 1, 1, 1.0), (4, -2, 0, 1, 1.0),
    (4, 0, -2, 1, 1.0), (4, 0, -2, 0, 1.0), (4, -2, 0, 0, 1.0), (5, -2, -1, 0, 1.0),
    (2, -2, 0, 0, 1.0), (3, -2, -1, 0, 1.0), (1, 0, 1, 0, 1.0), (5, -2, -1, 1, 1.0),
    (5, 0, -3, 0, 1.0), (2, -2, 0, -1, 1.0), (2, -2, 0, 1, 1.0), (2, -3, 0, 0, 1.0),
    (1, 2, -1, 0, 1.0), (0, 2, 0, 0, 1.0),
]


def year_length(year_startup: float, dynamic: float = 0.0) -> float:
    return YEAR_IN_DAYS + year_startup + dynamic


def doodson_periods(year_startup: float, dynamic: float) -> np.ndarray:
    yl = year_length(year_startup, dynamic)
    return np.array([
        1.0 / (s / TROPICAL + (h * hm) / yl + p / P_PERIOD + n / N_PERIOD)
        for (s, h, p, n, hm) in DOODSON_ARGS
    ])


def read_resp(path: Path) -> dict:
    """Parse <exe>.resp into {NAME: str}; later duplicates override;
    '-'-prefixed names are commented out. Returns {} if the file is
    missing (some indices, e.g. in the winding_scalogram layout, don't
    ship one -- Config's own defaults then apply)."""
    opts = {}
    if not path.is_file():
        return opts
    tokens = shlex.split(path.read_text())
    for i in range(0, len(tokens) - 1, 2):
        name = tokens[i]
        if name.startswith("-"):
            continue
        opts[name.upper()] = tokens[i + 1]
    return opts


class Config:
    def __init__(self, resp: dict, overrides: dict):
        self.resp = resp
        self.overrides = overrides

    def get(self, name: str, default):
        if name in self.overrides:
            v = self.overrides[name]
        elif name in self.resp:
            v = self.resp[name]
        elif name in os.environ:
            v = os.environ[name]
        else:
            return default
        if isinstance(default, bool):
            return str(v).strip().upper() in ("TRUE", "T", "1")
        if isinstance(default, int):
            return int(float(v))
        if isinstance(default, float):
            return float(v)
        return str(v)


def tide_sum(dates, ap, periods, year_len, scaling=0.0, integ=0.0):
    y = np.zeros_like(dates, dtype=float)
    for j in range(len(periods)):
        freq = year_len / periods[j]
        theta = 2.0 * math.pi * freq * dates + ap[j, 1]
        y += ap[j, 0] * np.cos(theta) - integ * freq * np.sin(theta)
    return y


def impulse_delta(dates, del_a, del_b, asym, sampling):
    trunc = np.rint((dates - np.floor(dates)) * sampling).astype(int)
    dpos = int(round(abs(del_b) * sampling))
    v = np.zeros_like(dates)
    v[trunc == dpos] = del_a
    m2 = (dpos + sampling // 2) % 12
    v[trunc == m2] = asym
    return v


def iir(raw, lag_a, lag_c, init, start_date, dates):
    """Seeded with `init` at the first sample with Date > start_date,
    then forward recursion + a backward pass reconstructing pre-history
    (this backward pass is what makes backward extrapolation before an
    index's own record start well-defined -- see gem-lte-core's
    brestexcl pre-1880 work)."""
    mem = 1.0 if lag_a > 1.0 else (0.0 if lag_a < 0.0 else lag_a)
    n = len(raw)
    idx = int(np.searchsorted(dates, start_date, side="right"))
    idx = min(max(idx, 0), n - 1)
    y = np.empty(n)
    y[idx] = init
    for i in range(idx + 1, n):
        y[i] = raw[i] + mem * y[i - 1] - math.copysign(lag_c, y[i - 1])
    for i in range(idx, 0, -1):
        y[i - 1] = -raw[i - 1] + mem * y[i] + math.copysign(lag_c, y[i])
    return y


def bessel(v, e_s, e_c, k, e_s2, e_c2):
    return (v + e_s * np.sin(2 * math.pi * k * v)
            + e_c * np.cos(2 * math.pi * k * v)
            + e_s2 * e_s * np.sin(4 * math.pi * k * v)
            + e_c2 * e_c * np.cos(4 * math.pi * k * v))


def filter9point(v):
    r = v.copy()
    for _ in range(1):
        pass
    r[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]
    return r


def regression_factors(t, fv, y, m, trend_on, third=0.0):
    nm = len(m)
    cols = [np.ones_like(fv), fv]
    for k in range(nm):
        e = np.exp(fv * third)
        cols.append(np.sin(2 * math.pi * m[k] * fv) * e)
        cols.append(np.cos(2 * math.pi * m[k] * fv) * e)
    if trend_on:
        cols += [np.cos(2 * math.pi * t), np.sin(2 * math.pi * t),
                 np.cos(4 * math.pi * t), np.sin(4 * math.pi * t),
                 t, (t - t[0]) ** 2.0]
    A = np.column_stack(cols)
    ncol = A.shape[1]
    coef = np.linalg.solve(A.T @ A, A.T @ y)
    level, k0 = float(coef[0]), float(coef[1])
    amp = np.empty(nm)
    phase = np.empty(nm)
    for k in range(nm):
        c_sin, c_cos = coef[2 + 2 * k], coef[3 + 2 * k]
        amp[k] = math.hypot(c_sin, c_cos)
        phase[k] = math.atan2(c_cos, c_sin)
    semi1, semi2 = float(coef[ncol - 3]), float(coef[ncol - 4])
    ann1, ann2 = float(coef[ncol - 5]), float(coef[ncol - 6])
    if trend_on:
        trend, accel = float(coef[ncol - 2]), float(coef[ncol - 1])
    else:
        trend, accel = 0.0, 0.0
    return level, k0, amp, phase, (semi1, semi2, ann1, ann2), trend, accel


def lte_response(fv, dates, m, amp, phase, level, k0, trend, accel,
                  nonlin, annual, third=0.0):
    semi1, semi2, ann1, ann2 = annual
    out = np.zeros_like(fv)
    for j in range(len(m)):
        sw = np.sin(2 * math.pi * m[j] * fv + phase[j]) * np.exp(fv * third)
        pw = np.abs(sw) ** nonlin
        out += amp[j] * np.where(sw < 0.0, -pw, pw)
    out += (level + k0 * fv + trend * dates + accel * (dates - dates[0]) ** 2.0)
    out += (ann1 * np.sin(2 * math.pi * dates) + ann2 * np.cos(2 * math.pi * dates)
            + semi1 * np.sin(4 * math.pi * dates) + semi2 * np.cos(4 * math.pi * dates))
    return out


def one_shot_lag_diff(x: np.ndarray, ir: float) -> np.ndarray:
    """x(t) - ir*x(t-12), first 12 rows unchanged. Used both as the
    standard-mode post-hoc IR decompensation AND, applied to every
    design-matrix column, as the UNCOMPENSATED-mode basis fold (see
    build_design / fit_global_uncompensated)."""
    out = x.copy()
    out[12:] = x[12:] - ir * x[:-12]
    return out


def build_design(t, fv, m, trend_on):
    cols = [np.ones_like(fv), fv]
    for k in range(len(m)):
        cols.append(np.sin(2 * math.pi * m[k] * fv))
        cols.append(np.cos(2 * math.pi * m[k] * fv))
    if trend_on:
        cols += [np.cos(2 * math.pi * t), np.sin(2 * math.pi * t),
                 np.cos(4 * math.pi * t), np.sin(4 * math.pi * t),
                 t, (t - t[0]) ** 2.0]
    return np.column_stack(cols)


# ===========================================================================
# Section 2 -- index discovery / loading (portable: no hardcoded repo path)
# ===========================================================================

def default_search_roots() -> list[Path]:
    roots = [Path.cwd(), Path(__file__).resolve().parent]
    home = Path.home()
    for rel in ("eval/gem-lte-core/experiments/Feb2026",
                "eval/winding_scalogram",
                "eval/winding_scalogram/milestone_2026_09",
                "eval/winding_scalogram/supplemental_2026_09",
                "github/pukpr/GEM-LTE/experiments/Feb2026"):
        p = home / rel
        if p not in roots:
            roots.append(p)
    return roots


def find_index_root(index: str, root_arg: str | None) -> Path:
    candidates = [Path(root_arg)] if root_arg else default_search_roots()
    for base in candidates:
        d = base / index
        if (d / f"{index}.dat").exists() and (d / "lt.exe.p").exists():
            return base
    tried = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"could not find {index}/{index}.dat + {index}/lt.exe.p under any "
        f"of: {tried}  (pass --root explicitly)")


def load_index(root: Path, index: str) -> dict:
    idx_dir = root / index
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    resp = read_resp(idx_dir / "lt.exe.resp")
    overrides = {}
    if "harm" in params:
        overrides["NH"] = " ".join(str(int(float(h))) for h in params["harm"])
    cfg = Config(resp, overrides)

    dat_path = idx_dir / f"{index}.dat"
    raw = np.loadtxt(dat_path)
    dates, data = raw[:, 0].copy(), raw[:, 1].copy()
    data_raw = data.copy()

    lpap = np.array(params["lpap"], dtype=float)
    lt = np.array(params["ltep"], dtype=float)
    year_cand = float(params.get("year", 0.0))
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
          "ann1", "ann2", "IR", "ma", "mp", "shfT", "init")}

    year_startup = cfg.get("YEAR", 0.0)
    jerk = cfg.get("JERK", 0.0)
    sampling = int(cfg.get("SAMPLING", 12.0))
    f9 = cfg.get("F9", 1)
    nonlin = cfg.get("NONLIN", 1.0)
    nm = cfg.get("NM", len(lt))
    nh_str = cfg.get("NH", "")
    harms = [int(float(h)) for h in str(nh_str).split()] if str(nh_str).strip() else []
    trend_on = cfg.get("TREND", True)
    uncompensated = cfg.get("UNCOMPENSATED", False)
    idate = cfg.get("IDATE", 0.0)

    if f9 > 0:
        for _ in range(f9):
            data = filter9point(data)

    yl = year_length(year_startup, year_cand)
    # periods MUST come from doodson_periods(), not lpap[:,0] -- see the
    # module-level comment on DOODSON_ARGS for why lpap[:,0] is wrong.
    periods = doodson_periods(year_startup, year_cand)
    ap = lpap[:, 1:3].copy()
    if jerk != 0.0:
        w_norm, w_der = 1.0 - jerk, jerk * JERK_REFERENCE_PERIOD / periods
        ap[:, 0] *= np.sqrt(w_norm ** 2 + w_der ** 2)
        ap[:, 1] += np.arctan2(w_der, w_norm)

    tf = tide_sum(dates, ap, periods, yl, scaling=0.0, integ=B["shfT"])
    forcing = iir(tf * impulse_delta(dates, B["delA"], B["delB"], B["asym"], sampling),
                  lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
                  start_date=idate, dates=dates)
    m = np.concatenate([lt[:nm], [lt[nm - 1] * h for h in harms]]) if harms else lt[:nm]
    # bessel's own k parameter is the BACKBONE term lt[nm-1], not the last
    # element of the full (base+harmonics) m array -- found by direct
    # comparison against cv_rolling_blocked.prepare(): using m[-1] here
    # silently fed the FM step the wrong winding entirely whenever harms
    # is non-empty (m[-1] is then the last harmonic, e.g. 3.1x rather
    # than the ~0.2 backbone).
    forcing = bessel(forcing, B["impA"], B["impB"], lt[nm - 1], B["offs"], B["bg"])

    return dict(index=index, root=root, dates=dates, data=data, data_raw=data_raw,
                forcing=forcing, m=m, nonlin=nonlin, trend_on=trend_on, f9=f9,
                ir=B["IR"], uncompensated=uncompensated, params=params, resp=resp,
                yl=yl, year_startup=year_startup, year_cand=year_cand, idate=idate,
                nm=nm)


def build_forcing_at(prep: dict, dates_new: np.ndarray) -> np.ndarray:
    """Re-evaluate this index's forcing pipeline at an arbitrary (e.g.
    extended-backward) date array -- pure function of dates/params, no
    real data needed (see gem-lte-core's brestexcl pre-1880 work)."""
    params = prep["params"]
    lpap = np.array(params["lpap"], dtype=float)
    periods = doodson_periods(prep["year_startup"], prep["year_cand"])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp",
          "shfT", "init")}
    k1 = float(np.array(params["ltep"], dtype=float)[prep["nm"] - 1])
    tf = tide_sum(dates_new, lpap[:, 1:3], periods, prep["yl"], 0.0, B["shfT"])
    comb = impulse_delta(dates_new, B["delA"], B["delB"], B["asym"], 12)
    R = iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
            start_date=prep["idate"], dates=dates_new)
    return bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])


# ===========================================================================
# Section 3 -- GLOBAL mode: production-matching regression + genuine holdout
# ===========================================================================

def decompensate_f9(x, ir, f9):
    out = x.copy()
    if ir != 0.0:
        out = one_shot_lag_diff(out, ir)
    for _ in range(f9):
        out = filter9point(out)
    return out


def fit_global(prep: dict, holdout=(2000.0, 2010.0)) -> dict:
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    nz = data != 0.0
    h0, h1 = holdout
    test_idx = np.nonzero((dates >= h0) & (dates < h1) & nz)[0]
    train_idx = np.nonzero(((dates < h0) | (dates >= h1)) & nz)[0]
    if len(train_idx) < 2 * (2 * len(m) + 7) or len(test_idx) < 3:
        # record too short/sparse for this holdout window -- fall back to
        # an 80/20 chronological split so short records (pna, noi, kap...)
        # still get a genuine out-of-sample number
        valid = np.nonzero(nz)[0]
        cut = valid[int(0.8 * len(valid))]
        train_idx = valid[valid < cut]
        test_idx = valid[valid >= cut]

    if prep["uncompensated"]:
        A = build_design(dates, forcing, m, prep["trend_on"])
        A_L = np.apply_along_axis(one_shot_lag_diff, 0, A, prep["ir"])
        beta = np.linalg.solve(A_L[train_idx].T @ A_L[train_idx],
                                A_L[train_idx].T @ data[train_idx])
        raw_modulation = A @ beta
        model = one_shot_lag_diff(raw_modulation, prep["ir"]) if prep["ir"] != 0.0 else raw_modulation
        for _ in range(prep["f9"]):
            model = filter9point(model)
    else:
        level, k0, amp, phase, annual, trend, accel = regression_factors(
            dates[train_idx], forcing[train_idx], data[train_idx], m, prep["trend_on"])
        model = lte_response(forcing, dates, m, amp, phase, level, k0, trend,
                              accel, prep["nonlin"], annual)
        model = decompensate_f9(model, prep["ir"], prep["f9"])

    train_r = float(np.corrcoef(model[train_idx], data[train_idx])[0, 1])
    val_r = float(np.corrcoef(model[test_idx], data[test_idx])[0, 1])
    return dict(train_r=train_r, val_r=val_r, model=model,
                n_train=len(train_idx), n_test=len(test_idx),
                holdout=(float(dates[test_idx[0]]), float(dates[test_idx[-1]]))
                if len(test_idx) else None)


# ===========================================================================
# Section 4 -- LOCAL mode: kernel-weighted regression + mandatory honesty
# ===========================================================================

def local_design(F, t, Ms, nmax, trend_on=True):
    cols = [np.ones_like(F), F]
    for M in Ms:
        th = 2 * math.pi * M
        for n in range(1, nmax + 1):
            cols += [np.sin(n * th * F), np.cos(n * th * F)]
    if trend_on:
        t0 = t - t[0]
        cols += [np.sin(2 * math.pi * t), np.cos(2 * math.pi * t),
                 np.sin(4 * math.pi * t), np.cos(4 * math.pi * t),
                 t0, t0 ** 2]
    return np.column_stack(cols)


def local_fit(F, t, x, sigma_years, Ms, nmax=3, trend_on=True,
              window_mult=3.0, stride=1, support_mask=None,
              return_params=False):
    """Gaussian-kernel-weighted regression re-fit at every sample point,
    on the data's OWN (possibly irregular/gapped) time grid -- no
    fabrication across real gaps, unlike re-gridding onto a uniform
    month axis first.

    Two speed tricks, both checked to be numerically negligible against
    the exact (window_mult=1e6, stride=1) computation before being
    relied on -- needed since the surrogate-honesty check reruns the
    whole fit --n-surrogates times, and an exact O(n^2) computation is
    ~8s for one ~140yr monthly record on its own:
      - truncate each point's window to +/- window_mult*sigma_years
        (weight beyond 3 sigma is <0.02) instead of the full record;
      - fit only every `stride`-th point (the local fit varies smoothly
        on a sigma_years-wide scale, so this is a cheap, safe
        subsample) and linearly interpolate the rest.
    """
    X = local_design(F, t, Ms, nmax, trend_on)
    n = len(t)
    n_M = len(Ms)
    fit_pts = np.empty(n)
    # per-winding, FIRST-harmonic (n=1) coefficient trajectory -- column
    # k*2*nmax + 2/3 (0-indexed within the M-block, offset by the 2
    # leading [1, F] columns) is that winding's own sin/cos(1*theta*F)
    # pair; amp/phase from these two coefficients is the physically
    # meaningful "how big, and what phase, is this winding RIGHT HERE"
    # trajectory -- see local_fit_params_from_coefs below.
    coefs = np.full((n, X.shape[1]), np.nan) if return_params else None
    half_width = window_mult * sigma_years
    idx_range = range(0, n, stride) if stride > 1 else range(n)
    computed = np.zeros(n, dtype=bool)

    def _one_point(i):
        lo = np.searchsorted(t, t[i] - half_width, side="left")
        hi = np.searchsorted(t, t[i] + half_width, side="right")
        wi = np.exp(-0.5 * ((t[lo:hi] - t[i]) / sigma_years) ** 2)
        if support_mask is not None:
            wi = wi * support_mask[lo:hi]
            if wi.sum() <= 0.0:
                return np.nan, None
        Xw = X[lo:hi] * wi[:, None]
        b, *_ = np.linalg.lstsq(Xw, x[lo:hi] * wi, rcond=None)
        return X[i] @ b, b

    for i in idx_range:
        val, b = _one_point(i)
        fit_pts[i] = val
        if return_params and b is not None:
            coefs[i] = b
        computed[i] = True
    if stride > 1 and not computed[-1]:
        i = n - 1
        val, b = _one_point(i)
        fit_pts[i] = val
        if return_params and b is not None:
            coefs[i] = b
        computed[i] = True

    if stride > 1:
        done = np.nonzero(computed)[0]
        fit = np.interp(t, t[done], fit_pts[done])
    else:
        fit = fit_pts

    if not return_params:
        return fit

    params = {}
    done = np.nonzero(computed & ~np.isnan(coefs[:, 0]))[0]
    for k, M in enumerate(Ms):
        c_sin = coefs[:, 2 + k * 2 * nmax]
        c_cos = coefs[:, 2 + k * 2 * nmax + 1]
        amp_done = np.hypot(c_sin[done], c_cos[done])
        phase_done = np.arctan2(c_cos[done], c_sin[done])
        amp_t = np.interp(t, t[done], amp_done)
        # interpolate phase via its unit-circle components, not the raw
        # angle, to avoid a spurious jump artifact at the +/-pi wraparound
        phase_t = np.arctan2(np.interp(t, t[done], np.sin(phase_done)),
                              np.interp(t, t[done], np.cos(phase_done)))
        params[M] = dict(amp=amp_t, phase=phase_t)
    return fit, params


def quad_detrend(v, t):
    tc = t - t.mean()
    V = np.vander(tc, 3)
    b, *_ = np.linalg.lstsq(V, v, rcond=None)
    return v - V @ b


def iaaft(y: np.ndarray, seed: int) -> np.ndarray:
    """Iterative Amplitude Adjusted Fourier Transform surrogate: a
    phase-randomized clone matching both the power spectrum AND the
    empirical value distribution of y -- the standard null model for
    "could an unrelated series with the same texture score this well."
    """
    rng = np.random.default_rng(seed)
    n = len(y)
    amp = np.abs(np.fft.rfft(y - y.mean()))
    ph = rng.uniform(0, 2 * np.pi, len(amp))
    ph[0] = 0.0
    s = np.fft.irfft(amp * np.exp(1j * ph), n=n)
    y_sorted = np.sort(y)
    for _ in range(30):
        s = np.interp(np.argsort(np.argsort(s)), np.arange(n), y_sorted)
        f = np.fft.rfft(s - s.mean())
        s = np.fft.irfft(amp * np.exp(1j * np.angle(f)), n=n)
    return s


def fit_local(prep: dict, sigma_years=10.0, nmax=3, n_surrogates=12,
              seed0=4001, window_mult=3.0, stride=2) -> dict:
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    nz = data != 0.0
    t = dates[nz]
    x = data[nz]
    x = (x - x.mean()) / x.std()
    F = forcing[nz]
    Ms = sorted({round(abs(mm), 4) for mm in m if 0.002 < abs(mm) < 10})

    y = local_fit(F, t, x, sigma_years, Ms, nmax=nmax, trend_on=prep["trend_on"],
                  window_mult=window_mult, stride=stride)
    r_raw = float(np.corrcoef(x, y)[0, 1])
    xd, yd = quad_detrend(x, t), quad_detrend(y, t)
    r_var = float(np.corrcoef(xd, yd)[0, 1])
    dCC = float(np.corrcoef(np.diff(x), np.diff(y))[0, 1])

    surr_rs = []
    for j in range(n_surrogates):
        xs = iaaft(x, seed0 + j)
        xs = (xs - xs.mean()) / xs.std()
        ys = local_fit(F, t, xs, sigma_years, Ms, nmax=nmax, trend_on=prep["trend_on"],
                       window_mult=window_mult, stride=stride)
        surr_rs.append(float(np.corrcoef(xs, ys)[0, 1]))
    surr_rs = np.array(surr_rs)
    honest_r = r_raw - float(surr_rs.mean())

    return dict(r_raw=r_raw, r_var=r_var, dCC=dCC, honest_r=honest_r,
                surrogate_mean=float(surr_rs.mean()), surrogate_std=float(surr_rs.std()),
                Ms=Ms, nmax=nmax, sigma=sigma_years, t=t, x=x, y=y)


# ===========================================================================
# Section 5 -- CLI / driver
# ===========================================================================

KNOWN_INDICES = [
    "nino4", "nino34", "pdo", "amo", "tna", "nao", "iode", "emi", "npi",
    "pna", "noi", "tpi", "baltic", "darwin1880", "qbo30", "qbo50",
    "brestexcl", "kap10-10-20-30",
]


def run_one(index: str, root_arg: str | None, mode: str, sigma: float,
            nmax: int, n_surrogates: int, holdout, plot: bool) -> dict | None:
    try:
        root = find_index_root(index, root_arg)
        prep = load_index(root, index)
    except Exception as exc:  # noqa: BLE001 -- sweep mode must keep going
        print(f"[{index}] SKIPPED: {exc}")
        return None

    out = {"index": index, "root": str(root)}
    print(f"[{index}] root={root}  m={np.round(prep['m'], 4).tolist()}  "
          f"uncompensated={prep['uncompensated']}  ir={prep['ir']:.4f}  "
          f"f9={prep['f9']}  nonlin={prep['nonlin']}")

    if mode in ("global", "both"):
        g = fit_global(prep, holdout=holdout)
        out["global"] = {k: v for k, v in g.items() if k != "model"}
        print(f"    global: train_r={g['train_r']:+.4f}  val_r={g['val_r']:+.4f}  "
              f"holdout={g['holdout']}  n_train={g['n_train']} n_test={g['n_test']}")
        if plot:
            _plot_global(prep, g, root)

    if mode in ("local", "both"):
        loc = fit_local(prep, sigma_years=sigma, nmax=nmax, n_surrogates=n_surrogates)
        out["local"] = {k: v for k, v in loc.items() if k not in ("t", "x", "y")}
        print(f"    local:  r_raw={loc['r_raw']:+.4f}  r_var={loc['r_var']:+.4f}  "
              f"dCC={loc['dCC']:+.4f}  honest_r={loc['honest_r']:+.4f}  "
              f"(surrogate {loc['surrogate_mean']:+.3f}+/-{loc['surrogate_std']:.3f}, "
              f"n={n_surrogates})  Ms={loc['Ms']}")
        if plot:
            _plot_local(prep, loc, root)

    return out


def _plot_global(prep, g, root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dates, data = prep["dates"], prep["data"]
    nz = data != 0.0
    model = g["model"]
    fig, ax = plt.subplots(figsize=(11, 3.5))
    ax.plot(dates[nz], data[nz], color="0.35", lw=0.9, label="real data")
    ax.plot(dates[nz], model[nz], color="teal", lw=1.0, alpha=0.85,
            label="global (production-matching) fit")
    if g["holdout"] is not None:
        ax.axvspan(*g["holdout"], color="gold", alpha=0.15, label="held out")
    ax.set_title(f"{prep['index']}: GLOBAL mode  train_r={g['train_r']:+.3f}  "
                 f"val_r={g['val_r']:+.3f}", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_path = Path(root) / prep["index"] / f"{prep['index']}_ws_global.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"    [saved] {out_path}")


def _plot_local(prep, loc, root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 3.5))
    ax.plot(loc["t"], loc["x"], color="0.6", lw=1.0, label="obs (standardized)")
    ax.plot(loc["t"], loc["y"], color="crimson", lw=1.1, alpha=0.85,
            label=f"local fit (sigma={loc['sigma']:g}yr, n<={loc['nmax']})")
    ax.set_title(f"{prep['index']}: r_raw={loc['r_raw']:+.3f}  "
                  f"r_var={loc['r_var']:+.3f}  honest_r={loc['honest_r']:+.3f}",
                  fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_path = Path(root) / prep["index"] / f"{prep['index']}_ws_local.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"    [saved] {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", help="index name (subdirectory with <index>.dat + lt.exe.p)")
    ap.add_argument("--root", help="explicit root directory containing <index>/ "
                                    "(default: search CWD, this script's dir, and "
                                    "known repo locations under $HOME)")
    ap.add_argument("--mode", choices=["global", "local", "both"], default="both")
    ap.add_argument("--sigma", type=float, default=10.0, help="local mode: kernel width, years")
    ap.add_argument("--nmax", type=int, default=3, help="local mode: max harmonic order per winding")
    ap.add_argument("--n-surrogates", type=int, default=12, help="local mode: IAAFT surrogate count")
    ap.add_argument("--holdout", nargs=2, type=float, default=[2000.0, 2010.0],
                     metavar=("START", "END"), help="global mode: held-out window")
    ap.add_argument("--plot", action="store_true", help="save a local-mode overlay PNG per index")
    ap.add_argument("--all", action="store_true",
                     help="sweep every index in KNOWN_INDICES found under --root/known roots")
    ap.add_argument("--out-json", help="write the full result set to this JSON path")
    args = ap.parse_args()

    if not args.index and not args.all:
        ap.error("pass --index NAME or --all")

    indices = KNOWN_INDICES if args.all else [args.index]
    results = {}
    for idx in indices:
        r = run_one(idx, args.root, args.mode, args.sigma, args.nmax,
                    args.n_surrogates, tuple(args.holdout), args.plot)
        if r is not None:
            results[idx] = r

    if args.all and results:
        print("\n=== summary ===")
        header = f"{'index':14s}"
        if args.mode in ("global", "both"):
            header += f"{'train_r':>10s}{'val_r':>10s}"
        if args.mode in ("local", "both"):
            header += f"{'r_raw':>10s}{'r_var':>10s}{'honest_r':>10s}"
        print(header)
        for idx, r in results.items():
            line = f"{idx:14s}"
            if "global" in r:
                line += f"{r['global']['train_r']:+10.4f}{r['global']['val_r']:+10.4f}"
            if "local" in r:
                line += (f"{r['local']['r_raw']:+10.4f}{r['local']['r_var']:+10.4f}"
                         f"{r['local']['honest_r']:+10.4f}")
            print(line)

    if args.out_json:
        with open(args.out_json, "w") as f:
            json.dump(results, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        print(f"\n[written] {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
