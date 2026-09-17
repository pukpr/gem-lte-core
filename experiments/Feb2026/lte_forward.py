#!/usr/bin/env python3
"""lte_forward.py — sequential, non-optimizing equivalent of the combined
Ada (lt.exe / enso_opt) + Python (lte_gui.py) software program.

Rationale
---------
The full system runs an Ada random-descent *optimizer* (N threads, up to
MAXLOOPS restarts) wrapped around a deterministic *forward computation*:
given a parameter set (lt.exe.p) and a climate index (.dat), the pipeline
produces the forced climate response and dumps it to lte_results.csv.

This script keeps only the forward computation, as one Python pipeline, so
the salient transformations are visible without optimizer noise:

    1. Doodson periods            (GEM.LTE.Doodson / Year_Adjustment)
    2. Jerked tidal factors       (Jerked_Tidal_Factors)
    3. Tidal superposition        (Tide_Sum_Diff, Scaling=0)
    4. Monthly impulse comb       (Impulse_Delta, Amplify)
    5. IIR integrator             (IIR: memory + sign ramp + pre-history)
    6. Bessel forcing modulation  (Bessel, NM>1 branch)
    7. Least-squares regression   (Regression_Factors: level, k0, mode
                                   amp/phase, annual harmonics, trend)
    8. LTE standing-wave response (LTE: sin(2pi*m*F+ph) with NonLin power)
    9. Smoothing                  (Median + 2x FIR, or F9 Filter9Point)
   10. Results dump               (Save -> lte_results.csv layout)

There is NO random search, NO Monitor/thread coordination, NO metric-based
restart or reset.  One parameter set in -> one model output out.  The
optimizer's inner "evaluate" step is exactly this code path, so running it
with the optimizer's final parameters reproduces lte_results.csv;
--verify checks that equality column by column.

Every routine names its Ada counterpart under src/*.adb.

Usage
-----
    ./lte_forward.py amo --verify
    ./lte_forward.py amo --cv 1880 1885

Config precedence mirrors GEM.Getenv (src/gem.adb):
    <dir>/lt.exe.resp  >  environment  >  built-in default
('-'-prefixed resp lines are commented out; later duplicates override.)
Built-in defaults match the lte_gui.py / lte_run.py run configuration
(TRAIN_START=1880, TRAIN_END=1885, EXCLUDE=true, TREND=true, F9=1,
METRIC=CC, SAMPLING=12).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import sys
from pathlib import Path
from typing import overload

import numpy as np

ROOT = Path(__file__).resolve().parent  # .../experiments/Feb2026

# ---------------------------------------------------------------------------
# Configuration  (GEM.Getenv / Read_Response_File, src/gem.adb)
# ---------------------------------------------------------------------------

def read_resp(path: Path) -> dict:
    """Parse <exe>.resp into {NAME: str}; later duplicates override;
    '-'-prefixed names are commented out."""
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
    """cmd-line override > .resp file > environment > default.
    Result is coerced to the type of `default` (bool/int/float/str)."""

    def __init__(self, resp: dict, overrides: dict):
        self.resp = resp
        self.overrides = overrides

    @overload
    def get(self, name: str, default: bool) -> bool: ...
    @overload
    def get(self, name: str, default: int) -> int: ...
    @overload
    def get(self, name: str, default: float) -> float: ...
    @overload
    def get(self, name: str, default: str) -> str: ...
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


# ---------------------------------------------------------------------------
# GEM.LTE (src/gem-lte.ads/.adb): Doodson periods
# ---------------------------------------------------------------------------

DRACONIC = 27.212220815
TROPICAL = 27.321661554
ANOMALISTIC = 27.554549886
N_PERIOD = 1.0 / (1.0 / DRACONIC - 1.0 / TROPICAL)      # nodal  ~346.6 d
P_PERIOD = 1.0 / (1.0 / TROPICAL - 1.0 / ANOMALISTIC)   # perigee ~3232 d
YEAR_IN_DAYS = 365.2422484
JERK_REFERENCE_PERIOD = 13.66083077  # Mf fortnightly

# Doodson_Args: (s, h, p, N, Year_Multiplier)  -- 42 constituents
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
    """GEM.LTE.Year_Length: base + startup YEAR (+resp/env YEAR) + candidate."""
    return YEAR_IN_DAYS + year_startup + dynamic


def doodson_periods(year_startup: float, dynamic: float) -> np.ndarray:
    """GEM.LTE.Doodson / Year_Adjustment -> constituent periods (days)."""
    yl = year_length(year_startup, dynamic)
    return np.array([
        1.0 / (s / TROPICAL + (h * hm) / yl + p / P_PERIOD + n / N_PERIOD)
        for (s, h, p, n, hm) in DOODSON_ARGS
    ])


# ---------------------------------------------------------------------------
# Forward primitives (src/gem-lte-primitives.adb)
# ---------------------------------------------------------------------------

def tide_sum(dates, ap, periods, year_len, scaling=0.0, integ=0.0):
    """Tide_Sum_Diff (Cos_Phase branch):
        y = sum_j [ amp_j*cos(2pi*f_j*t + ph_j) - integ*f_j*sin(2pi*f_j*t + ph_j) ]
    with f_j = Year_Len/period_j (cycles/yr).  Calc_Forcing calls it with
    Scaling=0, so Ada's TF1*TF2 cross term is absent (and the partition into
    TF1/TF2 loops is then irrelevant)."""
    y = np.zeros_like(dates, dtype=float)
    for j in range(len(periods)):
        freq = year_len / periods[j]
        theta = 2.0 * math.pi * freq * dates + ap[j, 1]
        y += ap[j, 0] * np.cos(theta) - integ * freq * np.sin(theta)
    return y


def impulse_delta(dates, del_a, del_b, asym, sampling):
    """Impulse_Delta: monthly Dirac comb.
    Trunc = int(frac(t)*SAMPLING); DPos = int(|delB|*SAMPLING) — Ada's
    Integer(x) conversion ROUNDS to nearest, so use np.rint, not trunc."""
    trunc = np.rint((dates - np.floor(dates)) * sampling).astype(int)
    dpos = int(round(abs(del_b) * sampling))
    v = np.zeros_like(dates)
    v[trunc == dpos] = del_a
    m2 = (dpos + sampling // 2) % 12               # Ada mod 12 literal
    v[trunc == m2] = asym
    return v


def iir(raw, lag_a, lag_c, init, start_date, dates):
    """IIR integrator:  y[i] = x[i] + Mem*y[i-1] - copysign(lag_c, y[i-1]),
    seeded with `init` at the first sample with Date > start_date, then a
    backward pass y[i-1] = -x[i-1] + Mem*y[i] + copysign(lag_c, y[i])
    reconstructs pre-history.  Mem = lag_a clamped to [0, 1]."""
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
    """Bessel: v += eS sin(2pi k v) + eC cos(2pi k v) + eS2*eS sin(4pi k v)
    + eC2*eC cos(4pi k v).  In the NM>1 branch, the Ada solver supplies
    eS=impA, eC=impB, k=M(NM), eS2=offs, and eC2=bg."""
    return (v + e_s * np.sin(2 * math.pi * k * v)
            + e_c * np.cos(2 * math.pi * k * v)
            + e_s2 * e_s * np.sin(4 * math.pi * k * v)
            + e_c2 * e_c * np.cos(4 * math.pi * k * v))


def regression_factors(t, fv, y, m, trend_on, third=0.0):
    """Regression_Factors: OLS of data y onto design matrix
        [1, F, sin(2pi m_k F), cos(2pi m_k F) for k, (+6 annual/trend cols)]
    returning (level, k0, amp[NM], phase[NM], (Semi1,Semi2,Ann1,Ann2),
    trend, accel).  The Annual indices alias the top modulation columns when
    trend is off (Ada quirk) — reproduced so outputs match bit-for-bit."""
    nm = len(m)
    cols = [np.ones_like(fv), fv]
    for k in range(nm):
        e = np.exp(fv * third)
        cols.append(np.sin(2 * math.pi * m[k] * fv) * e)
        cols.append(np.cos(2 * math.pi * m[k] * fv) * e)
    if trend_on:
        # Ada col order (1-based): Num-5=cos(2pi t), Num-4=sin(2pi t),
        # Num-3=cos(4pi t), Num-2=sin(4pi t), Num-1=t, Num=(t-t0)^2
        cols += [np.cos(2 * math.pi * t), np.sin(2 * math.pi * t),
                 np.cos(4 * math.pi * t), np.sin(4 * math.pi * t),
                 t, (t - t[0]) ** 2.0]
    A = np.column_stack(cols)
    ncol = A.shape[1]
    # Ada Regression_Coefficients: MLR.Inverse(A^T A) * A^T * y (normal eqs).
    coef = np.linalg.solve(A.T @ A, A.T @ y)
    level, k0 = float(coef[0]), float(coef[1])
    amp = np.empty(nm)
    phase = np.empty(nm)
    for k in range(nm):                            # Ada: K starts at 4 (1-based)
        c_sin, c_cos = coef[2 + 2 * k], coef[3 + 2 * k]
        amp[k] = math.hypot(c_sin, c_cos)
        # GNAT: Arctan(Coefficients(K), Coefficients(K-1)) = atan2(Y=c_cos,
        # X=c_sin); matches A*sin(th+ph) = c_sin*sin(th) + c_cos*cos(th).
        phase[k] = math.atan2(c_cos, c_sin)
    semi1, semi2 = float(coef[ncol - 3]), float(coef[ncol - 4])  # Coef(Num-2/-3)
    ann1, ann2 = float(coef[ncol - 5]), float(coef[ncol - 6])    # Coef(Num-4/-5)
    if trend_on:
        trend, accel = float(coef[ncol - 2]), float(coef[ncol - 1])
    else:
        trend, accel = 0.0, 0.0
    return level, k0, amp, phase, (semi1, semi2, ann1, ann2), trend, accel


def lte_response(fv, dates, m, amp, phase, level, k0, trend, accel,
                 nonlin, annual, third=0.0):
    """LTE: y = sum_j amp_j * sgn(sin(2pi m_j F + ph_j)) * |sin(.)|^NonLin
        + level + k0*F + trend*date + accel*(date - date_first)^2
        + Ann1 sin(2pi d) + Ann2 cos(2pi d) + Sem1 sin(4pi d) + Sem2 cos(4pi d)"""
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


def median3(v):
    """Median: 3-point running median; first/last unchanged."""
    r = v.copy()
    for i in range(1, len(v) - 1):
        r[i] = sorted((v[i - 1], v[i], v[i + 1]))[1]
    return r


def fir(v, behind, current, ahead):
    """FIR: 3-point weighted filter; first/last unchanged."""
    r = v.copy()
    r[1:-1] = behind * v[:-2] + current * v[1:-1] + ahead * v[2:]
    return r


def filter9point(v):
    """Filter9Point: (1,2,1)/4 boxcar; first/last unchanged."""
    r = v.copy()
    r[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]
    return r


def me_power_spectrum(forcing_x, signal, f_start, f_step, n_out):
    """ME_Power_Spectrum loop body: spectrum[j] = S^2 + C^2 where
    S,C = sum over i>=First+8 of sin/cos(2pi F_j Forcing_i) * signal_i."""
    x = forcing_x[8:]
    y = signal[8:]
    spec = np.empty(n_out)
    f = f_start
    for j in range(n_out):
        a = 2 * math.pi * f * x
        s = float(np.dot(np.sin(a), y))
        c = float(np.dot(np.cos(a), y))
        spec[j] = s * s + c * c
        f += f_step
    return spec


def window_box(v, lobe):
    """Window: running mean over 2*lobe+1; edges unchanged."""
    r = v.copy()
    k = 2 * lobe + 1
    cs = np.concatenate(([0.0], np.cumsum(v)))
    for i in range(lobe, len(v) - lobe):
        r[i] = (cs[i + lobe + 1] - cs[i - lobe]) / k
    return r


def cc(x, y, dates, cc_start, cc_end):
    """CC: Pearson r with Ada's leading/trailing zero-value trimming."""
    start, stop = 0, len(x) - 1
    for i in range(len(x)):
        if x[i] == 0.0 or y[i] == 0.0:
            start = i + 1
        else:
            break
    for i in range(len(x) - 1, -1, -1):
        if x[i] == 0.0 or y[i] == 0.0:
            stop = i - 1
        else:
            break
    xs, ys = x[start:stop + 1], y[start:stop + 1]
    ds = dates[start:stop + 1]
    keep = (ds > cc_start) & (ds < cc_end) & (xs != 0.0) & (ys != 0.0)
    xs, ys = xs[keep], ys[keep]
    n = len(xs)
    if n == 0:
        return 0.0
    den = ((n * (xs * xs).sum() - xs.sum() ** 2)
           * (n * (ys * ys).sum() - ys.sum() ** 2))
    if den <= 0.0:
        return 0.0
    return float((n * (xs * ys).sum() - xs.sum() * ys.sum()) / math.sqrt(den))


def dtw_distance(x, y, window_size=3):
    """DTW_Distance: Sakoe-Chiba-banded Dynamic Time Warping distance,
    normalized to a CC-like scale (higher is better, 1.0 = perfect
    alignment). Faithful port of gem-lte-primitives.adb's DTW_Distance /
    nested Distance, including its quirks:
      - a y[j] == 0.0 point costs nothing (treated as a free move, e.g. for
        masked/invalid data points, matching Ada's `if Y(J).Value = 0.0
        then Cost := 0.0`)
      - the "diagonal"/"left" neighbor lookup clamps to index 0 instead of
        treating cells before the band as unreachable (Ada's
        `Integer'Max (J - 1, Y'First)`)
      - the row buffer (`curr`) is only reset at index 0 each row; any
        column outside the current band keeps its value from an earlier
        row, exactly as Ada's full-array `DTW_Previous := DTW_Current`
        copies whatever `DTW_Current` last held there — used only via the
        clamped neighbor lookups above, not as a stray global minimum.

    Less sensitive to timing/alignment than Pearson's cc(): a warped but
    correctly-shaped match still scores well.

    Distance(x, y) is normalized against Distance(-y, y) — the DTW cost of
    y against its own mirror image, used as a "worst case" reference scale:
        (Distance(-y, y) - Distance(x, y)) / Distance(-y, y)
    """
    def _distance(xs, ys, w):
        n = len(xs)
        inf = float("inf")
        prev = [inf] * n
        curr = [inf] * n
        prev[0] = 0.0
        last_j = 0
        for i in range(n):
            curr[0] = inf  # DTW_Current(X'First) reset every row
            lo = max(0, i - w)
            hi = min(n - 1, i + w)
            for j in range(lo, hi + 1):
                cost = 0.0 if ys[j] == 0.0 else abs(xs[i] - ys[j])
                j_prev = max(j - 1, 0)
                up = prev[j] if j > 0 else inf
                curr[j] = cost + min(prev[j_prev], curr[j_prev], up)
                last_j = j
            prev = curr[:]  # DTW_Previous := DTW_Current (full-array copy)
        return prev[last_j]

    xs = [float(v) for v in x]
    ys = [float(v) for v in y]
    if not xs or not ys:
        return 0.0
    w = max(1, int(window_size))
    max_d = _distance([-v for v in ys], ys, w)
    if max_d <= 0.0:
        return 0.0
    return (max_d - _distance(xs, ys, w)) / max_d


# ---------------------------------------------------------------------------
# The pipeline: Dipole_Model inner evaluation, executed exactly once
# ---------------------------------------------------------------------------

def forward(cfg: Config, params: dict):
    dat = Path(cfg.get("CLIMATE_INDEX", "amo.dat"))
    if not dat.is_absolute():
        dat = Path.cwd() / dat
    raw = np.loadtxt(dat)
    dates, data = raw[:, 0].copy(), raw[:, 1].copy()
    data0 = data.copy()
    n = len(dates)

    # --- parameters from lt.exe.p (Shared.Read_JSON) ---
    lpap = np.array(params["lpap"], dtype=float)   # rows: [period, amp, phase]
    lt = np.array(params["ltep"], dtype=float)
    year_cand = float(params.get("year", 0.0))
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
          "ann1", "ann2", "IR", "ma", "mp", "shfT", "init")}

    # --- configuration on the forward path ---
    year_startup = cfg.get("YEAR", 0.0)
    jerk = cfg.get("JERK", 0.0)
    sampling = int(cfg.get("SAMPLING", 12.0))
    filt = cfg.get("FILTER", 0.33333333)
    f9 = cfg.get("F9", 1)
    idate = cfg.get("IDATE", 0.0)
    nonlin = cfg.get("NONLIN", 1.0)
    cc_start = cfg.get("CC_START", 0.0)
    cc_end = cfg.get("CC_END", 999999999.0)
    # METRIC follows the same env-var convention as Ada's Min_Entropy flag
    # (GEM.Getenv("METRIC", "")): "DTW" swaps the fitted metric (the "cc"
    # result the optimizer maximizes) from Pearson's cc() to dtw_distance(),
    # which is less sensitive to alignment/timing than a plain correlation.
    metric = str(cfg.get("METRIC", "CC")).strip().upper()
    dtw_window = int(cfg.get("DTW_WINDOW", 3))
    f_start = cfg.get("FSTART", 0.01)
    f_step = cfg.get("FSTEP", 0.18)
    nm = cfg.get("NM", len(lt))
    nh_str = cfg.get("NH", "")           # "" (resp has NH "") -> no harmonics
    harms = [int(float(h)) for h in str(nh_str).split()]
    exclude = cfg.get("EXCLUDE", True)
    trend_on = cfg.get("TREND", True)

    # (Ada line ~815: Data_Records := Filter9Point(Data_Records) mutates the
    # working series BEFORE regression, and Save() writes that same filtered
    # series — verified: reference column 3 equals the 1x-filtered .dat, and
    # the recovered coefficients come from regressing on it.)
    if f9 > 0:
        for _ in range(f9):
            data = filter9point(data)
    data_fit = data

    # training window (Find_Index semantics; 1-based Ada -> 0-based here)
    ts = cfg.get("TRAIN_START", dates[0])
    te = cfg.get("TRAIN_END", dates[-1])
    first = int(np.searchsorted(dates, ts, side="right")) - 1
    last = int(np.searchsorted(dates, te, side="right")) - 1
    first = max(first, 0)
    last = max(last, 0)

    # 1) Doodson periods (Year_Adjustment with candidate year correction)
    periods = doodson_periods(year_startup, year_cand)
    assert np.allclose(np.abs(lpap[:, 0]), np.abs(periods), rtol=0.01), \
        "lt.exe.p lpap periods disagree with the Doodson set"

    # 2) JERK: rotate each constituent into a mix of tide and its derivative
    ap = lpap[:, 1:3].copy()                      # (Amplitude, Phase)
    if jerk != 0.0:
        if not 0.0 <= jerk <= 1.0:
            raise SystemExit("JERK must be between 0.0 and 1.0")
        w_norm, w_der = 1.0 - jerk, jerk * JERK_REFERENCE_PERIOD / periods
        ap[:, 0] *= np.sqrt(w_norm ** 2 + w_der ** 2)
        ap[:, 1] += np.arctan2(w_der, w_norm)

    # 3-5) tide sum -> impulse comb -> IIR integrator  (Calc_Forcing)
    tf = tide_sum(dates, ap, periods, year_length(year_startup, year_cand),
                  scaling=0.0, integ=B["shfT"])
    forcing = iir(tf * impulse_delta(dates, B["delA"], B["delB"], B["asym"],
                                     sampling),
                  lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
                  start_date=idate, dates=dates)

    # 6) Bessel modulation of the forcing manifold (NM>1, no Lock_Freq).
    # This is the exact argument mapping in Dipole_Model: the final LTE
    # modulation is the Bessel frequency and the second-harmonic coefficients
    # are Offset ("offs") and bg, not the fitted annual coefficients.
    m = np.concatenate([lt[:nm], [lt[nm - 1] * h for h in harms]])
    forcing = bessel(forcing, B["impA"], B["impB"], m[nm - 1],
                     B["offs"], B["bg"])

    # 7) least-squares regression over the (excluded) training set
    if exclude:
        sel = np.r_[0:first + 1, last:n]
    else:
        sel = np.arange(first, last + 1)
    # Filter out zero-valued data points (they cause singular regression)
    nonzero = data[sel] != 0.0
    sel = sel[nonzero]
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[sel], forcing[sel], data[sel], m, trend_on)
    # NOTE: Ada applies Monotonic_Increase (abs) only AFTER the Model is
    # computed — it affects the reported values, not the series.
    model = lte_response(forcing, dates, m, amp, phase, level, k0,
                         trend, accel, nonlin, annual)
    if cfg.get("SECULAR", True):
        trend, accel = abs(trend), abs(accel)
    if B["IR"] != 0.0:
        model = model.copy()
        for i in range(n - 1, 11, -1):
            model[i] -= B["IR"] * model[i - 12]

    # 9) smoothing
    if f9 > 0:
        for _ in range(f9):
            model = filter9point(model)
    else:
        model = median3(model)
        model = fir(fir(model, filt, 1.0 - 2.0 * filt, filt),
                    filt, 1.0 - 2.0 * filt, filt)

    # 10) metric + power-spectrum columns (Metric, Save)
    if metric == "DTW":
        corr = dtw_distance(model[sel], data[sel], dtw_window)
    else:
        corr = cc(model[sel], data[sel], dates[sel], cc_start, cc_end)
    finite = np.isfinite(model) & np.isfinite(data)
    pearson = float(np.corrcoef(model[finite], data[finite])[0, 1])
    fspec = f_start + f_step * np.arange(n)
    mspec = window_box(me_power_spectrum(forcing, model,
                                         f_start, f_step, n), 2)
    dspec = window_box(me_power_spectrum(forcing, data,
                                         f_start, f_step, n), 2)
    return dict(dates=dates, model=model, data=data, forcing=forcing,
                fspec=fspec, mspec=mspec, dspec=dspec, cc=corr,
                pearson=pearson,
                periods=periods, m=m, nm=nm, nh=len(harms), f9=f9,
                reg=dict(level=level, k0=k0, amp=amp.tolist(),
                         phase=phase.tolist(), annual=annual,
                         trend=trend, accel=accel))


def fspec_x(grid, forcing):
    """ME_Power_Spectrum is driven by 'Forcing' values on the model loop and
    by Forcing values on the data loop; the date grid is only used by FT_CC.
    For Save() both spectra use the *forcing* series as the phase variable —
    this identity helper keeps the call sites readable."""
    return forcing


def write_results(path: Path, r: dict):
    """Save(Safe.Save) layout: date, model, data, forcing, f, model_spec,
    data_spec — ', '-joined, 15 significant digits, Ada 'Img spacing."""
    with open(path, "w") as f:
        for i in range(len(r["dates"])):
            row = [r["dates"][i], r["model"][i], r["data"][i],
                   r["forcing"][i], r["fspec"][i], r["mspec"][i],
                   r["dspec"][i]]
            f.write(", ".join(f" {v: .14E}" for v in row) + "\n")


def verify(expected: Path, r: dict) -> bool:
    """Column-by-column agreement with the Ada program's lte_results.csv.
    Tolerances reflect accumulated 64-bit transcendental rounding through the
    pipeline (sin/cos of large arguments, regression), not logic differences."""
    exp = np.loadtxt(expected, delimiter=",")
    got = np.column_stack([r["dates"], r["model"], r["data"], r["forcing"],
                           r["fspec"], r["mspec"], r["dspec"]])
    if exp.shape != got.shape:
        print(f"shape mismatch: {exp.shape} vs {got.shape}")
        return False
    names = ["date", "model", "data", "forcing", "spec_freq",
             "spec_model", "spec_data"]
    ok = True
    print(f"{'column':<12}{'rows':>6}{'max|diff|':>12}{'rel':>10}"
          f"{'corr':>12}   status")
    for j, name in enumerate(names):
        e, g = exp[:, j], got[:, j]
        d = float(np.max(np.abs(e - g)))
        scale = max(float(np.max(np.abs(e))), 1e-30)
        rel = d / scale
        corr = float(np.corrcoef(e, g)[0, 1]) if np.std(g) > 0 else 1.0
        good = rel <= 1e-5 and corr > 1.0 - 1e-9
        ok &= good
        print(f"{name:<12}{len(e):>6}{d:>12.2e}{rel:>10.1e}{corr:>12.9f}"
              f"   {'MATCH' if good else 'DIFF'}")
    print(f"\nCC(train) = {r['cc']:.6f}")
    return ok


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Sequential (non-optimizing) LTE forward pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index", nargs="?", default="amo",
                    help="index subdirectory of Feb2026 (default amo)")
    ap.add_argument("--cv", nargs=2, metavar=("START", "END"),
                    default=("1880", "1885"),
                    help="TRAIN_START TRAIN_END (default: 1880 1885)")
    ap.add_argument("--param", default=None,
                    help="parameter JSON (default <dir>/lt.exe.p)")
    ap.add_argument("--resp", default=None,
                    help="resp file (default <dir>/lt.exe.resp)")
    ap.add_argument("--metric", default=None)
    ap.add_argument("--trend", choices=["true", "false"], default=None)
    ap.add_argument("--exclude", choices=["true", "false"], default=None)
    ap.add_argument("--f9", type=int, default=None)
    ap.add_argument("--out", default="lte_forward_results.csv")
    ap.add_argument("--verify", action="store_true",
                    help="compare against lte_results.csv in the index dir")
    args = ap.parse_args()

    run_dir = (ROOT / args.index).resolve()
    if not run_dir.is_dir():
        print(f"error: {run_dir} not found", file=sys.stderr)
        return 2
    os.chdir(run_dir)

    resp = read_resp(Path(args.resp) if args.resp else run_dir / "lt.exe.resp")
    resp.setdefault("CLIMATE_INDEX", f"{args.index}.dat")
    overrides = {}
    if args.metric:
        overrides["METRIC"] = args.metric
    if args.trend:
        overrides["TREND"] = args.trend
    if args.exclude:
        overrides["EXCLUDE"] = args.exclude
    if args.f9 is not None:
        overrides["F9"] = str(args.f9)
    overrides["TRAIN_START"], overrides["TRAIN_END"] = args.cv

    params = json.loads(Path(args.param or run_dir / "lt.exe.p").read_text())
    r = forward(Config(resp, overrides), params)

    write_results(Path(args.out), r)
    print(f"forward pass complete: {run_dir / args.out}  "
          f"({len(r['dates'])} samples, NM={r['nm']}, NH={r['nh']}, "
          f"F9={r['f9']})")
    reg = r["reg"]
    print(f"regression: level={reg['level']:.5f} k0={reg['k0']:.5f} "
          f"trend={reg['trend']:.5f} accel={reg['accel']:.4g}")
    print(f"amp   = {np.round(reg['amp'], 5).tolist()}")
    print(f"phase = {np.round(reg['phase'], 5).tolist()}")
    print(f"annual= {np.round(reg['annual'], 6).tolist()}")
    print(f"Pearson correlation (model, data) = {r['pearson']:.6f}")

    if args.verify:
        exp = run_dir / "lte_results.csv"
        if not exp.is_file():
            print("no lte_results.csv to verify against", file=sys.stderr)
            return 1
        return 0 if verify(exp, r) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
