"""Forecast engines for the GEM-LTE ledger.

lte          Frozen LTE fit: experiments/Feb2026/<dir>/{lt.exe.p, lt.exe.resp,
             <idx>.dat, lte_run.sh} run through lte_forward.forward() on a date
             grid extended past the end of the data.  Nothing is refit: the
             tidal forcing is deterministic, so the projection beyond the data
             is a pure out-of-sample prediction of the committed parameter set.

tidal_lod    LOD: least-squares fit of the Doodson tidal constituents carried
             in a donor lt.exe.p (the LOD-calibrated set) plus annual/semi-
             annual and a slow quadratic, on daily IERS LOD up to the issue
             date; projected forward.  Refit at each monthly issue.

alias_cw     Chandler wobble, forced-aliasing hypothesis (geoenergymath.com,
             2026-09-17): complex polar motion m = x - i y fitted with prograde
             terms at the alias period 365.242/(2*frac(365.242/27.2122)) days
             and the annual period, plus mean pole drift.  Refit monthly.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shlex
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
FEB = REPO / "experiments" / "Feb2026"

# GUI/CLI defaults from experiments/Feb2026/lte_run.py DEFAULTS
RUN_DEFAULTS = {"TRAIN_START": "1880", "TRAIN_END": "1885", "METRIC": "CC",
                "IDATE": "1920.9", "EXCLUDE": "true", "TREND": "true",
                "F9": "1", "DLOD_REF": "TRUE"}


def _load_lte_forward():
    spec = importlib.util.spec_from_file_location("lte_forward", FEB / "lte_forward.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_files(paths):
    h = hashlib.sha256()
    for p in paths:
        p = Path(p)
        h.update(p.name.encode())
        h.update(p.read_bytes() if p.exists() else b"<missing>")
    return h.hexdigest()


def run_config(run_dir: Path, index: str):
    """Reproduce the exact lte_run.sh / lte_gui.py run environment."""
    lf = _load_lte_forward()
    resp = lf.read_resp(run_dir / "lt.exe.resp")
    env = dict(RUN_DEFAULTS)
    for k in ("TRAIN_START", "TRAIN_END", "EXCLUDE", "TREND", "F9", "METRIC"):
        if k in resp:
            env[k] = resp[k]
    sh = run_dir / "lte_run.sh"
    if sh.exists():
        for line in sh.read_text().splitlines():
            if line.startswith("export "):
                k, _, v = line[7:].partition("=")
                env[k.strip()] = v.strip().strip('"')
    env.pop("TIMEOUT", None)
    env.setdefault("CLIMATE_INDEX", f"{index}.dat")
    return lf, resp, env


def lte_fingerprint(run_dir: Path, index: str):
    files = [run_dir / "lt.exe.p", run_dir / "lt.exe.resp", run_dir / f"{index}.dat",
             run_dir / "lte_run.sh", FEB / "lte_forward.py", Path(__file__)]
    return sha256_files(files)


def lte_project(run_dir: Path, index: str, horizon_years: float = 5.0):
    """Return dict(dates, model, data, segment, meta)."""
    lf, resp, env = run_config(run_dir, index)
    dat = run_dir / env["CLIMATE_INDEX"]
    raw = np.loadtxt(dat)
    dates, data = raw[:, 0], raw[:, 1]
    sampling = int(float(resp.get("SAMPLING", 12)))
    step = 1.0 / sampling
    n_fut = int(round(horizon_years * sampling))
    fut = dates[-1] + step * np.arange(1, n_fut + 1)
    ext = np.column_stack([np.r_[dates, fut], np.r_[data, np.zeros(n_fut)]])

    # Edge-hold the data smoother so the appended zeros (which the Ada
    # regression already ignores: data == 0 is dropped) do not bleed into
    # the last real samples via Filter9Point.
    n_real = len(dates)
    orig_f9 = lf.filter9point

    def f9_edge(v):
        if len(v) == len(ext) and np.all(v[n_real:] == 0.0):
            w = v.copy(); w[n_real:] = v[n_real - 1]
            r = orig_f9(w); r[n_real:] = 0.0; r[n_real - 1] = orig_f9(v[:n_real])[-1]
            return r
        return orig_f9(v)

    lf.filter9point = f9_edge
    params = json.loads((run_dir / "lt.exe.p").read_text())
    saved_env = {k: os.environ.get(k) for k in env}
    cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as td:
            ext_path = Path(td) / dat.name
            np.savetxt(ext_path, ext, fmt="%.6f\t%.8g")
            # Ada/GEM.Getenv precedence: lt.exe.resp > environment > default,
            # so the run environment goes into os.environ, not overrides.
            os.environ.update({k: str(v) for k, v in env.items()})
            overrides = {"CLIMATE_INDEX": str(ext_path)}
            os.chdir(run_dir)
            r_ext = lf.forward(lf.Config(resp, overrides), params)
            lf.filter9point = orig_f9
            overrides["CLIMATE_INDEX"] = str(dat)
            r_ref = lf.forward(lf.Config(resp, overrides), params)
    finally:
        os.chdir(cwd)
        lf.filter9point = orig_f9
        for k, v in saved_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    model = r_ext["model"]
    # Interior agreement check: the extension must not alter the fit
    # (ignore the last ~4 samples, which the symmetric model smoother
    # legitimately updates now that forcing exists beyond the edge).
    k = max(n_real - 4, 1)
    dev = float(np.max(np.abs(model[:k] - r_ref["model"][:k])))
    scale = float(np.std(r_ref["model"][:k])) or 1.0

    eff = {**env, **{k: v for k, v in resp.items() if k in env}}
    meta_env = eff
    ts, te = float(eff["TRAIN_START"]), float(eff["TRAIN_END"])
    exclude = str(eff["EXCLUDE"]).lower() in ("true", "t", "1")
    all_d = ext[:, 0]
    seg = np.where(all_d > dates[-1], "post_data",
                   np.where((all_d >= ts) & (all_d <= te),
                            "holdout" if exclude else "train",
                            "train" if exclude else "outside_train"))
    # holdout skill on the committed data (cross-validation interval)
    data_f = data.copy()
    for _ in range(int(float(eff.get("F9", 1)))):
        data_f = orig_f9(data_f)
    hmask = (dates >= ts) & (dates <= te)
    tmask = ~hmask if exclude else hmask
    def _r(m):
        return float(np.corrcoef(model[:n_real][m], data_f[m])[0, 1]) if m.sum() > 3 else None
    meta = dict(engine="lte", index=index, run_dir=str(run_dir.relative_to(REPO)),
                env=meta_env, sampling=sampling, data_start=float(dates[0]),
                data_end=float(dates[-1]), n_data=int(n_real),
                horizon_end=float(all_d[-1]),
                r_train=_r(tmask), r_holdout=_r(hmask), r_all=_r(np.ones(n_real, bool)),
                extension_max_dev_over_std=dev / scale)
    return dict(dates=all_d, model=model, data=np.r_[data, np.full(n_fut, np.nan)],
                data_filtered=np.r_[data_f, np.full(n_fut, np.nan)],
                segment=seg, meta=meta)


# ---------------------------------------------------------------------------
# Harmonic engines (LOD, Chandler wobble)
# ---------------------------------------------------------------------------

def donor_tidal_periods(donor_index="amo"):
    p = json.loads((FEB / donor_index / "lt.exe.p").read_text())
    return np.abs(np.array(p["lpap"], float)[:, 0])  # days


def _design(t_years, periods_days, t0, poly=2):
    cols = [np.ones_like(t_years)]
    tt = t_years - t0
    for k in range(1, poly + 1):
        cols.append(tt ** k)
    for P in periods_days:
        w = 2 * np.pi * 365.2425 / P
        cols += [np.cos(w * t_years), np.sin(w * t_years)]
    return np.column_stack(cols)


def rayleigh_dedup(periods_days, span_years):
    """Keep constituents resolvable over the window: frequency spacing
    >= 1/span (Rayleigh criterion) and period <= span/2.  Longer periods are
    left to the low-order polynomial."""
    span_d = span_years * 365.2425
    keep = []
    for P in sorted(periods_days):
        if P > span_d / 2:
            continue
        if all(abs(1 / P - 1 / Q) >= 1 / span_d for Q in keep):
            keep.append(P)
    return keep


def tidal_lod_project(d, v, issue, train_years=12.0, horizon_years=1.0,
                      donor_index="amo", poly=2, ridge=1e-6):
    tidal = [P for P in donor_tidal_periods(donor_index) if P > 5.0]
    periods = rayleigh_dedup(list(tidal) + [365.2425, 182.62125], train_years)
    sel = (d <= issue) & (d > issue - train_years)
    t0 = issue
    X = _design(d[sel], periods, t0, poly)
    sc = np.linalg.norm(X, axis=0)
    Xs = X / sc
    A = Xs.T @ Xs + ridge * len(Xs) * np.eye(Xs.shape[1])
    coef = np.linalg.solve(A, Xs.T @ v[sel]) / sc
    step = 1.0 / 365.2425
    fut = d[sel][-1] + step * np.arange(1, int(horizon_years * 365.2425) + 1)
    grid = np.r_[d[sel], fut]
    model = _design(grid, periods, t0, poly) @ coef
    seg = np.where(grid > d[sel][-1], "post_data", "train")
    meta = dict(engine="tidal_lod", donor=donor_index, n_constituents=len(periods),
                periods_days=[float(p) for p in periods], train_years=train_years,
                data_end=float(d[sel][-1]), horizon_end=float(grid[-1]),
                r_train=float(np.corrcoef(model[:sel.sum()], v[sel])[0, 1]))
    return dict(dates=grid, model=model, data=np.r_[v[sel], np.full(len(fut), np.nan)],
                data_filtered=np.r_[v[sel], np.full(len(fut), np.nan)],
                segment=seg, meta=meta)


def alias_period_days(year=365.242, draconic=27.2122):
    r = year / draconic
    return year / (r - np.floor(r)) / 2.0


def alias_cw_project(d, x, y, issue, train_years=24.0, horizon_years=2.0,
                     component="x"):
    """Fit complex prograde terms; return projection for one component."""
    Pcw = alias_period_days()
    sel = (d <= issue) & (d > issue - train_years)
    t = d[sel]; m = x[sel] - 1j * y[sel]
    def X(tt):
        cols = [np.ones_like(tt), tt - issue]
        for P in (Pcw, 365.242):
            cols.append(np.exp(1j * 2 * np.pi * 365.2425 / P * tt))
        return np.column_stack(cols).astype(complex)
    coef, *_ = np.linalg.lstsq(X(t), m, rcond=None)
    fut = t[-1] + (1 / 365.2425) * np.arange(1, int(horizon_years * 365.2425) + 1)
    grid = np.r_[t, fut]
    mm = X(grid) @ coef
    model = mm.real if component == "x" else -mm.imag
    obs = x[sel] if component == "x" else y[sel]
    seg = np.where(grid > t[-1], "post_data", "train")
    meta = dict(engine="alias_cw", component=component, chandler_period_days=float(Pcw),
                train_years=train_years, data_end=float(t[-1]), horizon_end=float(grid[-1]),
                amp_cw_arcsec=float(abs(coef[2])), amp_annual_arcsec=float(abs(coef[3])),
                r_train=float(np.corrcoef(model[:sel.sum()], obs)[0, 1]))
    return dict(dates=grid, model=model, data=np.r_[obs, np.full(len(fut), np.nan)],
                data_filtered=np.r_[obs, np.full(len(fut), np.nan)],
                segment=seg, meta=meta)
