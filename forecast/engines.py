"""Forecast engines for the GEM-LTE ledger.

lte          Frozen LTE fit: experiments/Feb2026/<dir>/{lt.exe.p, lt.exe.resp,
             <idx>.dat, lte_run.sh} run by the Ada solver itself (lt.exe, test
             mode, FORECAST_YEARS) on a date grid extended past the end of the data.  Nothing is refit: the
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
import subprocess
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


def lte_binary() -> Path:
    """The Ada solver: $LTE_EXE, else experiments/Feb2026/lt.exe, else the gprbuild output enso_opt."""
    for p in (os.environ.get("LTE_EXE"), FEB / "lt.exe", FEB / "enso_opt"):
        if p and Path(p).exists():
            return Path(p).resolve()
    raise FileNotFoundError("no LTE binary: build with `gprbuild -P lte.gpr enso_opt` or set LTE_EXE")


def read_resp(path: Path) -> dict:
    """lt.exe.resp as GEM.Read_Response_File sees it: KEY VALUE pairs, quoted values kept whole, later keys win."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-")):
            continue
        parts = shlex.split(line)
        if len(parts) >= 2:
            out[parts[0]] = " ".join(parts[1:])
    return out


def sha256_files(paths):
    h = hashlib.sha256()
    for p in paths:
        p = Path(p)
        h.update(p.name.encode())
        h.update(p.read_bytes() if p.exists() else b"<missing>")
    return h.hexdigest()


def run_config(run_dir: Path, index: str, extra_env: dict | None = None):
    """Reproduce the exact lte_run.sh / lte_gui.py run environment; extra_env (the
    target's `env:` in targets.yml, e.g. TIDE_LEGACY for fits made before the
    2026-10-06 tide fix) is applied last."""
    resp = read_resp(run_dir / "lt.exe.resp")
    env = dict(RUN_DEFAULTS)
    for k in ("TRAIN_START", "TRAIN_END", "EXCLUDE", "TREND", "F9", "METRIC"):
        if k in resp:
            env[k] = resp[k]
    sh = run_dir / "lte_run.sh"
    if sh.exists():
        for line in sh.read_text().splitlines():
            if line.startswith("export "):
                for item in shlex.split(line[7:]):
                    k, _, v = item.partition("=")
                    env[k.strip()] = v.strip()
    for k in ("TIMEOUT", "TEST_ONLY", "NUMBER_OF_PROCESSORS"):
        env.pop(k, None)
    env.update({k: str(v) for k, v in (extra_env or {}).items()})
    env.setdefault("CLIMATE_INDEX", f"{index}.dat")
    return resp, env


def lte_layout(t: dict):
    """Run directory, leg number (None for a single fit) and the data file scored, for an
    lte target.  `dir: <name>` is experiments/Feb2026/<name>/<name>.dat; `legs_dir: <path>`
    (repo-relative) is a multi-leg run (LEGS N, lt.exe.resp.1..N) scored on leg `leg`."""
    if "legs_dir" in t:
        run_dir = REPO / t["legs_dir"]
        leg = int(t.get("leg", 1))
        dat = run_dir / read_resp(run_dir / f"lt.exe.resp.{leg}")["CLIMATE_INDEX"]
        return run_dir, leg, dat
    run_dir = FEB / t["dir"]
    return run_dir, None, run_dir / f"{t['dir']}.dat"


def _run_files(run_dir: Path):
    """Files lt.exe reads from the run directory: parameters, settings (with every leg's),
    the data of every leg (or of the single fit), the dLOD reference."""
    resp = read_resp(run_dir / "lt.exe.resp")
    n_legs = int(float(resp.get("LEGS", 0)))
    files = [run_dir / "lt.exe.p", run_dir / "lt.exe.resp"]   # dlod_ref.dat is an lt.exe output
    for k in range(1, n_legs + 1):
        rk = run_dir / f"lt.exe.resp.{k}"
        files += [rk, run_dir / read_resp(rk)["CLIMATE_INDEX"]]
    return n_legs, files


def lte_fingerprint(run_dir: Path, index: str | None, extra_env: dict | None = None):
    n_legs, files = _run_files(run_dir)
    if not n_legs:
        files.append(run_dir / f"{index}.dat")
    # the solver is identified by its source, not its bytes, so a CI build and a local
    # build of the same commit give the same issue id
    files += [run_dir / "lte_run.sh", Path(__file__), REPO / "lte.gpr"]
    files += sorted((REPO / "src").glob("*.ad[bs]"))
    h = sha256_files(files)
    if extra_env:
        h = hashlib.sha256((h + json.dumps(extra_env, sort_keys=True)).encode()).hexdigest()
    return h


def _run_lt(run_dir: Path, env: dict, dat: Path | None, forecast_years: float) -> dict:
    """One TEST_ONLY run of lt.exe on a scratch copy of run_dir.  Returns {leg: array} of
    lte_results(.N).csv (date, model, data as fitted, manifold, ...); leg 0 = single fit."""
    n_legs, files = _run_files(run_dir)
    if dat is not None:
        files.append(dat)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td); wd = td / "run"; wd.mkdir()
        for f in files:
            if f.exists():
                (wd / f.name).write_bytes(f.read_bytes())
        (td / "lt.exe").symlink_to(lte_binary())
        if (FEB / "dlod3.dat").exists():
            for where in (td, wd):
                (where / "dlod3.dat").symlink_to((FEB / "dlod3.dat").resolve())
        full = {**os.environ, **{k: str(v) for k, v in env.items()},
                "TEST_ONLY": "TRUE", "TIMEOUT": "120",
                "NUMBER_OF_PROCESSORS": str(max(1, n_legs)), "FORECAST_YEARS": f"{forecast_years:g}"}
        if dat is not None:
            full["CLIMATE_INDEX"] = dat.name
        r = subprocess.run("ulimit -s 65536 2>/dev/null; ../lt.exe", shell=True, cwd=wd, env=full,
                           stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=1800)
        names = {k: f"lte_results.{k}.csv" for k in range(1, n_legs + 1)} if n_legs else {0: "lte_results.csv"}
        out = {}
        for k, n in names.items():
            if not (wd / n).exists():
                raise RuntimeError(f"lt.exe produced no {n} for {run_dir.name}:\n{r.stdout[-2000:]}")
            out[k] = np.loadtxt(wd / n, delimiter=",")
        return out


def lte_project(run_dir: Path, index: str | None, horizon_years: float = 5.0,
                extra_env: dict | None = None, leg: int | None = None):
    """Return dict(dates, model, data, segment, meta), computed by the Ada lt.exe itself.

    lt.exe is run in test mode (nothing is fitted or searched) with FORECAST_YEARS set:
    it extends the date grid past the last sample on the data's own step, computes the
    manifold and the model there, and keeps the extension out of the regression, the
    data smoother and the scores. A second run without the extension checks that the
    fit on the real samples is unchanged.

    For a multi-leg run (leg = N) all legs are run together, as they were fitted; leg N is
    the one scored, and every leg's projection is returned in `legs` for the record."""
    resp, env = run_config(run_dir, index or "", extra_env)
    if leg is not None:
        env.pop("CLIMATE_INDEX", None)
        lresp = read_resp(run_dir / f"lt.exe.resp.{leg}")
        dat = run_dir / lresp["CLIMATE_INDEX"]
        sampling = float(lresp.get("SAMPLING", resp.get("SAMPLING", 12)))
        ext_all = _run_lt(run_dir, env, None, horizon_years)
        ref_all = _run_lt(run_dir, env, None, 0.0)
        ext, ref = ext_all[leg], ref_all[leg]
    else:
        dat = run_dir / env["CLIMATE_INDEX"]
        sampling = float(resp.get("SAMPLING", env.get("SAMPLING", 12)))
        ext_all = _run_lt(run_dir, env, dat, horizon_years); ref_all = _run_lt(run_dir, env, dat, 0.0)
        ext, ref = ext_all[0], ref_all[0]
    raw = np.loadtxt(dat)
    dates, data = raw[:, 0], raw[:, 1]
    n_real = len(dates)
    if len(ref) != n_real or not np.allclose(ext[:n_real, 0], dates, atol=1e-6):
        raise RuntimeError(f"{run_dir.name}: lt.exe rows do not match the data file ({len(ref)} vs {n_real})")
    all_d, model = ext[:, 0], ext[:, 1]
    data_f = ref[:, 2]                      # the data as lt.exe fitted it (after its F9 passes)
    # the extension must not alter the fit; the last few samples may legitimately move,
    # since the model's own smoother now sees values beyond the edge
    devs = []
    for k in ref_all:
        n_k = len(ref_all[k]); kk = max(n_k - 4, 1)
        devs.append(float(np.max(np.abs(ext_all[k][:kk, 1] - ref_all[k][:kk, 1])) /
                          (float(np.std(ref_all[k][:kk, 1])) or 1.0)))

    eff = {**env, **{kk: v for kk, v in resp.items() if kk in env}}
    ts, te = float(eff["TRAIN_START"]), float(eff["TRAIN_END"])
    exclude = str(eff["EXCLUDE"]).lower() in ("true", "t", "1")
    seg = np.where(all_d > dates[-1], "post_data",
                   np.where((all_d >= ts) & (all_d <= te),
                            "holdout" if exclude else "train",
                            "train" if exclude else "outside_train"))
    hmask = (dates >= ts) & (dates <= te)
    tmask = ~hmask if exclude else hmask
    ok = data_f != 0                        # lt.exe drops zero-valued samples
    def _r(m):
        m = m & ok
        return float(np.corrcoef(model[:n_real][m], data_f[m])[0, 1]) if m.sum() > 3 else None
    meta = dict(engine="lte", solver=str(lte_binary().name), index=dat.stem,
                run_dir=str(run_dir.relative_to(REPO)), leg=leg, data_file=dat.name,
                env=eff, sampling=sampling,
                data_start=float(dates[0]), data_end=float(dates[-1]), n_data=int(n_real),
                horizon_end=float(all_d[-1]),
                r_train=_r(tmask), r_holdout=_r(hmask), r_all=_r(np.ones(n_real, bool)),
                extension_max_dev_over_std=max(devs))
    n_fut = len(all_d) - n_real
    out = dict(dates=all_d, model=model, data=np.r_[data, np.full(n_fut, np.nan)],
               data_filtered=np.r_[data_f, np.full(n_fut, np.nan)],
               segment=seg, meta=meta)
    if leg is not None:
        out["legs"] = {k: v[:, :2] for k, v in ext_all.items()}
    return out


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
