#!/usr/bin/env python3
"""GEM-LTE forecast ledger.

    python3 forecast/ledger.py all            # observe + issue + score + report
    python3 forecast/ledger.py issue [-t soi qbo30 ...]
    python3 forecast/ledger.py observe | score | report
    python3 forecast/ledger.py list

Layout (all committed to git, so git history is the timestamp):

    forecast/targets.yml                       registry
    forecast/observations/<target>.csv         live obs (date, raw, calibrated)
    forecast/ledger/<target>/<issue_id>/       one frozen projection
        projection.csv                         date, model, segment
        manifest.json                          fingerprint, config, fit stats
    forecast/scores/<target>.csv               per-issue scores
    forecast/REPORT.md                         summary table

Segments scored for each frozen projection
    holdout     the CV interval excluded from the committed fit
    post_data   after the last sample in the committed .dat (never fitted)
    post_issue  after the ledger issue time (strictly ex-ante)
Baselines, issued from the same origin with only data available then:
    persistence, damped persistence (AR1 on anomalies), monthly climatology.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import fnmatch
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import engines as E  # noqa: E402
import sources as S  # noqa: E402

LEDGER, OBS, SCORES = HERE / "ledger", HERE / "observations", HERE / "scores"


# ---------------------------------------------------------------------------
def now_year():
    t = dt.datetime.now(dt.timezone.utc)
    start = dt.datetime(t.year, 1, 1, tzinfo=dt.timezone.utc)
    return t.year + (t - start).total_seconds() / (365.2425 * 86400)


def git_rev():
    try:
        return subprocess.check_output(["git", "-C", str(E.REPO), "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:  # noqa: BLE001
        return None


def load_targets():
    cfg = yaml.safe_load((HERE / "targets.yml").read_text())
    defaults = cfg.get("defaults", {})
    out = {}
    for name, t in cfg["targets"].items():
        t = {**defaults, **t}
        if t.get("enabled", True) is False:
            continue
        if "glob" in t:  # expand to every matching fitted directory
            for d in sorted(E.FEB.iterdir()):
                if d.is_dir() and fnmatch.fnmatch(d.name, t["glob"]) and \
                        re.fullmatch(t.get("regex", ".*"), d.name) and \
                        (d / "lt.exe.p").exists() and (d / f"{d.name}.dat").exists():
                    tt = {k: v for k, v in t.items() if k not in ("glob", "regex")}
                    tt["dir"] = d.name
                    tt["group"] = name
                    if tt["source"]["kind"] in ("kaplan_box", "ersst_box"):
                        m = re.match(r"k([NS])(\d{3})_([EW])(\d{3})", d.name)
                        lat = int(m[2]) * (1 if m[1] == "N" else -1)
                        lon = int(m[4]) * (1 if m[3] == "E" else -1)
                        tt["source"] = {**tt["source"], "lat_c": lat, "lon_c": lon}
                    out[d.name] = tt
        else:
            t.setdefault("group", name)
            out[name] = t
    return out


def f9(v, n):
    for _ in range(n):
        r = v.copy(); r[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]; v = r
    return v


# ---------------------------------------------------------------------------
# observe
# ---------------------------------------------------------------------------
CAL_FILE = HERE / "observations" / "_calibration.json"


def observe(targets, only=None):
    OBS.mkdir(exist_ok=True)
    cal = json.loads(CAL_FILE.read_text()) if CAL_FILE.exists() else {}
    for name, t in targets.items():
        if only and name not in only and t["group"] not in only:
            continue
        try:
            d, v = S.fetch(t["source"])
        except Exception as e:  # noqa: BLE001
            print(f"[observe] {name}: FAILED {e}")
            continue
        if t.get("monthly_from_daily"):
            d, v = S.to_monthly(d, v)
        if t.get("deseasonalize"):
            # .dat for MSL sites has annual/semiannual removed upstream
            mon = np.round((d % 1) * 12).astype(int) % 12
            v = v.copy()
            for mo in range(12):
                s_ = mon == mo
                if s_.any():
                    v[s_] -= np.mean(v[s_])
        gain, offs, n_cal, r_cal = 1.0, 0.0, 0, float("nan")
        if t["engine"] == "lte" and t.get("calibrate", True):
            # map live obs onto the committed .dat units (anomaly base,
            # scaling) using overlap that predates the committed data end
            dat = np.loadtxt(E.FEB / t["dir"] / f"{t['dir']}.dat")
            key = {round(a, 3): b for a, b in dat}
            pairs = [(x, key[round(a, 3)]) for a, x in zip(d, v) if round(a, 3) in key]
            if len(pairs) >= 24:
                a = np.array(pairs)
                gain, offs = np.polyfit(a[:, 0], a[:, 1], 1)
                n_cal = len(pairs)
                r_cal = float(np.corrcoef(a[:, 0], a[:, 1])[0, 1])
        with open(OBS / f"{name}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([f"# fetched_utc={dt.datetime.now(dt.timezone.utc).isoformat()} "
                        f"calibration: value = {gain:.6g}*raw + {offs:.6g} (n={n_cal}, r={r_cal:.4f})"])
            w.writerow(["date", "raw", "value"])
            for a, x in zip(d, v):
                w.writerow([f"{a:.6f}", f"{x:.6g}", f"{gain * x + offs:.6g}"])
        cal[name] = dict(gain=gain, offset=offs, n=n_cal,
                         r=None if not np.isfinite(r_cal) else round(r_cal, 4))
        print(f"[observe] {name}: {len(d)} obs to {d[-1]:.3f} (cal n={n_cal}, r={r_cal:.3f})")


    CAL_FILE.write_text(json.dumps(cal, indent=1, sort_keys=True))


def read_obs(name):
    p = OBS / f"{name}.csv"
    if not p.exists():
        return None, None
    a = np.loadtxt(p, delimiter=",", skiprows=2, usecols=(0, 2), ndmin=2)
    return a[:, 0], a[:, 1]


# ---------------------------------------------------------------------------
# issue
# ---------------------------------------------------------------------------
def _write_issue(name, issue_id, proj, extra, forecast_only=False):
    out = LEDGER / name / issue_id
    if out.exists():
        return False
    if forecast_only:  # refit engines: keep the projected segment only
        k = proj["segment"] == "post_data"
        proj = {**proj, "dates": proj["dates"][k], "model": proj["model"][k],
                "segment": proj["segment"][k]}
    out.mkdir(parents=True)
    with open(out / "projection.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "model", "segment"])
        for a, m, s in zip(proj["dates"], proj["model"], proj["segment"]):
            w.writerow([f"{a:.6f}", f"{m:.8g}", s])
    meta = {**proj["meta"], **extra,
            "issued_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "issued_year": now_year(), "git_rev": git_rev()}
    (out / "manifest.json").write_text(json.dumps(meta, indent=1, default=str))
    return True


def issue(targets, only=None):
    for name, t in targets.items():
        if only and name not in only and t["group"] not in only:
            continue
        try:
            if t["engine"] == "lte":
                run_dir = E.FEB / t["dir"]
                if not (run_dir / "lt.exe.p").exists():
                    print(f"[issue] {name}: no committed fit at {run_dir} (pending)")
                    continue
                fp = E.lte_fingerprint(run_dir, t["dir"])
                proj = E.lte_project(run_dir, t["dir"], t.get("horizon_years", 5))
                new = _write_issue(name, fp[:12], proj, {"fingerprint": fp})
            elif t["engine"] in ("tidal_lod", "alias_cw"):
                ym = dt.date.today().strftime("%Y-%m")  # refit once per month
                if (LEDGER / name / ym).exists():
                    continue
                if t["engine"] == "tidal_lod":
                    d, v = S.fetch(t["source"])
                    proj = E.tidal_lod_project(d, v, d[-1], t.get("train_years", 12),
                                               t.get("horizon_years", 1),
                                               t.get("donor", "amo"))
                else:
                    d, x = S.fetch({**t["source"], "column": "x"})
                    _, y = S.fetch({**t["source"], "column": "y"})
                    proj = E.alias_cw_project(d, x, y, d[-1], t.get("train_years", 24),
                                              t.get("horizon_years", 2), t["component"])
                new = _write_issue(name, ym, proj, {}, forecast_only=True)
            else:
                continue
            print(f"[issue] {name}: {'NEW ' if new else 'unchanged '}"
                  f"r_train={proj['meta'].get('r_train')}")
        except Exception as e:  # noqa: BLE001
            print(f"[issue] {name}: FAILED {type(e).__name__}: {e}")


def backfill(targets, only=None, years=10.0, step_months=3):
    """Rolling-origin hindcasts for the refit engines (LOD, Chandler):
    at each origin only data <= origin are used, exactly as a live issue."""
    for name, t in targets.items():
        if t["engine"] not in ("tidal_lod", "alias_cw"):
            continue
        if only and name not in only and t["group"] not in only:
            continue
        if t["engine"] == "tidal_lod":
            d, v = S.fetch(t["source"])
        else:
            d, x = S.fetch({**t["source"], "column": "x"})
            _, y = S.fetch({**t["source"], "column": "y"})
        horizon = t.get("horizon_years", 1)
        last = d[-1] - horizon
        origins = np.arange(last - years, last + 1e-9, step_months / 12.0)
        n = 0
        for o in origins:
            y0 = int(o); m0 = int(round((o - y0) * 12)) + 1
            iid = f"hindcast-{y0:04d}-{min(m0, 12):02d}"
            if (LEDGER / name / iid).exists():
                continue
            if t["engine"] == "tidal_lod":
                proj = E.tidal_lod_project(d, v, o, t.get("train_years", 12), horizon,
                                           t.get("donor", "amo"))
            else:
                proj = E.alias_cw_project(d, x, y, o, t.get("train_years", 24), horizon,
                                          t["component"])
            n += _write_issue(name, iid, proj, {"hindcast_origin": float(o)},
                              forecast_only=True)
        print(f"[backfill] {name}: {n} new rolling-origin hindcasts")


# ---------------------------------------------------------------------------
# baselines + score
# ---------------------------------------------------------------------------
def baselines(d_hist, v_hist, d_eval, step):
    """Persistence, AR1-damped persistence, calendar-month climatology from
    history only (all inputs strictly at/before the origin)."""
    last_d, last = d_hist[-1], v_hist[-1]
    w = d_hist > last_d - 30
    mu = float(np.mean(v_hist[w]))
    a = v_hist[w] - mu
    phi = float(np.clip(np.corrcoef(a[:-1], a[1:])[0, 1], 0, 0.999)) if len(a) > 3 else 0.0
    lead = np.maximum(np.round((d_eval - last_d) / step), 1)
    persist = np.full(len(d_eval), last)
    ar1 = mu + phi ** lead * (last - mu)
    mon = lambda x: np.round((x % 1) * 12).astype(int) % 12  # noqa: E731
    clim_tab = {m: float(np.mean(v_hist[w][mon(d_hist[w]) == m]))
                for m in range(12) if np.any(mon(d_hist[w]) == m)}
    clim = np.array([clim_tab.get(m, mu) for m in mon(d_eval)])
    return dict(persistence=persist, ar1=ar1, climatology=clim), phi


def _metrics(model, obs, base):
    if len(obs) < 3:
        return None
    mse = lambda p: float(np.mean((p - obs) ** 2))  # noqa: E731
    out = dict(n=int(len(obs)), r=float(np.corrcoef(model, obs)[0, 1]),
               rmse=mse(model) ** 0.5)
    for k, b in base.items():
        out[f"rmse_{k}"] = mse(b) ** 0.5
        out[f"ss_{k}"] = 1 - mse(model) / mse(b) if mse(b) > 0 else None
    return out


def score(targets, only=None):
    SCORES.mkdir(exist_ok=True)
    for name, t in targets.items():
        if only and name not in only and t["group"] not in only:
            continue
        od, ov = read_obs(name)
        rows = []
        for issue_dir in sorted((LEDGER / name).glob("*")) if (LEDGER / name).exists() else []:
            meta = json.loads((issue_dir / "manifest.json").read_text())
            pr = np.genfromtxt(issue_dir / "projection.csv", delimiter=",", names=True,
                               dtype=None, encoding=None)
            pd_, pm = pr["date"].astype(float), pr["model"].astype(float)
            step = float(np.median(np.diff(pd_))) if len(pd_) > 1 else 1 / 12
            n_f9 = int(float(meta.get("env", {}).get("F9", 0))) if meta["engine"] == "lte" else 0
            if "hindcast_origin" in meta:
                meta["issued_year"] = float("inf")  # no ex-ante segment
            row = dict(target=name, issue=issue_dir.name, engine=meta["engine"],
                       issued=meta["issued_utc"], data_end=round(meta["data_end"], 4),
                       r_train=meta.get("r_train"), r_holdout=meta.get("r_holdout"))
            if od is not None and len(od) > 10:
                obs_f = f9(ov.copy(), n_f9)
                idx = np.clip(np.searchsorted(pd_, od - step / 2), 0, len(pd_) - 1)
                match = np.abs(pd_[idx] - od) < step / 2
                # the last n_f9 obs are not final under the centred smoother
                usable = match & (np.arange(len(od)) < len(od) - n_f9)
                for seg, origin in (("post_data", meta["data_end"]),
                                    ("post_issue", meta["issued_year"])):
                    m = usable & (od > origin)
                    hist = od <= origin
                    if hist.sum() < 24 or m.sum() < 3:
                        row[f"{seg}_n"] = int(m.sum())
                        continue
                    base, phi = baselines(od[hist], obs_f[hist], od[m], step)
                    res = _metrics(pm[idx[m]], obs_f[m], base)
                    for k, v in (res or {}).items():
                        row[f"{seg}_{k}"] = v
            rows.append(row)
        if rows:
            keys = sorted({k for r in rows for k in r},
                          key=lambda k: (k not in ("target", "issue", "engine", "issued",
                                                   "data_end"), k))
            with open(SCORES / f"{name}.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                for r in rows:
                    w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                                for k, v in r.items()})


# ---------------------------------------------------------------------------
def report(targets):
    lines = ["# GEM-LTE forecast ledger", "",
             f"Updated {dt.datetime.now(dt.timezone.utc):%Y-%m-%d %H:%M} UTC. "
             "Latest frozen projection per target. `ss_*` = 1 - MSE(model)/MSE(baseline); "
             "positive beats the baseline. post_data = after the committed .dat ends "
             "(never fitted); post_issue = after the ledger timestamp (strictly ex-ante).", "",
             "`obs cal r` = correlation of the live feed with the committed .dat over "
             "their overlap; LOW (< 0.9) means the live feed is not the same series "
             "and post-data scores for that target are not yet meaningful.", "",
             "| target | issue | data end | obs cal r | r holdout | post-data n | r | SS persist | SS AR1 | SS clim | post-issue n | r | SS persist |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    groups, hindcasts = {}, {}
    cal = json.loads(CAL_FILE.read_text()) if CAL_FILE.exists() else {}
    for name, t in targets.items():
        p = SCORES / f"{name}.csv"
        if not p.exists():
            continue
        with open(p) as f:
            rows = list(csv.DictReader(f))
        if not rows:
            continue
        live = [x for x in rows if not x["issue"].startswith("hindcast-")]
        hind = [x for x in rows if x["issue"].startswith("hindcast-")]
        if hind:
            hindcasts[name] = hind
        if not live:
            continue
        r = live[-1]
        c = cal.get(name, {}).get("r")
        r["cal"] = "" if c is None else (f"{c:.2f}" + (" LOW" if c < 0.9 else ""))
        groups.setdefault(t["group"], []).append(r)
        g = lambda k: r.get(k) or ""  # noqa: E731
        lines.append(f"| {name} | {g('issue')} | {g('data_end')} | {g('cal')} | {g('r_holdout')} | "
                     f"{g('post_data_n')} | {g('post_data_r')} | {g('post_data_ss_persistence')} | "
                     f"{g('post_data_ss_ar1')} | {g('post_data_ss_climatology')} | "
                     f"{g('post_issue_n')} | {g('post_issue_r')} | {g('post_issue_ss_persistence')} |")
    lines += ["", "## Group summary (median over members)", "",
              "| group | members | median post-data r | median SS AR1 | share SS AR1 > 0 |",
              "|---|---|---|---|---|"]
    for gname, rs in groups.items():
        def col(k):
            return np.array([float(r[k]) for r in rs if r.get(k) not in (None, "")])
        rr, ss = col("post_data_r"), col("post_data_ss_ar1")
        if len(rr):
            lines.append(f"| {gname} | {len(rs)} | {np.median(rr):.3f} | "
                         f"{np.median(ss):.3f} | {np.mean(ss > 0):.2f} |")
    if hindcasts:
        lines += ["", "## Rolling-origin hindcasts (refit engines)", "",
                  "Each origin uses only data available at that origin; scored over the "
                  "full horizon that followed.", "",
                  "| target | origins | median r | median SS persist | median SS AR1 | "
                  "median SS clim | share SS AR1 > 0 |", "|---|---|---|---|---|---|---|"]
        for name, rs in hindcasts.items():
            def hc(k):
                return np.array([float(r[k]) for r in rs if r.get(k) not in (None, "")])
            r_, sp, sa, sc = (hc("post_data_r"), hc("post_data_ss_persistence"),
                              hc("post_data_ss_ar1"), hc("post_data_ss_climatology"))
            if len(r_):
                lines.append(f"| {name} | {len(rs)} | {np.median(r_):.3f} | {np.median(sp):.3f} | "
                             f"{np.median(sa):.3f} | {np.median(sc):.3f} | {np.mean(sa > 0):.2f} |")
    (HERE / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(f"[report] {HERE / 'REPORT.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["all", "observe", "issue", "backfill", "score", "report",
                                    "list"])
    ap.add_argument("--years", type=float, default=10.0, help="backfill span")
    ap.add_argument("--step-months", type=int, default=3, help="backfill origin spacing")
    ap.add_argument("-t", "--targets", nargs="*", help="target names or group names")
    a = ap.parse_args()
    T = load_targets()
    if a.cmd == "list":
        for n, t in T.items():
            print(f"{n:16s} {t['engine']:10s} {t['group']:14s} {t['source']['kind']}")
        return
    if a.cmd in ("all", "observe"):
        observe(T, a.targets)
    if a.cmd in ("all", "issue"):
        issue(T, a.targets)
    if a.cmd == "backfill":
        backfill(T, a.targets, a.years, a.step_months)
        score(T, a.targets)
    if a.cmd in ("all", "score"):
        score(T, a.targets)
    if a.cmd in ("all", "report", "score", "backfill"):
        report(T)


if __name__ == "__main__":
    main()
