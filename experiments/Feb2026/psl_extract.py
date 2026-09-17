#!/usr/bin/env python3
"""psl_extract.py — create an analysis subdirectory from a NOAA PSL gridded
timeseries extraction (Kaplan SST V2 / NCEP Reanalysis / etc.), ready for
lte_run.py + winding_rank.py / winding_scalogram.py.

Reproduces the saved-page workflow used for kap10-10-20-30 (the form at
https://psl.noaa.gov/cgi-bin/data/timeseries/timeseries1.pl) directly and
converts the year-x-month HTML table into the bare "<year> <value>" .dat
layout the Ada tool expects. Verified: re-fetching the kap10-10-20-30 box
(Kaplan SST, 20..30N, -10..10E) reproduces the original .dat at r=1.000000,
max diff 0.0000 — the month convention is start-of-month (year + m/12).

    # Kaplan SST V2 (ntype=4), Atlantic Nino region:
    ./psl_extract.py atl3n3s --ntype 4 --lat 3 -3 --lon -20 0
    # NCEP Reanalysis SLP (ntype=1):
    ./psl_extract.py nao_slp --ntype 1 --var "Sea Level Pressure" \
        --lat 37 70 --lon -170 -10
    # then:
    ./lte_run.py atl3n3s            # (or the GUI) to fit the manifold
    ./winding_rank.py --index atl3n3s

--lon follows the form's degrees-east convention; 0..360 gives a zonal mean.
"""
import argparse
import gzip
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = "https://psl.noaa.gov/cgi-bin/data/timeseries/timeseries.pl"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
          "Aug", "Sep", "Oct", "Nov", "Dec"]

DATASETS = {1: "NCEP Reanalysis", 2: "NOAA/NCEI Climate Division",
            4: "Kaplan SST V2", 6: "U of Delaware Precipitation"}


def fetch(ntype, lat1, lat2, lon1, lon2, var=None, level=None,
          seasonal=False, mon1="Jan", mon2="Jan", area=False,
          retries=3, pause=2.0):
    params = {"ntype": ntype, "lat1": lat1, "lat2": lat2,
              "lon1": lon1, "lon2": lon2,
              "iseas": 1 if seasonal else 0,
              "iarea": 1 if area else 0,
              "typeout": 1, "mon1": mon1, "mon2": mon2,
              "Submit": "Create Timeseries"}
    if var:
        params["var"] = var
    if level:
        params["level"] = level
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent":
                                               "Mozilla/5.0 (lte extract)"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read()
                if raw[:2] == b"\x1f\x8b":
                    raw = gzip.decompress(raw)
                return raw.decode(errors="replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(pause * (attempt + 1))
    raise SystemExit(f"fetch failed after {retries} tries: {last}")


def parse_table(html):
    """-> (t, x) start-of-month decimal years, values. Handles both monthly
    (year + 12 cols) and seasonal (year + 1 col) layouts."""
    pre = re.search(r"<pre>(.*?)</pre>", html, re.S)
    body = pre.group(1) if pre else html
    body = body.replace("&nbsp;", " ")
    lines = body.splitlines()
    grid = tuple(l.strip() for l in lines
                 if l.startswith(("Latitude Range used",
                                  "Longitude Range used")))
    rows = []
    for line in lines:
        p = line.split()
        if len(p) >= 13 and re.match(r"^\d{4}$", p[0]):
            vals = [float(v) if v != "---" else float("nan")
                    for v in p[1:13]]
            vals = [np.nan if v <= -999.0 else v for v in vals]
            rows.append((int(p[0]), vals))
        elif len(p) == 2 and re.match(r"^\d{4}$", p[0]):
            continue  # header row "first_year last_year"
    if not rows:
        # single-column seasonal layout: year + one value
        for line in lines:
            p = line.split()
            if len(p) == 2 and re.match(r"^\d{4}$", p[0]):
                try:
                    rows.append((int(p[0]), [float(p[1])]))
                except ValueError:
                    pass
    if not rows:
        raise SystemExit("no numeric year-rows found in response — page "
                         "format may differ; inspect with --dump")
    ncol = 12 if len(rows[0][1]) >= 12 else len(rows[0][1])
    ts, xs = [], []
    for y, vals in rows:
        for m, v in enumerate(vals[:ncol]):
            ts.append(y + m / 12.0)
            xs.append(v)
    ts, xs = np.array(ts), np.array(xs)
    keep = np.isfinite(xs)
    return ts[keep], xs[keep], ncol, grid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=
                                 argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="subdirectory to create under "
                                 "experiments/Feb2026 (e.g. atl3n3s)")
    ap.add_argument("--ntype", type=int, default=4,
                    help="PSL dataset: 1=NCEP Reanalysis, 2=NCEI Climate "
                         "Division, 4=Kaplan SST V2, 6=UDel Precip "
                         "(default 4)")
    ap.add_argument("--lat", type=float, nargs=2, required=True,
                    metavar=("FROM", "TO"), help="N to S, e.g. 20 30")
    ap.add_argument("--lon", type=float, nargs=2, required=True,
                    metavar=("FROM", "TO"),
                    help="W to E, degrees east; 0 360 = zonal mean")
    ap.add_argument("--var", default=None,
                    help="NCEP variable name (only for --ntype 1), e.g. "
                         "'Sea Level Pressure'")
    ap.add_argument("--level", default=None,
                    help="NCEP pressure level, e.g. 1000 (only --ntype 1)")
    ap.add_argument("--seasonal", action="store_true")
    ap.add_argument("--months", nargs=2, default=["Jan", "Jan"],
                    metavar=("M1", "M2"), help="season month range")
    ap.add_argument("--area-weight", action="store_true")
    ap.add_argument("--dump", default=None,
                    help="also save the raw HTML page here")
    args = ap.parse_args()

    outdir = ROOT / args.name
    if outdir.exists():
        raise SystemExit(f"{outdir} already exists — refusing to overwrite")
    html = fetch(args.ntype, args.lat[0], args.lat[1],
                 args.lon[0], args.lon[1], var=args.var, level=args.level,
                 seasonal=args.seasonal, mon1=args.months[0],
                 mon2=args.months[1], area=args.area_weight)
    if args.dump:
        Path(args.dump).write_text(html)
    t, x, ncol, grid = parse_table(html)
    print(f"{DATASETS.get(args.ntype, args.ntype)}  "
          f"lat {args.lat[0]}..{args.lat[1]} lon {args.lon[0]}..{args.lon[1]}"
          f"{' var=' + args.var if args.var else ''}: "
          f"{len(t)} points, {t[0]:.0f}..{t[-1]:.2f}, "
          f"{'monthly' if ncol == 12 else f'seasonal ({ncol} col)'}")
    for g in grid:
        print(f"  {g}")

    outdir.mkdir(parents=True)
    dat = outdir / f"{args.name}.dat"
    dat.write_text("\n".join(f"{y:.6f} {v:.3f}" for y, v in zip(t, x)) + "\n")
    print(f"  wrote {dat}")
    txt = outdir / f"{args.name}.txt"
    yr0, yr1 = int(t[0]), int(t[-1])
    with txt.open("w") as f:
        for y in range(yr0, yr1 + 1):
            m = (t >= y) & (t < y + 1)
            if m.sum() == 12:
                f.write(f"{y:5d}" + "".join(f"{v:9.3f}" for v in x[m]) + "\n")
    print(f"  wrote {txt} (year x month layout, matches the saved-page "
          f"format)")
    print(f"\nNext:  ./lte_run.py {args.name}   "
          f"then  ./winding_rank.py --index {args.name}")


if __name__ == "__main__":
    main()
