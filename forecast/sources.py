"""Live observation fetchers for the GEM-LTE forecast ledger.

Every fetcher returns (dates, values) as float numpy arrays, dates in the
project convention: year + (month-1)/12 for monthly series, and
year + (day_of_year-1)/365.2425 for daily series.  Missing-value sentinels
are dropped.  Fetchers are registered in SOURCES by name and referenced from
targets.yml as  source: {kind: <name>, ...kwargs}.
"""
from __future__ import annotations

import datetime as dt
import io
import time
import urllib.request
from pathlib import Path

import numpy as np

CACHE = Path(__file__).resolve().parent / ".cache"
UA = {"User-Agent": "gem-lte-forecast-ledger/1.0 (+https://github.com/pukpr/gem-lte-core)"}


def _get(url: str, binary: bool = False, tries: int = 3):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                b = r.read()
            return b if binary else b.decode("latin-1")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"fetch failed: {url}: {last}")


def _get_any(urls, binary=False):
    errs = []
    for u in urls:
        try:
            return _get(u, binary=binary)
        except RuntimeError as e:
            errs.append(str(e))
    raise RuntimeError(" | ".join(errs))


def _ym(y, m):
    return y + (m - 1) / 12.0


def _yd(y, doy):
    return y + (doy - 1) / 365.2425


def _clean(d, v, bad=lambda x: False):
    d, v = np.asarray(d, float), np.asarray(v, float)
    ok = np.isfinite(v) & ~np.vectorize(bad)(v) if len(v) else np.zeros(0, bool)
    return d[ok], v[ok]


def year_table(url, missing=(-999.9, 99.99, -99.99, -9.99), skip_nonnumeric=True,
               urls=None):
    """Rows 'YEAR v1 .. v12' (PSL/NCEI/CPC style).  Handles CPC's fused
    negative numbers ('-999.90-999.90')."""
    txt = _get_any(urls or [url])
    d, v = [], []
    for line in txt.splitlines():
        line = line.replace("-", " -").replace("e -", "e-").replace("E -", "E-")
        parts = line.split()
        if len(parts) < 13:
            continue
        try:
            y = int(parts[0])
            vals = [float(x) for x in parts[1:13]]
        except ValueError:
            continue
        if not 1700 < y < 2200:
            continue
        for m, x in enumerate(vals, 1):
            d.append(_ym(y, m)); v.append(x)
    return _clean(d, v, lambda x: any(abs(x - s) < 1e-6 for s in missing)
                  or abs(x) > 900)


def cpc_sstoi(url="https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices",
              column="NINO4_ANOM"):
    cols = {"NINO12": 2, "NINO12_ANOM": 3, "NINO3": 4, "NINO3_ANOM": 5,
            "NINO4": 6, "NINO4_ANOM": 7, "NINO34": 8, "NINO34_ANOM": 9}
    j = cols[column]
    d, v = [], []
    for line in _get(url).splitlines()[1:]:
        p = line.split()
        if len(p) >= 10:
            d.append(_ym(int(p[0]), int(p[1]))); v.append(float(p[j]))
    return _clean(d, v)


def longpaddock_soi_monthly(url="https://data.longpaddock.qld.gov.au/SeasonalClimateOutlook/"
                            "SouthernOscillationIndex/SOIDataFiles/MonthlySOIPhase1887-1989Base.txt"):
    d, v = [], []
    for line in _get(url).splitlines():
        p = line.split()
        try:
            y, m, s = int(p[0]), int(p[1]), float(p[2])
        except (ValueError, IndexError):
            continue
        d.append(_ym(y, m)); v.append(s)
    return _clean(d, v)


def longpaddock_soi_daily(url="https://data.longpaddock.qld.gov.au/SeasonalClimateOutlook/"
                          "SouthernOscillationIndex/SOIDataFiles/DailySOI1933-1992Base.txt",
                          column="SOI"):
    """Daily Tahiti/Darwin pressures and SOI.  column in {SOI, Tahiti, Darwin}."""
    j = {"Tahiti": 2, "Darwin": 3, "SOI": 4}[column]
    d, v = [], []
    for line in _get(url).splitlines():
        p = line.split()
        try:
            y, doy, x = int(p[0]), int(p[1]), float(p[j])
        except (ValueError, IndexError):
            continue
        d.append(_yd(y, doy)); v.append(x)
    return _clean(d, v, lambda x: abs(x) > 2000)


def psmsl_stations_in_box(lat0, lat1, lon0, lon1,
                          url="https://psmsl.org/data/obtaining/rlr.monthly.data/filelist.txt"):
    ids = []
    for line in _get(url).splitlines():
        p = [x.strip() for x in line.split(";")]
        if len(p) >= 3:
            la, lo = float(p[1]), float(p[2])
            if lat0 <= la <= lat1 and lon0 <= lo <= lon1:
                ids.append(int(p[0]))
    return ids


def psmsl_rlr(stations=None, bbox=None,
              url="https://psmsl.org/data/obtaining/rlr.monthly.data/{id}.rlrdata",
              min_stations=1):
    """Monthly RLR mean sea level (mm); a list of station ids is averaged
    as anomalies (each station's own mean removed) then re-offset by the
    pooled mean, matching a composite index like experiments/Feb2026/baltic."""
    if bbox is not None:
        stations = psmsl_stations_in_box(*bbox)
    if isinstance(stations, (int, str)):
        stations = [stations]
    series = {}
    for sid in stations:
        try:
            txt = _get(url.format(id=sid))
        except RuntimeError:
            continue
        for line in txt.splitlines():
            p = [x.strip() for x in line.split(";")]
            if len(p) < 2:
                continue
            t, h = float(p[0]), float(p[1])
            if h <= -99999:
                continue
            y = int(t); m = int(round((t - y) * 12 + 0.5))
            key = round(_ym(y, m), 6)
            series.setdefault(sid, {})[key] = h
    if not series:
        raise RuntimeError("no PSMSL stations fetched")
    means = {s: np.mean(list(v.values())) for s, v in series.items()}
    pooled = float(np.mean(list(means.values())))
    allk = sorted(set().union(*[v.keys() for v in series.values()]))
    d, v = [], []
    for k in allk:
        xs = [series[s][k] - means[s] for s in series if k in series[s]]
        if len(xs) >= min_stations:
            d.append(k); v.append(pooled + float(np.mean(xs)))
    return _clean(d, v)


def kaplan_box(lat_c, lon_c, half=10.0,
               url="https://downloads.psl.noaa.gov/Datasets/kaplan_sst/sst.mean.anom.nc"):
    """cos(lat)-weighted box mean of Kaplan Extended SST v2 anomalies, using
    the same 20-degree box convention as build_k_sst.extract_regions()."""
    import netCDF4 as nc  # local import: only needed for the quads
    CACHE.mkdir(exist_ok=True)
    f = CACHE / "kaplan_sst.mean.anom.nc"
    if not f.exists() or time.time() - f.stat().st_mtime > 3 * 86400:
        f.write_bytes(_get(url, binary=True))
    ds = nc.Dataset(f)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    tv = ds.variables["time"]
    cal = nc.num2date(tv[:], tv.units)
    dates = np.array([c.year + (c.month - 1) / 12.0 for c in cal])
    lon360 = lon_c % 360.0
    la = np.nonzero((lat >= lat_c - half) & (lat < lat_c + half))[0]
    lo = np.nonzero((lon >= lon360 - half) & (lon < lon360 + half))[0]
    box = ds.variables["sst"][:, la[0]:la[-1] + 1, lo[0]:lo[-1] + 1]
    w = np.cos(np.deg2rad(lat[la]))[:, None] * np.ones((1, len(lo)))
    m = ~np.ma.getmaskarray(box)
    num = np.sum(np.where(m, box.filled(0.0) * w, 0.0), axis=(1, 2))
    den = np.sum(np.where(m, w, 0.0), axis=(1, 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        v = num / den
    return _clean(dates, v)


def ersst_box(lat_c, lon_c, half=10.0, base=(1951, 1981),
              url="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"):
    """Live substitute for the Kaplan quads (Kaplan v2 ends Jan 2023):
    ERSSTv5 cos(lat)-weighted box mean, anomaly vs the Kaplan 1951-1980
    calendar-month base; same 20-degree box convention."""
    import netCDF4 as nc
    CACHE.mkdir(exist_ok=True)
    f = CACHE / "ersst.v5.sst.mnmean.nc"
    if not f.exists() or time.time() - f.stat().st_mtime > 3 * 86400:
        f.write_bytes(_get(url, binary=True))
    ds = nc.Dataset(f)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    tv = ds.variables["time"]
    cal = nc.num2date(tv[:], tv.units)
    years = np.array([c.year for c in cal]); months = np.array([c.month for c in cal])
    dates = years + (months - 1) / 12.0
    lon360 = lon_c % 360.0
    la = np.nonzero((lat >= lat_c - half) & (lat < lat_c + half))[0]
    lo = np.nonzero((lon >= lon360 - half) & (lon < lon360 + half))[0]
    box = ds.variables["sst"][:, la.min():la.max() + 1, lo.min():lo.max() + 1]
    w = np.cos(np.deg2rad(lat[la.min():la.max() + 1]))[:, None] * np.ones((1, box.shape[2]))
    m = ~np.ma.getmaskarray(box)
    num = np.sum(np.where(m, box.filled(0.0) * w, 0.0), axis=(1, 2))
    den = np.sum(np.where(m, w, 0.0), axis=(1, 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        v = num / den
    for mo in range(1, 13):
        s = (months == mo); b = s & (years >= base[0]) & (years < base[1])
        v[s] = v[s] - np.nanmean(v[b])
    return _clean(dates, v)


_EOP_CACHE = {}


def iers_eop(column="LOD", url="https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now"):
    """Daily IERS EOP 20 C04.  column in {x, y, LOD, UT1-UTC}.  LOD in ms,
    x/y in arcsec."""
    if url not in _EOP_CACHE:
        rows = []
        for line in _get(url).splitlines():
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 13:
                continue
            y, mo, da = int(p[0]), int(p[1]), int(p[2])
            doy = dt.date(y, mo, da).timetuple().tm_yday
            rows.append((_yd(y, doy), float(p[5]), float(p[6]), float(p[7]),
                         float(p[12]) * 1000.0))
        _EOP_CACHE[url] = np.array(rows)
    a = _EOP_CACHE[url]
    j = {"x": 1, "y": 2, "UT1-UTC": 3, "LOD": 4}[column]
    return a[:, 0].copy(), a[:, j].copy()


def kit_qbo(level=30, url="https://www.atmohub.kit.edu/data/qbo.dat"):
    """FU Berlin / KIT Canton-Gan-Singapore monthly zonal wind, m/s."""
    cols = {70: (10, 18), 50: (18, 25), 40: (25, 32), 30: (32, 39),
            20: (39, 46), 15: (46, 53), 10: (53, 60)}
    a, b = cols[int(level)]
    d, v = [], []
    for line in _get(url).splitlines():
        if not line[:5].isdigit():
            continue
        try:
            yy, mm = int(line[6:8]), int(line[8:10])
            x = float(line[a:b][:-1])
        except ValueError:
            continue
        d.append(_ym(1900 + yy if yy >= 50 else 2000 + yy, mm)); v.append(x / 10.0)
    return _clean(d, v)


SOURCES = {
    "kit_qbo": kit_qbo,
    "year_table": year_table,
    "cpc_sstoi": cpc_sstoi,
    "longpaddock_soi_monthly": longpaddock_soi_monthly,
    "longpaddock_soi_daily": longpaddock_soi_daily,
    "psmsl_rlr": psmsl_rlr,
    "kaplan_box": kaplan_box,
    "ersst_box": ersst_box,
    "iers_eop": iers_eop,
}


def fetch(spec: dict):
    spec = dict(spec)
    kind = spec.pop("kind")
    return SOURCES[kind](**spec)


def to_monthly(d, v):
    """Calendar-month means of a daily series, project date convention."""
    y = np.floor(d).astype(int)
    doy = np.round((d - y) * 365.2425).astype(int) + 1
    mon = np.array([(dt.date(int(a), 1, 1) + dt.timedelta(int(b) - 1)).month
                    for a, b in zip(y, doy)])
    key = y * 12 + (mon - 1)
    out_d, out_v = [], []
    for k in np.unique(key):
        s = key == k
        out_d.append(k // 12 + (k % 12) / 12.0); out_v.append(float(np.mean(v[s])))
    return np.array(out_d), np.array(out_v)
