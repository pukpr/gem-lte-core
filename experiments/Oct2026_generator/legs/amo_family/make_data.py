#!/usr/bin/env python3
"""AMO-family legs (2026-11-21), collected like nao_family: every series monthly, calendar-month means removed (each over
its own record), NOT smoothed (amo.dat, the Feb2026 AMO, is unsmoothed monthly), dates as year + (month-1)/12, missing
values dropped (the HadSST series have gaps: lt.exe skips missing rows rather than filling them). Scale is left as
published (CC is scale-free). Left out as duplicates of amo_kaplan.dat: lteMod/run/amo_1870_extend.dat (CC 1.000) and
PSL's detrended 1948- version raw/amon.us.data (1.000).
  amo.dat              experiments/Feb2026/amo/amo.dat, copied: an earlier PSL Kaplan download, x10 (CC 0.993 with
                       amo_kaplan.dat)
  amo_kaplan.dat       NOAA PSL AMO unsmoothed, Kaplan SST V2, 0-70N 80W-0, 1856- (raw/amon.us.long.data)
  amo_kaplan_sm.dat    NOAA PSL AMO 121-month smoothed, Kaplan SST V2, 1856- (raw/amon.sm.long.data)
  amo_ersst_vo.dat     KNMI Climate Explorer: ERSST v5 25-60N 7-75W minus regression on global mean T
                       (van Oldenborgh et al 2009), 1880- (raw/iamo_ersst.dat)
  amo_hadsst_vo.dat    KNMI: the same from HadSST 4.2, 1850- (raw/iamo_hadsst.dat)
  amo_ersst_ts.dat     KNMI: ERSST v5 EQ-60N 0-80W minus 60S-60N (Trenberth and Shea 2006), 1854- (raw/iamo_ersst_ts.dat)
  amo_hadsst_ts.dat    KNMI: the same from HadSST 4.2 (raw/iamo_hadsst_ts.dat)
  amo_kaplan_ts.dat    Trenberth-Shea from the Kaplan SST V2 grid (Feb2026/k_sst/sst.mean.anom.nc), area-weighted
  spg_kaplan.dat       subpolar gyre SST, 45-65N 60W-10W, Kaplan SST V2 grid, area-weighted
  amo_ucar.dat         experiments/Feb2026/amo/amo_ucar.dat (NCAR, provenance not recorded), 1880-2021
  natl_psl.dat         experiments/Feb2026/amo/amo_psl.dat: absolute North Atlantic SST (deg C, provenance not recorded),
                       here as anomalies
  tna.dat              experiments/Feb2026/tna/tna.dat, copied (tropical North Atlantic)
  tna_psl.dat          NOAA PSL TNA, HadISST 5.5-23.5N 15-57.5W, 1948- (raw/tna.data)
  tsa_psl.dat          NOAA PSL TSA, tropical South Atlantic, 1948- (raw/tsa.data): the Atlantic counterpart, not an AMO
raw/ holds the downloads as fetched (2026-11-21).   python3 -I make_data.py"""
import os, shutil, numpy as np
H = os.path.expanduser("~/eval"); HERE = os.path.dirname(os.path.abspath(__file__))
RAW, OUT = os.path.join(HERE, "raw"), os.path.join(HERE, "data"); os.makedirs(OUT, exist_ok=True)
FEB = os.path.join(H, "gem-lte-core/experiments/Feb2026")
def anomalies(d):
    ks = np.array(sorted(d)); v = np.array([d[k] for k in ks])
    for m in range(12): v[ks % 12 == m] -= v[ks % 12 == m].mean()
    return dict(zip(ks, v))
def dat(f, scale=1.0):  # two-column year.frac value
    d = np.loadtxt(f, usecols=(0, 1)); return {int(round(t * 12)): v * scale for t, v in d}
def table(f, bad):  # year + 12 monthly values per row (PSL .data with a "y0 y1" first line, KNMI with # comments)
    out = {}
    for l in open(f):
        x = l.split()
        if l.startswith("#") or len(x) != 13: continue
        try: y = int(x[0]); v = [float(s) for s in x[1:]]
        except ValueError: continue
        for m, s in enumerate(v):
            if abs(s - bad) > 1e-3 and abs(s) < 90: out[y * 12 + m] = s
    return out
def kaplan(boxes):  # boxes: [(weight, lat0, lat1, lon0, lon1)], lon in 0..360; area-weighted box means, combined linearly
    import netCDF4
    d = netCDF4.Dataset(os.path.join(FEB, "k_sst/sst.mean.anom.nc"))
    lat, lon, t = d["lat"][:], d["lon"][:], d["time"][:]; sst = d["sst"][:]
    tot = 0.0
    for w, a0, a1, o0, o1 in boxes:
        la = (lat >= a0) & (lat <= a1); lo = ((lon >= o0) & (lon <= o1)) if o0 <= o1 else ((lon >= o0) | (lon <= o1))
        x = sst[:, la][:, :, lo]; cw = np.cos(np.radians(lat[la]))[None, :, None] * np.ones(x.shape)
        cw = np.ma.masked_array(cw, mask=np.ma.getmaskarray(x)); tot = tot + w * (x * cw).sum(axis=(1, 2)) / cw.sum(axis=(1, 2))
    import datetime as dt
    out = {}
    for ti, v in zip(t, tot):
        if np.ma.is_masked(v): continue
        day = dt.datetime(1800, 1, 1) + dt.timedelta(days=float(ti)); out[day.year * 12 + day.month - 1] = float(v)
    return out
def write(name, d, note):
    d = anomalies(d); ks = sorted(d)
    with open(os.path.join(OUT, name), "w") as fh:
        for k in ks: fh.write(f"{k // 12 + (k % 12) / 12:.6f}\t{d[k]:.6f}\n")
    print(f"{name:20s} {ks[0] // 12}-{ks[-1] // 12}  {len(ks):5d} months  {note}")
shutil.copy(os.path.join(FEB, "amo/amo.dat"), os.path.join(OUT, "amo.dat")); print("amo.dat              copied")
write("amo_kaplan.dat", table(os.path.join(RAW, "amon.us.long.data"), -99.99), "PSL Kaplan unsmoothed")
write("amo_kaplan_sm.dat", table(os.path.join(RAW, "amon.sm.long.data"), -99.99), "PSL Kaplan 121-month smoothed")
write("amo_ersst_vo.dat", table(os.path.join(RAW, "iamo_ersst.dat"), -999.9), "ERSST v5, minus global-T regression")
write("amo_hadsst_vo.dat", table(os.path.join(RAW, "iamo_hadsst.dat"), -999.9), "HadSST 4.2, minus global-T regression")
write("amo_ersst_ts.dat", table(os.path.join(RAW, "iamo_ersst_ts.dat"), -999.9), "ERSST v5, Trenberth-Shea")
write("amo_hadsst_ts.dat", table(os.path.join(RAW, "iamo_hadsst_ts.dat"), -999.9), "HadSST 4.2, Trenberth-Shea")
write("amo_kaplan_ts.dat", kaplan([(1.0, 0, 60, 280, 360), (-1.0, -60, 60, 0, 360)]), "Kaplan grid, Trenberth-Shea")
write("spg_kaplan.dat", kaplan([(1.0, 45, 65, 300, 350)]), "Kaplan grid, subpolar gyre 45-65N 60-10W")
write("amo_ucar.dat", dat(os.path.join(FEB, "amo/amo_ucar.dat")), "NCAR copy (provenance not recorded)")
write("natl_psl.dat", dat(os.path.join(FEB, "amo/amo_psl.dat")), "absolute N Atlantic SST as anomalies")
shutil.copy(os.path.join(FEB, "tna/tna.dat"), os.path.join(OUT, "tna.dat")); print("tna.dat              copied")
write("tna_psl.dat", table(os.path.join(RAW, "tna.data"), -99.99), "PSL TNA (HadISST)")
write("tsa_psl.dat", table(os.path.join(RAW, "tsa.data"), -99.99), "PSL TSA (HadISST), Atlantic counterpart")
