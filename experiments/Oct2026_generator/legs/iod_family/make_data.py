#!/usr/bin/env python3
"""IOD-family legs (2026-11-21), collected like nao_family / amo_family: every series monthly, calendar-month means removed
(each over its own record), then a centred 3-month boxcar (full windows only; 1 month dropped at each end and around
gaps), dates as year + (month-1)/12. That is iode.dat's footing: Feb2026/iode/iode.dat is the 3-month boxcar of
lteMod/run/iod_e.dat (CC 0.997; 0.977 unsmoothed), which is PSL's HadISST SETIO (0.999). Scale is left as published.
Regions: WTIO 10S-10N 50-70E (west pole), SETIO 10S-0 90-110E (east pole), DMI = WTIO - SETIO.
  iode.dat           experiments/Feb2026/iode/iode.dat, copied (SETIO, HadISST, 1880-2022)
  iodw.dat           experiments/Feb2026/iodw/iodw.dat, copied (WTIO, HadISST, 1880-2022)
  setio_had.dat      NOAA PSL SETIO, HadISST 1.1, 1870- (raw/dmieast.had.long.data): iode.dat's source, longer
  wtio_had.dat       NOAA PSL WTIO, HadISST 1.1, 1870- (raw/dmiwest.had.long.data): iodw.dat's source, longer
  dmi_had.dat        NOAA PSL DMI, HadISST 1.1, 1870- (raw/dmi.had.long.data)
  setio_ersst.dat    KNMI Climate Explorer SETIO, ERSST v5 (in situ only), 1854- (raw/iseio_ersst.dat)
  wtio_ersst.dat     KNMI WTIO, ERSST v5, 1854- (raw/iwtio_ersst.dat)
  dmi_ersst.dat      KNMI DMI, ERSST v5, 1854- (raw/idmi_ersst.dat)
  setio_kaplan.dat   SETIO from the Kaplan SST V2 grid (Feb2026/k_sst/sst.mean.anom.nc), area-weighted, 1856-2023
  wtio_kaplan.dat    WTIO from the Kaplan grid
  dmi_kaplan.dat     their difference
  iod.dat            experiments/Feb2026/iod/iod.dat (a DMI, 1900-; provenance not recorded; CC 0.97 with KNMI's ERSST DMI)
Left out as duplicates of dmi_had.dat: GEM-LTE/experiments/data/dmi.dat and lteMod/run/iod_knmi_1870.dat (CC 0.98, the
same series); Feb2026/dmi_west_minus_east.dat = iodw - iode (derivable). Left out as unidentified: GeoEnergyMath/iod/iod.dat
(CC <= 0.5 with every other series). raw/ holds the downloads as fetched (2026-11-21).   python3 -I make_data.py"""
import os, shutil, numpy as np
H = os.path.expanduser("~/eval"); HERE = os.path.dirname(os.path.abspath(__file__))
RAW, OUT = os.path.join(HERE, "raw"), os.path.join(HERE, "data"); os.makedirs(OUT, exist_ok=True)
FEB = os.path.join(H, "gem-lte-core/experiments/Feb2026")
def anomalies(d):
    ks = np.array(sorted(d)); v = np.array([d[k] for k in ks])
    for m in range(12): v[ks % 12 == m] -= v[ks % 12 == m].mean()
    return dict(zip(ks, v))
def box3(d):
    return {k: (d[k - 1] + d[k] + d[k + 1]) / 3 for k in d if k - 1 in d and k + 1 in d}
def dat(f):  # two-column year.frac value (KNMI's fractional dates are year + (month-1)/12 to 1e-4)
    out = {}
    for l in open(f):
        if l.startswith("#"): continue
        x = l.split()
        if len(x) < 2: continue
        v = float(x[1])
        if abs(v) < 90: out[int(round(float(x[0]) * 12))] = v
    return out
def table(f):  # year + 12 monthly values per row (PSL .data with a "y0 y1" first line, KNMI with # comments)
    out = {}
    for l in open(f):
        x = l.split()
        if l.startswith("#") or len(x) != 13: continue
        try: y = int(x[0]); v = [float(s) for s in x[1:]]
        except ValueError: continue
        for m, s in enumerate(v):
            if abs(s) < 90: out[y * 12 + m] = s
    return out
def kaplan(boxes):  # boxes: [(weight, lat0, lat1, lon0, lon1)], lon 0..360; area-weighted box means, combined linearly
    import netCDF4, datetime as dt
    d = netCDF4.Dataset(os.path.join(FEB, "k_sst/sst.mean.anom.nc"))
    lat, lon, t = d["lat"][:], d["lon"][:], d["time"][:]; sst = d["sst"][:]
    tot = 0.0
    for w, a0, a1, o0, o1 in boxes:
        la = (lat >= a0) & (lat <= a1); lo = (lon >= o0) & (lon <= o1)
        x = sst[:, la][:, :, lo]; cw = np.cos(np.radians(lat[la]))[None, :, None] * np.ones(x.shape)
        cw = np.ma.masked_array(cw, mask=np.ma.getmaskarray(x)); tot = tot + w * (x * cw).sum(axis=(1, 2)) / cw.sum(axis=(1, 2))
    out = {}
    for ti, v in zip(t, tot):
        if np.ma.is_masked(v): continue
        day = dt.datetime(1800, 1, 1) + dt.timedelta(days=float(ti)); out[day.year * 12 + day.month - 1] = float(v)
    return out
def write(name, d, note):
    d = box3(anomalies(d)); ks = sorted(d)
    with open(os.path.join(OUT, name), "w") as fh:
        for k in ks: fh.write(f"{k // 12 + (k % 12) / 12:.6f}\t{d[k]:.6f}\n")
    print(f"{name:18s} {ks[0] // 12}-{ks[-1] // 12}  {len(ks):5d} months  {note}")
for f in ("iode", "iodw"):
    shutil.copy(os.path.join(FEB, f, f + ".dat"), os.path.join(OUT, f + ".dat")); print(f"{f}.dat           copied")
write("setio_had.dat", table(os.path.join(RAW, "dmieast.had.long.data")), "PSL HadISST east pole")
write("wtio_had.dat", table(os.path.join(RAW, "dmiwest.had.long.data")), "PSL HadISST west pole")
write("dmi_had.dat", table(os.path.join(RAW, "dmi.had.long.data")), "PSL HadISST DMI")
write("setio_ersst.dat", dat(os.path.join(RAW, "iseio_ersst.dat")), "KNMI ERSST v5 east pole")
write("wtio_ersst.dat", dat(os.path.join(RAW, "iwtio_ersst.dat")), "KNMI ERSST v5 west pole")
write("dmi_ersst.dat", table(os.path.join(RAW, "idmi_ersst.dat")), "KNMI ERSST v5 DMI")
E, W = [(1.0, -10, 0, 90, 110)], [(1.0, -10, 10, 50, 70)]
write("setio_kaplan.dat", kaplan(E), "Kaplan grid east pole")
write("wtio_kaplan.dat", kaplan(W), "Kaplan grid west pole")
write("dmi_kaplan.dat", kaplan(W + [(-1.0,) + E[0][1:]]), "Kaplan grid DMI")
write("iod.dat", dat(os.path.join(FEB, "iod/iod.dat")), "Feb2026 iod.dat (DMI, provenance not recorded)")
