"""ORAS5 depth of the 20 C isotherm (so20chgt, monthly, native ORCA025 grid) over the tropical Pacific (2026-10-21).
For each year: retrieve the 12 monthly global fields from the Copernicus CDS (reanalysis-oras5, single_level;
'consolidated' up to 2014, 'operational' from 2015), average the native cells into 1 x 1 degree boxes over
15S-15N, 120E-285E (box centres at half degrees), save d20_<year>.npz (time, lat, lon, d20 [m], NaN = land/no data),
then delete the raw zip and netCDF. Resumable: years with an npz are skipped. A year whose request fails is logged and
retried on the next run.
Run with the QBO environment (cdsapi >= 0.7.2):  nohup ../QBO/.venv/bin/python fetch_d20.py 1958 2025 &
"""
import glob, os, shutil, sys, time, zipfile
import cdsapi
import netCDF4 as nc
import numpy as np

LAT = np.arange(-14.5, 15.0, 1.0); LON = np.arange(120.5, 285.0, 1.0)


def boxes(path):
    ds = nc.Dataset(path); lat = np.asarray(ds.variables["nav_lat"][:]); lon = np.asarray(ds.variables["nav_lon"][:]) % 360
    d = np.ma.filled(np.ma.asarray(ds.variables["so20chgt"][:]).squeeze().astype(float), np.nan)
    tv = ds.variables["time_counter"]; t = nc.num2date(tv[:], tv.units)[0]
    sel = (np.abs(lat) < 15) & (lon >= 120) & (lon < 285)
    i = np.floor(lat[sel] + 15).astype(int); j = np.floor(lon[sel] - 120).astype(int); v = d[sel]
    ok = np.isfinite(v) & (v > 0)
    s = np.zeros((len(LAT), len(LON))); c = np.zeros_like(s)
    np.add.at(s, (i[ok], j[ok]), v[ok]); np.add.at(c, (i[ok], j[ok]), 1)
    with np.errstate(invalid="ignore"):
        return t.year + (t.month - 1) / 12.0, np.where(c > 0, s / c, np.nan)


def year(y, cl):
    out = f"d20_{y}.npz"
    if os.path.exists(out):
        return
    z, tmp = f"raw_{y}.zip", f"raw_{y}"
    cl.retrieve("reanalysis-oras5", {"product_type": ["consolidated" if y <= 2014 else "operational"], "vertical_resolution": "single_level",
                                     "variable": ["depth_of_20_c_isotherm"], "year": [str(y)], "month": [f"{m:02d}" for m in range(1, 13)]}, z)
    zipfile.ZipFile(z).extractall(tmp)
    rows = sorted(boxes(f) for f in glob.glob(f"{tmp}/*.nc"))
    np.savez_compressed(out, time=np.array([r[0] for r in rows]), lat=LAT, lon=LON, d20=np.array([r[1] for r in rows]))
    shutil.rmtree(tmp); os.remove(z)
    print(time.strftime("%H:%M:%S"), y, len(rows), "months", flush=True)


if __name__ == "__main__":
    y0, y1 = int(sys.argv[1]), int(sys.argv[2]); cl = cdsapi.Client(quiet=True, progress=False)
    for y in range(y0, y1 + 1):
        try:
            year(y, cl)
        except Exception as e:  # noqa: BLE001
            print(time.strftime("%H:%M:%S"), y, "FAILED:", str(e)[:200], flush=True)
    print("done", flush=True)
