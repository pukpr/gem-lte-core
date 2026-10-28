"""ERA5 HOURLY stratospheric winds over the tropical band 20S-20N, 5 degree grid (2026-10-26).
Companion to fetch_sp.py: same grid and band, so the wavenumber / symmetric-antisymmetric tidal analysis of surface
pressure can be repeated on the QBO levels. One request per variable, level and half year: CDS refuses 2 or 4 levels x
a half year ("cost limits exceeded", >= ~8 700 fields).  Writes data/era5_<var><levels>_hourly_<year><a|b>.nc,
e.g. era5_u70_hourly_2016a.nc.
Levels 70, 50, 30, 10 hPa.  Run with the QBO environment (cdsapi >= 0.7.2):
    ../../experiments/Oct2026_generator/QBO/.venv/bin/python fetch_uv.py [y0] [y1] [u|v|uv] [levels, e.g. 30,50]"""
import os, sys
import cdsapi

Y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1999
Y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2016
VARS = {"u": "u_component_of_wind", "v": "v_component_of_wind"}
WHICH = sys.argv[3] if len(sys.argv) > 3 else "uv"
LEVELS = sys.argv[4].split(",") if len(sys.argv) > 4 else ["70", "50", "30", "10"]
HERE = os.path.dirname(os.path.abspath(__file__))
c = cdsapi.Client()
for yr in range(Y0, Y1 + 1):
    for half, months in (("a", range(1, 7)), ("b", range(7, 13))):
        for v in WHICH:
          for lev in ([L] for L in LEVELS):
            fn = os.path.join(HERE, "data", f"era5_{v}{'-'.join(lev)}_hourly_{yr}{half}.nc")
            if os.path.exists(fn):
                continue
            c.retrieve("reanalysis-era5-pressure-levels",
                       {"product_type": ["reanalysis"], "variable": [VARS[v]], "pressure_level": lev,
                        "year": [str(yr)], "month": [f"{m:02d}" for m in months], "day": [f"{d:02d}" for d in range(1, 32)],
                        "time": [f"{h:02d}:00" for h in range(24)], "area": [20, 0, -20, 355], "grid": ["5", "5"],
                        "data_format": "netcdf", "download_format": "unarchived"}, fn + ".part")
            os.replace(fn + ".part", fn)
            print(yr, half, v, lev, "done", flush=True)
