"""ERA5 HOURLY surface pressure over the tropical band 20S-20N on a 5 degree grid, one file per year (2026-10-25).
For the eastward zonal-wavenumber-1 spectrum (Sakazaki & Hamilton 2020, JAS, "An Array of Ringing Global Free Modes
Discovered in Tropical Surface Pressure Data"; discussed at geoenergymath.com/the-just-so-story-narrative/).
A 5 degree grid is ample for wavenumber 1: 9 latitudes x 72 longitudes x 8760 hours, about 25 MB a year.
Run with the QBO environment (cdsapi >= 0.7.2):  ../../experiments/Oct2026_generator/QBO/.venv/bin/python fetch_sp.py [y0] [y1]
Credentials in ~/.cdsapirc."""
import os, sys
import cdsapi

Y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1979
Y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2024
HERE = os.path.dirname(os.path.abspath(__file__))
c = cdsapi.Client()
for yr in range(Y0, Y1 + 1):
    for half, months in (("a", range(1, 7)), ("b", range(7, 13))):   # half years keep each request under the field limit
        fn = os.path.join(HERE, "data", f"era5_sp_hourly_{yr}{half}.nc")
        if os.path.exists(fn):
            continue
        c.retrieve("reanalysis-era5-single-levels",
                   {"product_type": ["reanalysis"], "variable": ["surface_pressure"], "year": [str(yr)],
                    "month": [f"{m:02d}" for m in months], "day": [f"{d:02d}" for d in range(1, 32)],
                    "time": [f"{h:02d}:00" for h in range(24)], "area": [20, 0, -20, 355], "grid": ["5", "5"],
                    "data_format": "netcdf", "download_format": "unarchived"}, fn + ".part")
        os.replace(fn + ".part", fn)
        print(yr, half, "done", flush=True)
