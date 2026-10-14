"""ERA5 u at one pressure level and one grid point, 6-hourly, 1979-present, in multi-year requests (few requests, tiny
files): writes era_qbo_u<L>_point<LON>E_daily.dat and _5day.dat (73-per-year grid, >= 3 days per bin).
usage: .venv/bin/python get_era5_point.py LEVEL [LAT 0] [LON 0] [YEARS_PER_REQUEST 8] [TIMES 00]
The request cost is the number of fields read (area does not reduce it), so fewer analysis times per day is what
makes it cheap: TIMES 00 = one analysis a day (00 UTC); TIMES 00,06,12,18 = four (daily means)."""
import sys, os, datetime as dt, cdsapi, numpy as np, netCDF4 as nc
L = sys.argv[1]; LAT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0; LON = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0; NY = int(sys.argv[4]) if len(sys.argv) > 4 else 8; TIMES = [f'{h}:00' for h in (sys.argv[5] if len(sys.argv) > 5 else '00').split(',')]
c = cdsapi.Client(quiet=True); last = dt.date.today().year; day = {}
for y0 in range(1979, last + 1, NY):
    yrs = [str(y) for y in range(y0, min(y0 + NY, last + 1))]; fn = f"era5_u{L}_point_{LAT:g}_{LON:g}_{yrs[0]}_{yrs[-1]}_{len(TIMES)}x.nc"
    if not os.path.exists(fn):
        c.retrieve('reanalysis-era5-pressure-levels', {'product_type': ['reanalysis'], 'variable': ['u_component_of_wind'], 'pressure_level': [L], 'year': yrs,
                   'month': [f'{m:02d}' for m in range(1, 13)], 'day': [f'{d:02d}' for d in range(1, 32)], 'time': TIMES,
                   'area': [LAT + 2.5, LON - 2.5, LAT - 2.5, LON + 2.5], 'grid': ['2.5', '2.5'], 'data_format': 'netcdf', 'download_format': 'unarchived'}, fn)
    ds = nc.Dataset(fn); la = np.asarray(ds.variables["latitude"][:]); lo = np.asarray(ds.variables["longitude"][:]); j = int(np.argmin(np.abs(la - LAT))); i = int(np.argmin(np.abs(((lo - LON + 180) % 360) - 180)))
    u = np.asarray(ds.variables["u"][:, 0, j, i]).reshape(-1); tv = ds.variables["valid_time"]   # a box of 3 x 3 points is requested: a zero-size area is treated by the server as the whole globe
    for d, v in zip(nc.num2date(tv[:], tv.units), u): day.setdefault((d.year, d.month, d.day), []).append(float(v))
    print(yrs[0], "-", yrs[-1], "done", flush=True)
rows = [(y + (dt.date(y, m, d) - dt.date(y, 1, 1)).days / ((dt.date(y + 1, 1, 1) - dt.date(y, 1, 1)).days), np.mean(v)) for (y, m, d), v in sorted(day.items())]
A = np.array(rows); tag = f"era_qbo_u{L}_point{LON:g}E"; np.savetxt(f"{tag}_daily.dat", A, fmt=["%.6f", "%.5f"], delimiter="\t"); y0 = np.floor(A[0, 0]); idx = np.floor((A[:, 0] - y0) * 73 + 1e-9).astype(int)
F = np.array([(y0 + k / 73.0, A[idx == k, 1].mean()) for k in range(idx.min(), idx.max() + 1) if (idx == k).sum() >= 3]); np.savetxt(f"{tag}_5day.dat", F, fmt=["%.6f", "%.4f"], delimiter="\t"); print(f"{len(A)} days, {len(F)} five-day means, {A[0, 0]:.2f}-{A[-1, 0]:.2f}", flush=True)
