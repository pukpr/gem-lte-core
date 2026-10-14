"""ERA5 DAILY-mean equatorial zonal wind at one pressure level (default 50 hPa), 1979-present, for the 5-day leg.
Writes  era_qbo_u<L>_daily.dat  (calendar decimal year, m/s; zonal mean over 5N-5S) and
        era_qbo_u<L>_5day.dat   (means on the 73-per-year grid used for the SOI, >= 3 days per bin).
Needs Copernicus credentials in ~/.cdsapirc and cdsapi >= 0.7.2 (see getERA5_30_50_70.py). Daily means are formed here from the 00, 06, 12 and 18 UTC analyses on a 2.5 degree grid.
Run with the local environment: .venv/bin/python getERA5_daily_50.py
usage: python3 getERA5_daily_50.py [level, default 50] [first year, default 1979]"""
import sys, os, datetime as dt, cdsapi, numpy as np, netCDF4 as nc
L = sys.argv[1] if len(sys.argv) > 1 else "50"; Y0 = int(sys.argv[2]) if len(sys.argv) > 2 else 1979; c = cdsapi.Client(); rows = []
for yr in range(Y0, dt.date.today().year + 1):
    fn = f"era5_u{L}_daily_{yr}.nc"
    if not os.path.exists(fn):
        # 6-hourly analyses on a 2.5 degree grid (5 latitudes x 144 longitudes): about 4 MB a year
        c.retrieve('reanalysis-era5-pressure-levels',
                   {'product_type': ['reanalysis'], 'variable': ['u_component_of_wind'], 'pressure_level': [L], 'year': [str(yr)],
                    'month': [f'{m:02d}' for m in range(1, 13)], 'day': [f'{d:02d}' for d in range(1, 32)],
                    'time': ['00:00', '06:00', '12:00', '18:00'], 'area': [5, -180, -5, 180], 'grid': ['2.5', '2.5'],
                    'data_format': 'netcdf', 'download_format': 'unarchived'}, fn)
    ds = nc.Dataset(fn); u = np.asarray(ds.variables["u"][:]).squeeze(); tv = ds.variables["valid_time" if "valid_time" in ds.variables else "time"]; lat = np.asarray(ds.variables["latitude"][:]); w = np.cos(np.deg2rad(lat))
    zm = (u.mean(axis=-1) * w).sum(axis=-1) / w.sum(); day = {}
    for d, v in zip(nc.num2date(tv[:], tv.units), zm): day.setdefault((d.year, d.month, d.day), []).append(float(v))
    for (y_, m_, d_), vv in sorted(day.items()):
        y0 = dt.date(y_, 1, 1); rows.append((y_ + (dt.date(y_, m_, d_) - y0).days / ((dt.date(y_ + 1, 1, 1) - y0).days), float(np.mean(vv))))
    print(yr, len(day), "days", flush=True)
A = np.array(rows); np.savetxt(f"era_qbo_u{L}_daily.dat", A, fmt=["%.6f", "%.5f"], delimiter="\t"); y0 = np.floor(A[0, 0]); idx = np.floor((A[:, 0] - y0) * 73 + 1e-9).astype(int)
F = np.array([(y0 + i / 73.0, A[idx == i, 1].mean()) for i in range(idx.min(), idx.max() + 1) if (idx == i).sum() >= 3]); np.savetxt(f"era_qbo_u{L}_5day.dat", F, fmt=["%.6f", "%.4f"], delimiter="\t"); print(f"{len(A)} days, {len(F)} five-day means, {A[0, 0]:.2f}-{A[-1, 0]:.2f}")
