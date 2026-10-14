"""ERA5 equatorial zonal wind at 30, 50 and 70 hPa (from getERA5.py, levels reduced; 2026-10-13).
Retrieves monthly means 1979-present over 5N-5S, then writes one zonal-mean series per level:
    era5_qbo_u30.dat, era5_qbo_u50.dat, era5_qbo_u70.dat      (decimal year, m/s)
Needs Copernicus credentials in ~/.cdsapirc:
    url: https://cds.climate.copernicus.eu/api
    key: <your personal access token>
(the current CDS needs cdsapi >= 0.7.2: pip install -U cdsapi). The licence for the dataset must have been accepted
once on the CDS web site."""
import cdsapi, numpy as np, netCDF4 as nc
OUT = "era5_qbo_equatorial_u_30_50_70.nc"; LEVELS = ["30", "50", "70"]
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'product_type': ['monthly_averaged_reanalysis'],
        'data_format': 'netcdf', 'format': 'netcdf',
        'variable': 'u_component_of_wind',
        'pressure_level': LEVELS,
        'year': [str(yr) for yr in range(1979, 2027)],
        'month': [f'{m:02d}' for m in range(1, 13)],
        'time': '00:00',
        'area': [5, -180, -5, 180],  # Equator: 5N to 5S zonal band
    },
    OUT)
ds = nc.Dataset(OUT); u = ds.variables["u"]; tv = ds.variables["valid_time" if "valid_time" in ds.variables else "time"]; lev = ds.variables["pressure_level" if "pressure_level" in ds.variables else "level"][:]
dates = nc.num2date(tv[:], tv.units); dec = np.array([d.year + (d.month - 1) / 12.0 for d in dates]); lat = ds.variables["latitude"][:]; w = np.cos(np.deg2rad(lat)); U = np.asarray(u[:])          # (time, level, lat, lon)
zm = (U.mean(axis=3) * w[None, None, :]).sum(axis=2) / w.sum()
for j, L in enumerate(lev):
    np.savetxt(f"era5_qbo_u{int(L)}.dat", np.column_stack([dec, zm[:, j]]), fmt=["%.6f", "%.4f"], delimiter="\t"); print(f"era5_qbo_u{int(L)}.dat: {len(dec)} months {dec[0]:.2f}-{dec[-1]:.2f}, sd {zm[:, j].std():.1f} m/s")
