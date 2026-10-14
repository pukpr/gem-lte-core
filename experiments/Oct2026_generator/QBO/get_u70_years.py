"""ERA5 u at 70 hPa, 6-hourly, 2.5 degree grid, 5S-5N, for a few years only (comparison with qbo70era0.dat)."""
import sys, os, cdsapi
c = cdsapi.Client(quiet=True)
for yr in sys.argv[1:]:
    fn = f"era5_u70_daily_{yr}.nc"
    if not os.path.exists(fn):
        c.retrieve('reanalysis-era5-pressure-levels', {'product_type': ['reanalysis'], 'variable': ['u_component_of_wind'], 'pressure_level': ['70'], 'year': [yr],
                   'month': [f'{m:02d}' for m in range(1, 13)], 'day': [f'{d:02d}' for d in range(1, 32)], 'time': ['00:00', '06:00', '12:00', '18:00'],
                   'area': [5, -180, -5, 180], 'grid': ['2.5', '2.5'], 'data_format': 'netcdf', 'download_format': 'unarchived'}, fn)
    print(yr, "done", flush=True)
