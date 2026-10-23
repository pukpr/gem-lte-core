"""ERA5 monthly-mean zonal wind, 100-200 hPa, 5S-5N, 2.5 deg grid, 1979-2025 (2026-10-22). -> era5_u_eq_monthly_low.nc
Run with the QBO environment: ../QBO/.venv/bin/python fetch_sao.py"""
import cdsapi
cdsapi.Client(quiet=True, progress=False).retrieve("reanalysis-era5-pressure-levels-monthly-means", {
    "product_type": ["monthly_averaged_reanalysis"], "variable": ["u_component_of_wind"],
    "pressure_level": ["100", "125", "150", "200"],
    "year": [str(y) for y in range(1979, 2026)], "month": [f"{m:02d}" for m in range(1, 13)], "time": ["00:00"],
    "area": [5, -180, -5, 180], "grid": [2.5, 2.5], "data_format": "netcdf", "download_format": "unarchived"},
    "era5_u_eq_monthly_low.nc")
