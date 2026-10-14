import cdsapi

c = cdsapi.Client()

# Downloads monthly-averaged equatorial zonal wind profiles (1979–present)
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'variable': 'u_component_of_wind',
        'pressure_level': [
            '10', '15', '20', '30', '40', '50', '70',
        ],
        'year': [str(yr) for yr in range(1979, 2027)],
        'month': [f'{m:02d}' for m in range(1, 13)],
        'time': '00:00',
        'area': [5, -180, -5, 180],  # Equator: 5N to 5S zonal band
    },
    'era5_qbo_equatorial_u_1979_present.nc'
)
