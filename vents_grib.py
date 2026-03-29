import matplotlib.pyplot as plt
from math import *
import numpy as np
import xarray as xr
import pandas as pd

fichier = "grib_vent/20260327_143145_AROME_P025_06.grb2"
# fichier = "grib_vent/gfs_0_25.2025-11-29T16-50-54Z.grb2"
# fichier = "grib_vent/20260304_105640_GFS_P25_06.grb2"
# fichier = "grib_vent/gfs_025_20260306_093631.grb2"

ds_u10 = xr.open_dataset(fichier, engine="cfgrib", backend_kwargs={"indexpath": "", "filter_by_keys": {"shortName": "10u"}})
ds_v10 = xr.open_dataset(fichier , engine="cfgrib", backend_kwargs={"indexpath": "", "filter_by_keys": {"shortName": "10v"}})

df_u10 = ds_u10.to_dataframe().reset_index()
df_v10 = ds_v10.to_dataframe().reset_index()

df = pd.merge(df_u10, df_v10, on=["valid_time", "latitude", "longitude", "step", "time"], how="outer")

columns_to_drop = ['step', 'time', 'surface', 'heightAboveGround_x', 'heightAboveGround_y']
df = df.drop(columns=columns_to_drop, errors='ignore')

df['force'] = np.sqrt(df['u10']**2 + df['v10']**2) * 1.94384
df['direction'] = np.degrees(np.arctan2(df['u10'], df['v10'])) + 180
df['direction'] = df['direction'] % 360

columns_to_drop = ['u10', 'v10']
df = df.drop(columns=columns_to_drop, errors='ignore')

df['temps'] = (df['valid_time'].apply(lambda x: x.timestamp()) - df['valid_time'][0].timestamp()) / 3600.0

# columns_to_drop = ['valid_time']
# df = df.drop(columns=columns_to_drop, errors='ignore')

# df['longitude'] = 360 - df['longitude']

df['x_data'] = df['longitude'] * 60 * 0.7
df['y_data'] = df['latitude'] * 60

# R = 3443.9184665
# df['x_data'] = R * np.deg2rad(df['longitude'])
# df['y_data'] = R * np.log(np.tan(np.pi/4 + np.deg2rad(df['latitude'])/2))


# columns_to_drop = ['longitude', 'latitude']
# df = df.drop(columns=columns_to_drop, errors='ignore')







from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator

points = np.array(df[['x_data', 'y_data', 'temps']])
values = np.array(df[['direction', 'force']])

interp_func = NearestNDInterpolator(points, values)

# def gfs(p, t):
#     x, y = p[0], p[1]
#     if x > 15125 or x < 14700 :
#         return np.array([0, 0])
#     if y > 2930 or y < 2520 :
#         return np.array([0, 0])
#     return interp_func(x, y, t)

points_lat = np.array(df[['longitude', 'latitude', 'temps']])
interp_func_lat = NearestNDInterpolator(points_lat, values)

# def gfs_lat(p, t):
#     x, y = p[0], p[1]
#     return interp_func_lat(x, y, t)

import numpy as np
import cartopy.io.shapereader as shpreader
from shapely.ops import unary_union
from shapely import contains_xy, prepare

# Charge les polygones de terre Natural Earth
land_shp = shpreader.natural_earth(
    resolution='10m',
    category='physical',
    name='land'
)

reader = shpreader.Reader(land_shp)
land_geom = unary_union(list(reader.geometries()))
prepare(land_geom)   # utile pour des tests répétés sur la même géométrie

def valide_nm(x, y):
    lat = y / 60.0
    lon = x / (60.0 * 0.7) - 360.0
    return not contains_xy(land_geom, lon, lat)

def valide_deg(lon, lat):
    return not contains_xy(land_geom, lon - 360.0, lat)

def gfs(p, t, ref):
    x, y = p[0], p[1]

    if ref == "nm":
        return interp_func(x, y, t) if valide_nm(x, y) else np.array([0.0, 0.0])

    if ref == "deg":
        return interp_func_lat(x, y, t) if valide_deg(x, y) else np.array([0.0, 0.0])

    raise ValueError("ref doit être 'nm' ou 'deg'")

# def gfs(p, t, ref):
#     x, y = p[0], p[1]
#     if ref == 'nm':
#         return interp_func(x, y, t)
#     if ref == 'deg':
#         return interp_func_lat(x, y, t)
    

# interp_func_1 = LinearNDInterpolator(points, values)
# interp_func_lat_1 = LinearNDInterpolator(points_lat, values)

# def gfs_lin(p, t, ref):
#     x, y = p[0], p[1]
#     if ref == 'nm':
#         return interp_func_1(x, y, t)
#     if ref == 'deg':
#         return interp_func_lat_1(x, y, t)