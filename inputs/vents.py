import matplotlib.pyplot as plt
from math import *

# VENTS THÉORIQUES

def vent_constant(p, t):
	return 0, 5

def vent_circulaire(p, t):
    x, y = p[0], p[1]

    r = hypot(x, y)
    if r == 0:
        return 0.0, 0.0  # Direction indéfinie au centre

    # Force du vent : 10 kts proche du centre, asymptote à 20 kts
    # r exprimé en NM (supposé)
    force = 20 * (1 - exp(-r / 40))

    # Vent tangent enroulement horaire : (y, -x)
    vx = force * (y / r)
    vy = -force * (x / r)

    # Angle mathématique du vecteur (0 rad = +x, CCW)
    angle_rad = atan2(vy, vx)

    # Conversion en angle "rose des vents" (0° = Nord venant du haut, sens horaire)
    direction_deg = (90 - degrees(angle_rad)) % 360

    return direction_deg, force

# VENTS GRIB

import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import NearestNDInterpolator

# transforme un fichier grib en table pandas propre
def table(path):
    ds_u10 = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": "", "filter_by_keys": {"shortName": "10u"}})
    ds_v10 = xr.open_dataset(path , engine="cfgrib", backend_kwargs={"indexpath": "", "filter_by_keys": {"shortName": "10v"}})
    df_u10 = ds_u10.to_dataframe().reset_index()
    df_v10 = ds_v10.to_dataframe().reset_index()
    df = pd.merge(df_u10, df_v10, on=["valid_time", "latitude", "longitude", "step", "time"], how="outer")

    columns_to_drop = ['step', 'time', 'surface', 'heightAboveGround_x', 'heightAboveGround_y']
    df = df.drop(columns=columns_to_drop, errors='ignore')

    # m/s en noeuds
    df['force'] = np.sqrt(df['u10']**2 + df['v10']**2) * 1.94384
    df['direction'] = np.degrees(np.arctan2(df['u10'], df['v10'])) + 180
    df['direction'] = df['direction'] % 360

    columns_to_drop = ['u10', 'v10']
    df = df.drop(columns=columns_to_drop, errors='ignore')

    df['temps'] = (df['valid_time'].apply(lambda x: x.timestamp()) - df['valid_time'][0].timestamp()) / 3600.0
    
    # transformation plane avec latitude moyenne
    df['x_data'] = df['longitude'] * 60 * 0.7
    df['y_data'] = df['latitude'] * 60

    return df

# renvoie un point s'il est bien en mer et non à terre

import cartopy.io.shapereader as shpreader
from shapely.ops import unary_union
from shapely import contains_xy, prepare

land_shp = shpreader.natural_earth(
    resolution='10m',
    category='physical',
    name='land'
)

reader = shpreader.Reader(land_shp)
land_geom = unary_union(list(reader.geometries()))
prepare(land_geom)   # utile pour des tests répétés sur la même géométrie


# renvoie la fonction vent correspondant à un grib dont le chemin est path
def vent_grib_nm(path):
    df = table(path)
    points = np.array(df[['x_data', 'y_data', 'temps']])
    values = np.array(df[['direction', 'force']])

    interp_func = NearestNDInterpolator(points, values)

    def vent(p, t):
        x, y = p[0], p[1]
        return interp_func(x, y, t) if (not contains_xy(land_geom, x / (60.0 * 0.7) - 360.0, y / 60.0)) else np.array([0.0, 0.0])

    return vent

def vent_grib_deg(path):
    df = table(path)
    points = np.array(df[['longitude', 'latitude', 'temps']])
    values = np.array(df[['direction', 'force']])

    interp_func = NearestNDInterpolator(points, values)

    def vent(p, t):
        x, y = p[0], p[1]
        return interp_func(x, y, t) if (not contains_xy(land_geom, x - 360.0, y)) else np.array([0.0, 0.0])

    return vent



# AFFICHAGE

def aff_vent_static(V, t, intx, inty, n):
	x = (n + 1) * [intx[0] + k * (intx[1] - intx[0]) / n for k in range(n + 1)]
	y = []
	for k in range(n + 1):
		y = y + (n + 1) * [inty[0] + k * (inty[1] - inty[0]) / n]
	u = [- V([x[k], y[k]], t)[1] * sin(radians(V([x[k], y[k]], t)[0])) for k in range((n + 1)**2)]
	v = [- V([x[k], y[k]], t)[1] * cos(radians(V([x[k], y[k]], t)[0])) for k in range((n + 1)**2)]
	c = []
	for k in range((n + 1)**2):
		coul = V([x[k], y[k]], t)[1] / 30
		c.append((1 - coul, 1 - coul, 1 - coul))
	plt.quiver(x, y, u, v, color = c)
	plt.axis('equal')


