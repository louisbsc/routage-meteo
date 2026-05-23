import matplotlib.pyplot as plt
from math import *

# VENTS THÉORIQUES

def vent_constant(p, t):
	return 0, 5

def vent_uniforme(direction, force):
	dir_, force_ = float(direction), float(force)
	def vent(p, t):
		p = np.asarray(p)
		batch = p.ndim == 2
		xs = p[:, 0] if batch else p[0:1]
		ys = p[:, 1] if batch else p[1:2]
		dirs   = np.full(len(xs), dir_)
		forces = np.full(len(xs), force_)
		land = contains_xy(land_geom, xs / (60.0 * 0.7), ys / 60.0)
		dirs[land]   = 0.0
		forces[land] = 0.0
		result = np.column_stack([dirs, forces])
		return result if batch else result[0]
	return vent

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
from scipy.interpolate import RegularGridInterpolator

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

    df['temps'] = (df['valid_time'].apply(lambda x: x.timestamp()) - df['valid_time'][0].timestamp()) / 3600.0

    df['longitude'] = df['longitude'].apply(lambda x: x - 360 if x > 180 else x)

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
def _build_uv_interpolators(df, x_col, y_col):
    xs = np.sort(df[x_col].unique())
    ys = np.sort(df[y_col].unique())
    ts = np.sort(df['temps'].unique())

    xi = np.searchsorted(xs, df[x_col].values)
    yi = np.searchsorted(ys, df[y_col].values)
    ti = np.searchsorted(ts, df['temps'].values)

    u_grid = np.zeros((len(xs), len(ys), len(ts)))
    v_grid = np.zeros((len(xs), len(ys), len(ts)))
    u_grid[xi, yi, ti] = df['u10'].values
    v_grid[xi, yi, ti] = df['v10'].values

    kw = dict(method='linear', bounds_error=False, fill_value=0.0)
    return (
        RegularGridInterpolator((xs, ys, ts), u_grid, **kw),
        RegularGridInterpolator((xs, ys, ts), v_grid, **kw),
    )

def vent_grib_nm(path):
    df = table(path)
    t_max_grib = float(df['temps'].max())
    interp_u, interp_v = _build_uv_interpolators(df, 'x_data', 'y_data')

    def vent(p, t):
        if t > t_max_grib:
            raise ValueError("grib trop court en date")
        p = np.asarray(p)
        batch = p.ndim == 2
        xs = p[:, 0] if batch else p[0:1]
        ys = p[:, 1] if batch else p[1:2]

        pts = np.column_stack([xs, ys, np.full(len(xs), t)])
        u = interp_u(pts)
        v = interp_v(pts)
        land = contains_xy(land_geom, xs / (60.0 * 0.7), ys / 60.0)
        u[land] = 0.0
        v[land] = 0.0

        force = np.sqrt(u**2 + v**2) * 1.94384
        direction = (np.degrees(np.arctan2(u, v)) + 180) % 360
        result = np.column_stack([direction, force])
        return result if batch else result[0]

    return vent

def vent_grib_deg(path):
    df = table(path)
    interp_u, interp_v = _build_uv_interpolators(df, 'longitude', 'latitude')

    def vent(p, t):
        p = np.asarray(p)
        batch = p.ndim == 2
        xs = p[:, 0] if batch else p[0:1]
        ys = p[:, 1] if batch else p[1:2]

        pts = np.column_stack([xs, ys, np.full(len(xs), t)])
        u = interp_u(pts)
        v = interp_v(pts)
        land = contains_xy(land_geom, xs, ys)
        u[land] = 0.0
        v[land] = 0.0

        force = np.sqrt(u**2 + v**2) * 1.94384
        direction = (np.degrees(np.arctan2(u, v)) + 180) % 360
        result = np.column_stack([direction, force])
        return result if batch else result[0]

    return vent



def vent_filtre_route(V, route, seuil):
    """
    Enveloppe une fonction vent V(p, t) en mettant direction et force à zéro
    pour tout point dont la distance à la polyligne route dépasse seuil.

    V      : fonction vent — même signature que vent_grib_nm / vent_grib_deg
    route  : (K, ≥2) waypoints dans le même système de coordonnées que p
    seuil  : float  (même unité que les coordonnées de V)
    """
    from core.utils import point_proche_route
    route_ = np.asarray(route, dtype=float)

    def vent(p, t):
        p_arr  = np.asarray(p)
        batch  = p_arr.ndim == 2
        pts    = p_arr if batch else p_arr[np.newaxis, :]    # (M, ≥2)
        result = np.array(V(p, t), dtype=float)
        if not batch:
            result = result[np.newaxis, :]                   # (1, 2)
        mask           = point_proche_route(route_, pts, seuil)  # (M,) bool
        result[~mask]  = 0.0
        return result if batch else result[0]

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


