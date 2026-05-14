import numpy as np
import eccodes
import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator

def C0(p, t):
    return 0, 0

def C1(p, t):
    x, y = p[:2]
    if x < 0:
        return 180, 1
    else:
        return 180, 0
    
def C2(p, t):
    if t < 10:
        return 180, 1
    else:
        return 180, 0


# transforme un fichier grib en table pandas propre
def table(path):
    records = []

    with open(path, "rb") as f:
        while True:
            msg = eccodes.codes_grib_new_from_file(f)
            if msg is None:
                break

            param_num = eccodes.codes_get(msg, "parameterNumber")
            # 2 = U (courant Est), 3 = V (courant Nord)
            if param_num not in (2, 3):
                eccodes.codes_release(msg)
                continue

            step      = eccodes.codes_get(msg, "stepRange")
            date      = eccodes.codes_get(msg, "dataDate")   # ex: 20260514
            time      = eccodes.codes_get(msg, "dataTime")   # ex: 0 (heure UTC)
            ni        = eccodes.codes_get(msg, "Ni")          # nb points longitude
            nj        = eccodes.codes_get(msg, "Nj")          # nb points latitude
            lat_first = eccodes.codes_get(msg, "latitudeOfFirstGridPointInDegrees")
            lon_first = eccodes.codes_get(msg, "longitudeOfFirstGridPointInDegrees")
            d_lat     = eccodes.codes_get(msg, "jDirectionIncrementInDegrees")
            d_lon     = eccodes.codes_get(msg, "iDirectionIncrementInDegrees")

            values = eccodes.codes_get_values(msg)  # array 1D, taille Ni×Nj
            # Remplacer les valeurs égales à 9999 par 0
            values = np.where(values == 9999, 0, values)
            eccodes.codes_release(msg)

            # Reconstruire les grilles lat/lon
            lats = lat_first + np.arange(nj) * d_lat
            lons = lon_first + np.arange(ni) * d_lon
            # Normaliser les longitudes > 180° en négatif
            lons = np.where(lons > 180, lons - 360, lons)

            # Meshgrid → tableau 2D (nj, ni)
            LON, LAT = np.meshgrid(lons, lats)

            component = "u_current" if param_num == 2 else "v_current"

            df_step = pd.DataFrame({
                "date":      date,
                "step_h":    int(step),
                "lat":       LAT.ravel(),
                "lon":       LON.ravel(),
                component:   values,
            })
            records.append(df_step)

    # Assembler tous les steps
    df_all = pd.concat(records, ignore_index=True)

    # Pivoter pour avoir u et v sur la même ligne
    df_u = df_all[df_all["u_current"].notna()][["date","step_h","lat","lon","u_current"]]
    df_v = df_all[df_all["v_current"].notna()][["date","step_h","lat","lon","v_current"]]

    df = pd.merge(df_u, df_v, on=["date","step_h","lat","lon"])

    # Remplacer à nouveau au cas où il resterait des 9999 après fusion (sécurité)
    df["u_current"] = np.where(df["u_current"] == 9999, 0, df["u_current"])
    df["v_current"] = np.where(df["v_current"] == 9999, 0, df["v_current"])

    # Calculer intensité (m/s) et direction (° météo, 0°=Nord, sens horaire)
    df["speed_ms"]   = np.sqrt(df["u_current"]**2 + df["v_current"]**2)
    df["direction"]  = (np.degrees(np.arctan2(df["u_current"], df["v_current"])) + 360) % 360

    # transformation plane avec latitude moyenne
    df['x_data'] = df['lon'] * 60 * 0.7
    df['y_data'] = df['lat'] * 60

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
    ts = np.sort(df['step_h'].unique())

    xi = np.searchsorted(xs, df[x_col].values)
    yi = np.searchsorted(ys, df[y_col].values)
    ti = np.searchsorted(ts, df['step_h'].values)

    u_grid = np.zeros((len(xs), len(ys), len(ts)))
    v_grid = np.zeros((len(xs), len(ys), len(ts)))
    u_grid[xi, yi, ti] = df['u_current'].values
    v_grid[xi, yi, ti] = df['v_current'].values

    kw = dict(method='linear', bounds_error=False, fill_value=0.0)
    return (
        RegularGridInterpolator((xs, ys, ts), u_grid, **kw),
        RegularGridInterpolator((xs, ys, ts), v_grid, **kw),
    )

def courant_grib_nm(path):
    df = table(path)
    t_max_grib = float(df['step_h'].max())
    interp_u, interp_v = _build_uv_interpolators(df, 'x_data', 'y_data')

    def courant(p, t):
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

    return courant


def courant_grib_deg(path):
    df = table(path)
    interp_u, interp_v = _build_uv_interpolators(df, 'lon', 'lat')

    def courant(p, t):
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

    return courant


def courant_uniforme(direction, force):
    """Courant uniforme. Convention identique au vent : direction = provenance (°), force en nœuds."""
    dir_  = float(direction)
    force_ = float(force)

    def courant(p, t):
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

    return courant