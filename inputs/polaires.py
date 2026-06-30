import math as _math

import numba
import numpy as np
import matplotlib.pyplot as plt
from math import *

# POLAIRES THÉORIQUES

def P0(ang, f):
	ang_arr = np.asarray(ang)
	f_arr = np.asarray(f)
	shape = np.broadcast(ang_arr, f_arr).shape
	return np.full(shape, 5)

def P1(ang, f):
	return  15 * (1 - exp(-f/10))

def P2aux_vec(ang):
    ang = np.mod(ang, 360)

    ang_rad = np.radians(ang)
    cos_ang = np.cos(ang_rad)
    sin_ang = np.sin(ang_rad)

    out = np.zeros_like(ang, dtype=float)

    # Cas 1 :  ang <= 45   ou   135 ≤ ang ≤ 225   ou   ang ≥ 315
    mask1 = (ang <= 45) | ((ang >= 135) & (ang <= 225)) | (ang >= 315)
    out[mask1] = 3 / np.abs(cos_ang[mask1])

    # Cas 2 : 45 < ang < 90
    mask2 = (ang > 45) & (ang < 90)
    out[mask2] = 6 * sin_ang[mask2]

    # Cas 3 : 225 < ang < 270
    mask3 = (ang > 225) & (ang < 270)
    out[mask3] = -6 * sin_ang[mask3]

    # Cas 4 : 270 ≤ ang < 315
    mask4 = (ang >= 270) & (ang < 315)
    out[mask4] = -6 * sin_ang[mask4]

    # Cas 5 : 90 ≤ ang < 135
    mask5 = (ang >= 90) & (ang < 135)
    out[mask5] = 6 * sin_ang[mask5]

    return out


def P2(ang, f):
    # f peut être un scalaire ou un array broadcastable
    return P2aux_vec(ang) * 2 * (1 - np.exp(-f / 10))

def polaire_uniforme(c):
    def P(ang, f):
        ang_arr = np.asarray(ang)
        f_arr = np.asarray(f)
        shape = np.broadcast(ang_arr, f_arr).shape
        out = np.full(shape, float(c))
        out[np.broadcast_to(f_arr, shape) == 0] = 0.0
        return out
    P.v_max = float(2 * c)
    return P

# POLAIRES RÉELLES

import pandas as pd


@numba.njit(cache=True, parallel=True)
def _bilinear_interp(xa, xs, v, qa, qs):
	"""Interpolation bilinéaire compilée — remplace RegularGridInterpolator.

	xa: (Na,) angles grille triés  xs: (Ns,) vitesses grille triées
	v:  (Na, Ns) vitesses bateau   qa/qs: (N,) requêtes
	Retourne (N,) — 0 hors bornes (fill_value=0).
	"""
	Na  = xa.shape[0]
	Ns  = xs.shape[0]
	N   = qa.shape[0]
	out = np.zeros(N)

	for i in numba.prange(N):
		a = qa[i]
		s = qs[i]
		if a < xa[0] or a > xa[Na - 1] or s < xs[0] or s > xs[Ns - 1]:
			continue

		# Recherche binaire sur l'angle
		lo_a, hi_a = 0, Na - 1
		while lo_a < hi_a - 1:
			mid = (lo_a + hi_a) >> 1
			if xa[mid] <= a:
				lo_a = mid
			else:
				hi_a = mid
		dxa = xa[hi_a] - xa[lo_a]
		ta  = 0.0 if dxa < 1e-15 else (a - xa[lo_a]) / dxa

		# Recherche binaire sur la vitesse
		lo_s, hi_s = 0, Ns - 1
		while lo_s < hi_s - 1:
			mid = (lo_s + hi_s) >> 1
			if xs[mid] <= s:
				lo_s = mid
			else:
				hi_s = mid
		dxs = xs[hi_s] - xs[lo_s]
		ts  = 0.0 if dxs < 1e-15 else (s - xs[lo_s]) / dxs

		v00 = v[lo_a, lo_s]
		v10 = v[hi_a, lo_s]
		v01 = v[lo_a, hi_s]
		v11 = v[hi_a, hi_s]
		out[i] = (1.0 - ta) * ((1.0 - ts) * v00 + ts * v01) + \
		          ta          * ((1.0 - ts) * v10 + ts * v11)

	return out


def polaire(path):
	df          = pd.read_csv(path)
	wind_angles = np.ascontiguousarray(df.iloc[:, 0].astype(float).values)
	wind_speeds = np.ascontiguousarray(df.columns[1:].astype(float).values)
	boat_speeds = np.ascontiguousarray(df.iloc[:, 1:].astype(float).values)

	# Warmup JIT à la définition — pas de latence au premier appel routage
	_bilinear_interp(wind_angles, wind_speeds, boat_speeds,
	                 np.array([90.0]), np.array([10.0]))

	def polaire_func(ang, f):
		a     = np.abs(np.asarray(ang, dtype=np.float64))
		f_col = (np.full(len(a), float(f), dtype=np.float64)
		         if np.ndim(f) == 0 else np.asarray(f, dtype=np.float64))
		return _bilinear_interp(wind_angles, wind_speeds, boat_speeds, a, f_col)

	polaire_func.v_max = float(boat_speeds.max())
	return polaire_func

# AFFICHAGE

def grille(a):
	for k in range(4):
		b = ((k + 1) / 4) * a
		x = [ b * cos(i * 2 * pi / 100) for i in range(101) ]
		y = [ b * sin(i * 2 * pi / 100) for i in range(101) ]
		plt.axis("equal")
		plt.plot(x, y, color='grey', lw='0.5')
	for k in range(8):
		ang = k * 2 * pi / 8
		x = [a * cos(ang), -a * cos(ang)]
		y = [a * sin(ang), -a * sin(ang)]
		plt.axis("equal")
		plt.plot(x, y, color='grey', lw='0.5')

def aff_polaire(P, l, n):
	grille(20)
	for f in l:
		X = []
		Y = []
		for k in range(n + 1):
			ang = k * 180 / n          # 0° (face au vent) → 180° (vent arrière)
			v = P([ang], [f])[0]
			X.append(v * sin(radians(ang)))
			Y.append(v * cos(radians(ang)))
		plt.plot(X, Y)
		plt.scatter(X, Y)  # symétrie tribord/bâbord
	plt.axis('equal')
	plt.show()

