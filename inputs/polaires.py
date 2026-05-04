import matplotlib.pyplot as plt
from math import *
import numpy as np

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

# POLAIRES RÉELLES

import pandas as pd
from scipy.interpolate import RegularGridInterpolator

def polaire(path):
	df = pd.read_csv(path)
	wind_angles = df.iloc[1:, 0].astype(float).values
	wind_speeds = df.columns[1:].astype(float)
	boat_speeds = df.iloc[1:, 1:].astype(float).values
	interpolate = RegularGridInterpolator((wind_angles, wind_speeds), boat_speeds, bounds_error=False, fill_value=0)

	def polaire_func(ang, f):
		ang = np.asarray(ang)                  # transforme en tableau si nécessaire
		a = np.abs(ang % 360 - 180)            # transformation des angles
		pts = np.column_stack((a, np.full_like(a, f)))  # on garde f constant
		return interpolate(pts)

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
			ang = k * 360 / n
			v = P(ang,f)
			X.append( v * sin(radians(ang)))
			Y.append(- v * cos(radians(ang)))
		plt.plot(X, Y)
	plt.axis('equal')
	plt.show()

