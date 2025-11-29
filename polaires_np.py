import matplotlib.pyplot as plt
from math import *
import numpy as np

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






from scipy.interpolate import RegularGridInterpolator

Pvector = [[0,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00],
	[0,  	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00],
	[0,		0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00],
	[0,		0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00],
	[0,		2.70,	3.50,	4.00,	5.00,	5.20,	2.70,	2.30],
	[0,		3.90,	4.60,	5.50,	6.20,	6.60,	5.00,	3.80],
	[0,		4.70,	5.70,	6.50,	7.20,	7.40,	7.00,	6.20],
	[0,		5.00,	6.10,	7.10,	7.90,	8.10,	8.40,	8.00],
	[0,		5.20,	6.30,	7.80,	9.00,	9.70,	10.10,	9.60],
	[0,		5.20,	6.40,	8.30,	9.90,	11.20,	11.60,	10.60],
	[0,		5.20,	6.60,	8.30,	10.50,	12.80,	13.10,	11.60],
	[0,		4.90,	6.50,	8.40,	10.70,	13.00,	13.40,	12.40],
	[0,		4.50,	6.20,	8.40,	10.80,	13.10,	13.60,	13.00],
	[0,		3.90,	5.70,	8.40,	10.50,	13.30,	14.30,	13.20],
	[0,		3.30,	4.80,	7.40,	9.50,	13.40,	14.80,	13.10],
	[0,		2.60,	3.70,	6.00,	7.90,	11.40,	13.80,	12.10],
	[0,		2.20,	3.10,	5.00,	6.50,	8.60,	10.80,	10.60],
	[0,		2.00,	3.00,	4.50,	6.00,	8.00,	9.00,	10.00],
	[0,		2.00,	3.00,	4.50,	6.00,	8.00,	9.00,	10.00]]

wind_speeds = [0, 6, 8, 12, 16, 20, 25, 30]
wind_angles = np.linspace(0, 180, 19, dtype=int)
interpolate = RegularGridInterpolator((wind_angles, wind_speeds), Pvector)

def Vector(ang, f):
	return interpolate((abs(ang % 360 - 180), f)).item()

def Vector_vect(ang, f):
	ang = np.asarray(ang)                  # transforme en tableau si nécessaire
	a = np.abs(ang % 360 - 180)            # transformation des angles
	pts = np.column_stack((a, np.full_like(a, f)))  # on garde f constant
	return interpolate(pts) 