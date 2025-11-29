from math import *
import numpy as np

import fonctions_utiles_np as f
import parametres_np as param
import env_concave_np as env


def isopoint_np(p, t, dt, n, V, P):
	x, y, index_iso, index_origine = p
	direction_vent, f = V(p, t)
	caps = np.linspace(0, 360, n, endpoint=False)
	angs = (caps - direction_vent) % 360 - 180
	vitesses = P(angs, f) 
	dx = vitesses * dt * np.cos(np.pi/2 - np.radians(caps))
	dy = vitesses * dt * np.sin(np.pi/2 - np.radians(caps))
	liste_index_iso = np.full(n, index_iso + 1, dtype=float)
	liste_index_origine = np.full(n, 0, dtype=float)
	return np.column_stack([dx + x, dy + y, liste_index_iso, liste_index_origine])
	
def nuageinfo_np(I, t, dt, n, V, P):
	blocs = []
	for i, p in enumerate(I):
		l = isopoint_np(p, t, dt, n, V, P)
		l[:, -1] = np.full(n, i, dtype=float)
		blocs.append(l)
	return np.vstack(blocs)

def isosuivante_np(I, t, dt, n, V, P):
	L = nuageinfo_np(I, t, dt, n, V, P)
	return env.enveloppe_np(L, param.r)

def niso_np(p, t, N, dt, n, V, P):
	I0 = np.array([p], dtype=float)
	I = isopoint_np(p, t, dt, n, V, P)
	L = [I0, I]
	for _ in range(N - 2):
		t += dt
		I = isosuivante_np(I, t, dt, n, V, P)
		L.append(I)
		print(f"nombre isochrones : {len(L)}, temps écoulé : {t:.2f} heures, nombre de points : {len(I)}")
	return np.vstack(L)


def estarrive_np(l, p_arr, e_arr):
	distances = f.distance_np(l, p_arr)
	mask = distances <= e_arr
	return mask.any()

def point_estarrive_np(l, p_arr, e_arr):
	distances = f.distance_np(l, p_arr)
	return l[np.argmin(distances)]



def toutesiso_np(p_dep, p_arr, t, dt, n, V, P, e_arr):
	p_arr = np.array(p_arr, dtype=float)
	I0 = np.array([p_dep], dtype=float)
	I = isopoint_np(p_dep, t, dt, n, V, P)
	L = [I0, I]
	while not estarrive_np(I, p_arr, e_arr):
		t += dt
		I = isosuivante_np(I, t, dt, n, V, P)
		L.append(I)
		print(f"nombre isochrones : {len(L)}, temps écoulé : {t:.2f} heures, nombre de points : {len(I)}")
	p_final = point_estarrive_np(I, p_arr, e_arr)
	L = np.vstack(L)
	route = [p_final]
	n_iter = int(p_final[2])
	for k in range(n_iter):
		p = L[L[:, 2] == int(p_final[2]) - 1][int(p_final[3])]
		route.append(p)
		p_final = p
	return L, np.array(route)
