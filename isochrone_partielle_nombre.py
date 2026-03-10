import numpy as np

import fonctions_utiles as f
import env_concave_partielle_nombre as env

def iso_point(p, t, dt, n, V, P):
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
	
def nuage_iso(I, t, dt, n, V, P):
	blocs = []
	for i, p in enumerate(I):
		l = iso_point(p, t, dt, n, V, P)
		l[:, -1] = np.full(n, i, dtype=float)
		blocs.append(l)
	return np.vstack(blocs)

def iso_suivante(p_dep, I, t, dt, n, n_i, V, P, r, dir, delta):
	L = nuage_iso(I, t, dt, n, V, P)
	return env.enveloppe(L, r, p_dep, dir, n_i, delta)

def n_iso(N, p_dep, t, dt, n, n_i, V, P, r, dir, delta):
	I0 = np.array([p_dep], dtype=float)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	for _ in range(N - 2):
		t += dt
		I = iso_suivante(p_dep, I, t, dt, n, n_i, V, P, r, dir, delta)
		L.append(I)
		print(f"nombre isochrones : {len(L) - 1}, temps écoulé : {t + dt:.2f} heures, nombre de points : {len(I)}")
	return np.vstack(L)

def iso_est_arrive(I, p_arr, e_arr):
	distances = f.distance_np(I, p_arr)
	mask = distances <= e_arr
	return mask.any()

def point_qui_est_arrive(l, p_arr, e_arr):
	distances = f.distance_np(l, p_arr)
	return l[np.argmin(distances)]

def toutes_iso(p_dep, p_arr, t, dt, n, n_i, V, P, e_arr, r, delta):
	dir = f.angle_direction(p_dep, p_arr)
	time_list = np.array([t], dtype=int)
	p_arr = np.array(p_arr, dtype=float)
	I0 = np.array([p_dep], dtype=float)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	while not iso_est_arrive(I, p_arr, e_arr):
		t += dt
		time_list = np.append(time_list, t)
		I = iso_suivante(p_dep, I, t, dt, n, n_i, V, P, r, dir, delta)
		L.append(I)
		print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	p_final = point_qui_est_arrive(I, p_arr, e_arr)
	L = np.vstack(L)
	route = [p_final]
	n_iter = int(p_final[2])
	for k in range(n_iter):
		p = L[L[:, 2] == int(p_final[2]) - 1][int(p_final[3])]
		route.append(p)
		p_final = p
	return L, np.array(route)[::-1], time_list