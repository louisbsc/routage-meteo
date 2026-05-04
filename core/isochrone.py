import numpy as np

import core.utils as f
import core.enveloppe as env

import inputs.courants

# def iso_point(p, t, dt, n, V, P, C = courants.C0):
# 	x, y, index_iso, index_origine = p
# 	dir_vent, vit_vent = V(p, t)
# 	dir_courant, vit_courant = C(p, t)
# 	dx_courant = vit_courant * np.cos(np.pi/2 - np.radians(dir_courant))
# 	dy_courant = vit_courant * np.sin(np.pi/2 - np.radians(dir_courant))
# 	cap = np.linspace(0, 360, n, endpoint=False)
# 	ang_au_vent = (cap - dir_vent) % 360 - 180
# 	vit_bateau = P(ang_au_vent, vit_vent)
# 	dx = vit_bateau * dt * np.cos(np.pi/2 - np.radians(cap)) + dx_courant
# 	dy = vit_bateau * dt * np.sin(np.pi/2 - np.radians(cap)) + dy_courant
# 	liste_index_iso = np.full(n, index_iso + 1, dtype=float)
# 	liste_index_origine = np.full(n, 0, dtype=float)
# 	return np.column_stack([dx + x, dy + y, liste_index_iso, liste_index_origine])
	
def iso_point(p, t, dt, n, V, P):
	x, y, index_iso, index_origine = p
	dir_vent, vit_vent = V(p, t)
	cap = np.linspace(0, 360, n, endpoint=False)
	ang_au_vent = (cap - dir_vent) % 360 - 180
	vit_bateau = P(ang_au_vent, vit_vent)
	dx = vit_bateau * dt * np.cos(np.pi/2 - np.radians(cap))
	dy = vit_bateau * dt * np.sin(np.pi/2 - np.radians(cap))
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

def iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, delta):
	L = nuage_iso(I, t, dt, n, V, P)
	return env.enveloppe(L, r, p_dep, p_arr, ang, delta)

def n_iso(N, p_dep, p_arr, t, dt, n, V, P, r, ang, dang, delta):
	I0 = np.array([p_dep], dtype=float)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	for _ in range(N - 1):
		t += dt
		ang -= dang
		I = iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, delta)
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

def toutes_iso(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, delta):
	time_list = np.array([t], dtype=int)
	p_arr = np.array(p_arr, dtype=float)
	I0 = np.array([p_dep], dtype=float)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	while not iso_est_arrive(I, p_arr, e_arr):
		t += dt
		ang -= dang
		time_list = np.append(time_list, t)
		I = iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, delta)
		L.append(I)
		print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	p_final = point_qui_est_arrive(I, p_arr, e_arr)
	time_list = np.append(time_list, t + dt)
	L = np.vstack(L)
	route = [p_final]
	n_iter = int(p_final[2])
	for k in range(n_iter):
		p = L[L[:, 2] == int(p_final[2]) - 1][int(p_final[3])]
		route.append(p)
		p_final = p
	return L, np.array(route)[::-1], time_list

def routage(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, delta):
	# p_dep[1] = 360 - p_dep[1]
	# p_arr[1] = 360 - p_arr[1]
	p_dep = [p_dep[1] * 60 * 0.7, p_dep[0] * 60, 0, 0]
	p_arr = [p_arr[1] * 60 * 0.7, p_arr[0] * 60]

	# R = 3443.9184665
	# p_dep = [R * p_dep[1] * np.pi / 180, R * np.log(np.tan(np.pi/4 + p_dep[0] * np.pi / 360)), 0, 0]
	# p_arr = [R * p_arr[1] * np.pi / 180, R * np.log(np.tan(np.pi/4 + p_arr[0] * np.pi / 360))]

	
	L, route, time_list = toutes_iso(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, delta)
	
	latitude = route[:, 1] / 60
	longitude = route[:, 0] / (60 * 0.7)

	# x = route[:, 0]
	# y = route[:, 1]
	# lon_rad = x / R
	# lat_rad = 2 * np.arctan(np.exp(y / R)) - np.pi / 2
	# longitude = np.degrees(lon_rad)
	# latitude = np.degrees(lat_rad)

	
	L[:, 0], L[:, 1] = L[:, 1] / 60, L[:, 0] / (60 * 0.7)
	# x = L[:, 0]
	# y = L[:, 1]
	# lon_rad = x / R
	# lat_rad = 2 * np.arctan(np.exp(y / R)) - np.pi / 2
	# L[:, 0] = np.degrees(lon_rad)
	# L[:, 1] = np.degrees(lat_rad)

	return latitude, longitude, time_list, L

	