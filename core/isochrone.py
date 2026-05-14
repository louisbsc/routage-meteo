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
	
def _arc_resample(dx_dense, dy_dense, n):
	"""Re-échantillonne une courbe polaire fermée à n points équirépartis en longueur d'arc."""
	dx_cl = np.append(dx_dense, dx_dense[0])
	dy_cl = np.append(dy_dense, dy_dense[0])
	seg_len = np.hypot(np.diff(dx_cl), np.diff(dy_cl))
	arc = np.concatenate([[0], np.cumsum(seg_len)])
	targets = np.linspace(0, arc[-1], n, endpoint=False)
	return np.interp(targets, arc, dx_cl), np.interp(targets, arc, dy_cl)

def iso_point(p, t, dt, n, V, P, n_dense=180):
	x, y, index_iso, index_origine = p
	dir_vent, vit_vent = V(p, t)

	cap_dense = np.linspace(0, 360, n_dense, endpoint=False)
	ang_dense = ((dir_vent - cap_dense + 180) % 360) - 180
	vit_dense = P(ang_dense, vit_vent)
	dx_dense = vit_dense * dt * np.cos(np.pi/2 - np.radians(cap_dense))
	dy_dense = vit_dense * dt * np.sin(np.pi/2 - np.radians(cap_dense))

	dx, dy = _arc_resample(dx_dense, dy_dense, n)

	liste_index_iso = np.full(n, index_iso + 1, dtype=float)
	liste_index_origine = np.full(n, 0, dtype=float)
	points = np.column_stack([dx + x, dy + y, liste_index_iso, liste_index_origine])
	_, idx_unique = np.unique(np.round(points[:, :2], decimals=3), axis=0, return_index=True)
	return points[idx_unique]
	

def nuage_iso(I, t, dt, n, V, P, n_dense=180):
	M = len(I)

	wind = V(I, t)
	dir_vent = wind[:, 0]  # (M,)
	vit_vent = wind[:, 1]  # (M,)

	# Échantillonnage dense vectorisé sur M × n_dense
	cap_dense = np.linspace(0, 360, n_dense, endpoint=False)          # (n_dense,)
	cos_dense = np.cos(np.pi / 2 - np.radians(cap_dense))
	sin_dense = np.sin(np.pi / 2 - np.radians(cap_dense))

	ang_dense = (cap_dense[None, :] - dir_vent[:, None] + 180) % 360 - 180  # (M, n_dense)
	vit_dense = P(ang_dense.ravel(), np.repeat(vit_vent, n_dense)).reshape(M, n_dense)

	dx_dense = vit_dense * dt * cos_dense[None, :]  # (M, n_dense)
	dy_dense = vit_dense * dt * sin_dense[None, :]  # (M, n_dense)

	# Longueurs d'arc cumulées vectorisées : (M, n_dense+1)
	dx_cl = np.hstack([dx_dense, dx_dense[:, :1]])
	dy_cl = np.hstack([dy_dense, dy_dense[:, :1]])
	seg_len = np.hypot(np.diff(dx_cl, axis=1), np.diff(dy_cl, axis=1))
	arc = np.hstack([np.zeros((M, 1)), np.cumsum(seg_len, axis=1)])

	# Re-échantillonnage à n points par courbe — vectorisé (searchsorted row-offset)
	n_arc     = arc.shape[1]                                                         # n_dense + 1
	arc_total = arc[:, -1]                                                           # (M,)
	t_tgt     = arc_total[:, None] * np.linspace(0, 1, n, endpoint=False)[None, :]  # (M, n)
	scale     = float(arc_total.max()) + 1.0
	row_off   = np.arange(M, dtype=float) * scale
	idx_g = np.searchsorted(
		(arc   + row_off[:, None]).ravel(),
		(t_tgt + row_off[:, None]).ravel(),
		side='right',
	)
	loc   = np.clip(idx_g - np.repeat(np.arange(M) * n_arc, n), 1, n_arc - 1).reshape(M, n)
	rows  = np.arange(M)[:, None]
	t0_   = arc[rows, loc - 1];  t1_ = arc[rows, loc]
	alpha = (t_tgt - t0_) / np.maximum(t1_ - t0_, 1e-15)
	all_x = dx_cl[rows, loc - 1] + alpha * (dx_cl[rows, loc] - dx_cl[rows, loc - 1]) + I[:, 0:1]
	all_y = dy_cl[rows, loc - 1] + alpha * (dy_cl[rows, loc] - dy_cl[rows, loc - 1]) + I[:, 1:2]

	index_iso     = np.repeat(I[:, 2] + 1, n)
	index_origine = np.repeat(np.arange(M, dtype=float), n)

	points = np.column_stack([all_x.ravel(), all_y.ravel(), index_iso, index_origine])

	_, idx_unique = np.unique(np.round(points[:, :2], decimals=3), axis=0, return_index=True)
	return points[idx_unique]

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

def toutes_iso(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, delta, progress_cb=None):
	time_list = np.array([t], dtype=int)
	p_arr = np.array(p_arr, dtype=float)
	I0 = np.array([p_dep], dtype=float)
	if progress_cb is not None:
		_dist_total = max(float(np.linalg.norm(p_arr[:2] - np.array(p_dep[:2]))), 1e-6)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P)
	I = env.enveloppe(I, r, p_dep, p_arr, ang, delta)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	while not iso_est_arrive(I, p_arr, e_arr):
		t += dt
		ang -= dang
		time_list = np.append(time_list, t)
		I = iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, delta)
		L.append(I)
		print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
		if progress_cb is not None:
			dist_to_arr = float(np.min(f.distance_np(I, p_arr)))
			progress_cb(max(1, min(99, int((1 - dist_to_arr / _dist_total) * 100))))
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

def routage(p_dep, p_arr, t, dt, n, V, P, ang, dang, delta, progress_cb=None):
	# p_dep[1] = 360 - p_dep[1]
	# p_arr[1] = 360 - p_arr[1]
	p_dep = [p_dep[1] * 60 * 0.7, p_dep[0] * 60, 0, 0]
	p_arr = [p_arr[1] * 60 * 0.7, p_arr[0] * 60]

	# R = 3443.9184665
	# p_dep = [R * p_dep[1] * np.pi / 180, R * np.log(np.tan(np.pi/4 + p_dep[0] * np.pi / 360)), 0, 0]
	# p_arr = [R * p_arr[1] * np.pi / 180, R * np.log(np.tan(np.pi/4 + p_arr[0] * np.pi / 360))]

	e_arr = P.v_max * dt / 2
	r = P.v_max * dt * 2 * np.pi / n
	L, route, time_list = toutes_iso(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, delta, progress_cb=progress_cb)
	
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

	