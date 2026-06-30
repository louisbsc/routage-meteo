import math

import numba
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

# Cache des tableaux trigonométriques pour éviter de les recalculer à chaque appel
_cap_cache: dict = {}

def _polar_offsets_batch(dir_vents, vit_vents, dt, P, n_dense=90):
	"""Déplacements (dx, dy) pour M points × n_dense caps. Retourne (M, n_dense), (M, n_dense)."""
	if n_dense not in _cap_cache:
		cap_dense = np.linspace(0, 360, n_dense, endpoint=False)
		_cap_cache[n_dense] = (
			cap_dense,
			np.cos(np.pi / 2 - np.radians(cap_dense)),
			np.sin(np.pi / 2 - np.radians(cap_dense)),
		)
	cap_dense, cos_cap, sin_cap = _cap_cache[n_dense]
	ang = (cap_dense[None, :] - dir_vents[:, None] + 180) % 360 - 180  # (M, n_dense)
	vit = P(ang.ravel(), np.repeat(vit_vents, n_dense)).reshape(len(dir_vents), n_dense)
	return vit * dt * cos_cap[None, :], vit * dt * sin_cap[None, :]


@numba.njit(cache=True, parallel=True)
def _arc_resample_batch(dx_dense, dy_dense, n):
	"""Ré-échantillonne M courbes polaires fermées à n points équirépartis en arc.

	Version compilée Numba parallèle : chaque courbe traitée indépendamment
	sur un cœur distinct via prange.
	Entrée : (M, n_dense). Retourne all_x, all_y : (M, n).
	"""
	M       = dx_dense.shape[0]
	n_dense = dx_dense.shape[1]
	all_x   = np.empty((M, n))
	all_y   = np.empty((M, n))

	for mi in numba.prange(M):
		arc = np.empty(n_dense + 1)   # alloué par thread via prange
		# Arc-length parameterization (courbe fermée : dernier → premier)
		arc[0] = 0.0
		for j in range(n_dense - 1):
			ddx = dx_dense[mi, j + 1] - dx_dense[mi, j]
			ddy = dy_dense[mi, j + 1] - dy_dense[mi, j]
			arc[j + 1] = arc[j] + math.sqrt(ddx * ddx + ddy * ddy)
		ddx = dx_dense[mi, 0] - dx_dense[mi, n_dense - 1]
		ddy = dy_dense[mi, 0] - dy_dense[mi, n_dense - 1]
		arc[n_dense] = arc[n_dense - 1] + math.sqrt(ddx * ddx + ddy * ddy)

		arc_total = arc[n_dense]

		if arc_total < 1e-15:
			for k in range(n):
				all_x[mi, k] = dx_dense[mi, 0]
				all_y[mi, k] = dy_dense[mi, 0]
			continue

		# Rééchantillonnage à n points équirépartis en arc
		for k in range(n):
			t_k = arc_total * k / n
			# Recherche binaire : arc[lo] <= t_k < arc[hi]
			lo, hi = 0, n_dense
			while lo < hi - 1:
				mid = (lo + hi) >> 1
				if arc[mid] <= t_k:
					lo = mid
				else:
					hi = mid
			# Interpolation linéaire entre lo et hi
			dt_ = arc[hi] - arc[lo]
			alpha = 0.0 if dt_ < 1e-15 else (t_k - arc[lo]) / dt_
			x0 = dx_dense[mi, lo]
			y0 = dy_dense[mi, lo]
			# hi == n_dense correspond à la fermeture (retour au point 0)
			if hi < n_dense:
				x1 = dx_dense[mi, hi]
				y1 = dy_dense[mi, hi]
			else:
				x1 = dx_dense[mi, 0]
				y1 = dy_dense[mi, 0]
			all_x[mi, k] = x0 + alpha * (x1 - x0)
			all_y[mi, k] = y0 + alpha * (y1 - y0)

	return all_x, all_y


def iso_point(p, t, dt, n, V, P, C=None, n_dense=90):
	return nuage_iso(np.array([p]), t, dt, n, V, P, C=C, n_dense=n_dense)


def nuage_iso(I, t, dt, n, V, P, C=None, n_dense=90):
	active_mask   = I[:, 3] >= 0
	inactive_mask = ~active_mask
	parts = []

	if active_mask.any():
		idx_act = np.where(active_mask)[0]
		I_act   = I[idx_act]
		wind    = V(I_act, t)

		if C is not None:
			current = C(I_act, t)
			d_w = np.radians(wind[:, 0])
			d_c = np.radians(current[:, 0])
			u_wind = -wind[:, 1] * np.sin(d_w)
			v_wind = -wind[:, 1] * np.cos(d_w)
			u_cur  = -current[:, 1] * np.sin(d_c)
			v_cur  = -current[:, 1] * np.cos(d_c)
			u_surf = u_wind - u_cur
			v_surf = v_wind - v_cur
			surf_speed = np.sqrt(u_surf**2 + v_surf**2)
			surf_dir   = (np.degrees(np.arctan2(u_surf, v_surf)) + 180) % 360
			dx_cur = u_cur * dt
			dy_cur = v_cur * dt
		else:
			surf_speed = wind[:, 1]
			surf_dir   = wind[:, 0]
			dx_cur = np.zeros(len(idx_act))
			dy_cur = np.zeros(len(idx_act))

		dx_dense, dy_dense = _polar_offsets_batch(surf_dir, surf_speed, dt, P, n_dense)
		dx_dense += dx_cur[:, None]
		dy_dense += dy_cur[:, None]
		all_x, all_y = _arc_resample_batch(dx_dense, dy_dense, n)
		all_x += I_act[:, 0:1]
		all_y += I_act[:, 1:2]
		index_iso     = np.repeat(I_act[:, 2] + 1, n)
		index_origine = np.repeat(idx_act.astype(float), n)
		parts.append(np.column_stack([all_x.ravel(), all_y.ravel(), index_iso, index_origine]))

	if inactive_mask.any():
		idx_inact = np.where(inactive_mask)[0]
		parts.append(np.column_stack([
			I[idx_inact, 0], I[idx_inact, 1],
			I[idx_inact, 2] + 1,
			-(idx_inact + 1).astype(float),
		]))

	points = np.vstack(parts)
	_, idx_unique = np.unique(np.round(points[:, :2], decimals=3), axis=0, return_index=True)
	return points[idx_unique]

def iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, C=None):
	L     = nuage_iso(I, t, dt, n, V, P, C=C)
	I_new = env.enveloppe(L, r, p_dep, p_arr, ang, r, I=I)
	wind_new = V(I_new, t + dt)
	inact = np.where(wind_new[:, 1] == 0.0)[0]
	if inact.size:
		I_new[inact, 3] = -(inact + 1).astype(float)
	return I_new

def n_iso(N, p_dep, p_arr, t, dt, n, V, P, r, ang, dang, C=None):
	I0 = np.array([p_dep], dtype=float)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P, C=C)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	for _ in range(N - 1):
		t += dt
		ang -= dang
		I = iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, C=C)
		L.append(I)
		print(f"nombre isochrones : {len(L) - 1}, temps écoulé : {t + dt:.2f} heures, nombre de points : {len(I)}")
	return np.vstack(L)

def iso_est_arrive(I, p_arr, e_arr):
	actifs = I[I[:, 3] >= 0]
	if actifs.shape[0] == 0:
		return False
	distances = f.distance_np(actifs, p_arr)
	mask = distances <= e_arr
	return mask.any()

def point_qui_est_arrive(l, p_arr, e_arr):
	actifs = l[l[:, 3] >= 0]
	distances = f.distance_np(actifs, p_arr)
	return actifs[np.argmin(distances)]

def toutes_iso(p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang, C=None, progress_cb=None):
	time_list = np.array([t], dtype=int)
	p_arr = np.array(p_arr, dtype=float)
	I0 = np.array([p_dep], dtype=float)
	if progress_cb is not None:
		_dist_total = max(float(np.linalg.norm(p_arr[:2] - np.array(p_dep[:2]))), 1e-6)
	print(f"nombre isochrones : 0, temps : {t:.2f} heures, nombre de points : {len(I0)}")
	I = iso_point(p_dep, t, dt, n, V, P, C=C)
	I = env.enveloppe(I, r, p_dep, p_arr, ang, r)
	wind_first = V(I, t + dt)
	inact_first = np.where(wind_first[:, 1] == 0.0)[0]
	if inact_first.size:
		I[inact_first, 3] = -(inact_first + 1).astype(float)
	L = [I0, I]
	print(f"nombre isochrones : {len(L) - 1}, temps : {t + dt:.2f} heures, nombre de points : {len(I)}")
	while not iso_est_arrive(I, p_arr, e_arr):
		t += dt
		ang -= dang
		time_list = np.append(time_list, t)
		I = iso_suivante(I, p_dep, p_arr, t, dt, n, V, P, r, ang, C=C)
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

def routage(p_dep, p_arr, t, dt, n, V, P, ang, dang, C=None, progress_cb=None):
	p_dep = [p_dep[1] * 60 * 0.7, p_dep[0] * 60, 0, 0]
	p_arr = [p_arr[1] * 60 * 0.7, p_arr[0] * 60]

	e_arr = P.v_max * dt / 2
	r = P.v_max * dt * 2 * np.pi / n
	L, route, time_list = toutes_iso(
		p_dep, p_arr, t, dt, n, V, P, e_arr, r, ang, dang,
		C=C, progress_cb=progress_cb,
	)

	latitude  = route[:, 1] / 60
	longitude = route[:, 0] / (60 * 0.7)

	L[:, 0], L[:, 1] = L[:, 1] / 60, L[:, 0] / (60 * 0.7)

	return latitude, longitude, time_list, L

def routage_raffine(route_lat, route_lon, p_dep, p_arr, t, dt, n, V, P, ang, seuil, facteur_raf, C=None, progress_cb=None):
	"""Raffinement d'une route existante : vent réduit au corridor, dt/facteur_raf, dang=0."""
	from inputs.vents import vent_filtre_route

	route_nm = np.column_stack([route_lon * 60 * 0.7, route_lat * 60])
	V_red = vent_filtre_route(V, route_nm, seuil)

	return routage(p_dep, p_arr, t, dt / facteur_raf, n, V_red, P, ang, 0, C=C, progress_cb=progress_cb)

