import numpy as np
import heapq


def isopoint(p, p_arr, t, dt, n, n1, V, P):
    x, y, index_iso, index_origine = p

    direction_vent, force = V(p, t)  # <- force au lieu de f
    caps = np.linspace(0, 360, n, endpoint=False)
    angs = (caps - direction_vent) % 360 - 180
    vitesses = P(angs, force)

    dx = vitesses * dt * np.cos(np.pi / 2 - np.radians(caps))
    dy = vitesses * dt * np.sin(np.pi / 2 - np.radians(caps))

    liste_index_iso = np.full(n, index_iso + 1, dtype=float)
    liste_index_origine = np.full(n, 0, dtype=float)

    points = np.column_stack([dx + x, dy + y, liste_index_iso, liste_index_origine])

    # --- pruning : garder n1 points les plus proches de l'arrivée ---
    px, py = p_arr[:2]
    dist2 = (points[:, 0] - px) ** 2 + (points[:, 1] - py) ** 2
    idx = np.argpartition(dist2, n1)[:n1]

    return points[idx]


def astar_continu(p_dep, p_arr, t0, dt, n, n1, V, P, e_arr,
                  V_moy=6.0, eps=0.2):

    time_list = np.array([t0], dtype=int)

    p_arr = np.array(p_arr, dtype=float)
    p_dep = np.array(p_dep, dtype=float)

    def dist_to_goal(x, y):
        return np.hypot(x - p_arr[0], y - p_arr[1])

    def key_xy(x, y):
        return (int(round(x / eps)), int(round(y / eps)))

    x0, y0 = float(p_dep[0]), float(p_dep[1])
    start_key = key_xy(x0, y0)

    g0 = 0.0
    h0 = dist_to_goal(x0, y0) / V_moy
    f0 = g0 + h0

    heap = []
    uid = 0
    heapq.heappush(heap, (f0, g0, uid, x0, y0, t0))
    uid += 1

    best_g = {start_key: g0}
    parent = {start_key: None}
    parent_pos = {start_key: (x0, y0)}

    while heap:

        fcur, gcur, _, x, y, t = heapq.heappop(heap)
        cur_key = key_xy(x, y)

        if gcur > best_g.get(cur_key, np.inf):
            continue

        if dist_to_goal(x, y) <= e_arr:
            k = cur_key
            route = []
            while k is not None:
                route.append(parent_pos[k])
                k = parent[k]
            route.reverse()
            return np.array(route), gcur, np.arange(t0, t0 + len(route) * dt, dt)

        p = np.array([x, y, 0.0, 0.0], dtype=float)
        voisins = isopoint(p, p_arr, t, dt, n, n1, V, P)

        for v in voisins:
            xv, yv = float(v[0]), float(v[1])
            tv = t + dt
            gv = gcur + dt

            kv = key_xy(xv, yv)

            if gv >= best_g.get(kv, np.inf):
                continue

            best_g[kv] = gv
            parent[kv] = cur_key
            parent_pos[kv] = (xv, yv)

            hv = dist_to_goal(xv, yv) / V_moy
            fv = gv + hv

            heapq.heappush(heap, (fv, gv, uid, xv, yv, tv))
            uid += 1

    return None, None, None
 


 