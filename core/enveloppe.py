import math

import numba
import numpy as np
from scipy.spatial import cKDTree

import core.utils as f


# ---------------------------------------------------------------------------
# Fonctions compilées Numba — remplacent les appels numpy dans la boucle
# chaude de `enveloppe`.  cache=True → compilation sauvegardée sur disque,
# pas de warmup lors des sessions suivantes.
# ---------------------------------------------------------------------------

@numba.njit(cache=True)
def bouclage(l, p, atol=1e-8):
    """l: (K, ≥2), p: (≥2,) — True si le segment l[-1]→p croise l'enveloppe."""
    n = l.shape[0]
    if n < 3:
        return False
    px, py = p[0], p[1]
    for i in range(1, n):
        if abs(l[i, 0] - px) < atol and abs(l[i, 1] - py) < atol:
            return True
    # Segment C→D : dernier point de l vers le candidat p
    cx, cy = l[n - 1, 0], l[n - 1, 1]
    dx, dy = px, py
    min_cx = cx if cx < dx else dx
    max_cx = cx if cx > dx else dx
    min_cy = cy if cy < dy else dy
    max_cy = cy if cy > dy else dy
    CDx = dx - cx
    CDy = dy - cy
    # Test intersection avec chaque segment l[i]→l[i+1] pour i = 0..n-3
    for i in range(n - 2):
        ax, ay = l[i, 0], l[i, 1]
        bx, by = l[i + 1, 0], l[i + 1, 1]
        # Préfiltre boîte englobante
        min_ax = ax if ax < bx else bx
        max_ax = ax if ax > bx else bx
        min_ay = ay if ay < by else by
        max_ay = ay if ay > by else by
        if max_ax < min_cx or min_ax > max_cx or max_ay < min_cy or min_ay > max_cy:
            continue
        ABx = bx - ax
        ABy = by - ay
        o1 = ABx * (cy - ay) - ABy * (cx - ax)
        o2 = ABx * (dy - ay) - ABy * (dx - ax)
        o3 = CDx * (ay - cy) - CDy * (ax - cx)
        o4 = CDx * (by - cy) - CDy * (bx - cx)
        if (o1 * o2 < 0) and (o3 * o4 < 0):
            return True
    return False


@numba.njit(cache=True)
def _angle_positif_batch(a, b, c):
    """a,b: (≥2,), c: (N, ≥2) → angles orientés positifs (N,).

    Équivalent JIT de angle_oriente_positif quand a et b sont des points 1D
    et c un tableau 2D de N points.
    """
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    norm_ba = math.sqrt(bax * bax + bay * bay)
    n = c.shape[0]
    out = np.empty(n)
    for i in range(n):
        bcx = c[i, 0] - b[0]
        bcy = c[i, 1] - b[1]
        norm_bc = math.sqrt(bcx * bcx + bcy * bcy)
        if norm_ba < 1e-15 or norm_bc < 1e-15:
            out[i] = 0.0
            continue
        dot = bax * bcx + bay * bcy
        det = bax * bcy - bay * bcx
        ang = math.atan2(det, dot)
        if abs(ang) < 1e-10:
            ang = 0.0
        if ang < 0.0:
            ang += 2.0 * math.pi
        out[i] = ang
    return out


@numba.njit(cache=True)
def _angle_negatif_scalar(a, b, c):
    """a,b,c: (≥2,) → angle orienté négatif scalaire.

    Équivalent JIT de angle_oriente_negatif quand tous les arguments sont
    des points 1D (cas exclusif dans la boucle de enveloppe).
    """
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    norm_ba = math.sqrt(bax * bax + bay * bay)
    bcx = c[0] - b[0]
    bcy = c[1] - b[1]
    norm_bc = math.sqrt(bcx * bcx + bcy * bcy)
    if norm_ba < 1e-15 or norm_bc < 1e-15:
        return 0.0
    dot = bax * bcx + bay * bcy
    det = bax * bcy - bay * bcx
    ang = math.atan2(-det, dot)
    if abs(ang) < 1e-10:
        ang = 0.0
    if ang < 0.0:
        ang += 2.0 * math.pi
    return ang


# ---------------------------------------------------------------------------
# Fonctions numpy conservées (rétrocompatibilité / appels hors enveloppe)
# ---------------------------------------------------------------------------

def intersection(A, B, C, D):
    # préfiltre bbox
    min_ax = np.minimum(A[:, 0], B[:, 0])
    max_ax = np.maximum(A[:, 0], B[:, 0])
    min_ay = np.minimum(A[:, 1], B[:, 1])
    max_ay = np.maximum(A[:, 1], B[:, 1])

    min_cx = min(C[0], D[0])
    max_cx = max(C[0], D[0])
    min_cy = min(C[1], D[1])
    max_cy = max(C[1], D[1])

    mask = (
        (max_ax >= min_cx) & (min_ax <= max_cx) &
        (max_ay >= min_cy) & (min_ay <= max_cy)
    )

    if not np.any(mask):
        return np.zeros(A.shape[0], dtype=bool)

    A2 = A[mask]
    B2 = B[mask]

    AB = B2 - A2
    AC = C - A2
    AD = D - A2

    CA = A2 - C
    CB = B2 - C
    CD = D - C

    o1 = AB[:, 0] * AC[:, 1] - AB[:, 1] * AC[:, 0]
    o2 = AB[:, 0] * AD[:, 1] - AB[:, 1] * AD[:, 0]
    o3 = CD[0] * CA[:, 1] - CD[1] * CA[:, 0]
    o4 = CD[0] * CB[:, 1] - CD[1] * CB[:, 0]

    inter2 = (o1 * o2 < 0) & (o3 * o4 < 0)

    out = np.zeros(A.shape[0], dtype=bool)
    out[mask] = inter2
    return out


def distance_points_droite(points, p0, theta):
    points = np.asarray(points, dtype=float)
    p0 = np.asarray(p0, dtype=float)

    n = np.array([-np.sin(theta), np.cos(theta)])

    d = (points[..., :2] - p0[:2]) @ n
    return np.abs(d)

def points_les_plus_a(N, ang, p_dep, delta):
    dist = distance_points_droite(N, p_dep, ang)
    mask = dist <= delta
    N_select = N[mask]

    if not np.any(mask):
        raise ValueError(f"Aucun point trouvé dans la bande delta.")

    projections = N_select[:, 0] * np.cos(ang) + N_select[:, 1] * np.sin(ang)
    idx_max = np.argmax(projections)
    p = N_select[idx_max]
    p_ref = p[:2] + np.array([np.cos(ang), np.sin(ang)])
    p_ref = np.concatenate([p_ref, [0.0, 0.0]])
    return p, p_ref


def enveloppe(N, r, p_dep, p_arr, ang, delta, I=None):

    p_dep = np.array(p_dep, dtype=float)
    p_arr = np.array(p_arr, dtype=float)
    dir = f.angle_direction(p_dep, p_arr)
    p1, p = points_les_plus_a(N, dir + ang, p_dep, delta)

    tree = cKDTree(N[:, :2])

    if I is not None:
        mask_inact = N[:, 3] < 0
        if mask_inact.any():
            inactive_pts    = N[mask_inact]
            inactive_lookup = {int(-pt[3]) - 1: pt for pt in inactive_pts}
            inactive_I_set  = set(inactive_lookup.keys())
        else:
            inactive_lookup = {}
            inactive_I_set  = set()
    else:
        inactive_lookup = {}
        inactive_I_set  = set()

    voisins1 = N[tree.query_ball_point(p1[:2], r)]

    angles1 = _angle_positif_batch(p, p1, voisins1)

    # gestion cas 2 angles égaux
    max_angle = np.max(angles1)
    eps = 1e-12
    candidats = np.where(np.abs(angles1 - max_angle) < eps)[0]
    if len(candidats) == 1:
        idx = candidats[0]
    else:
        dists = np.linalg.norm(voisins1[candidats, :2] - p1[:2], axis=1)
        idx = candidats[np.argmin(dists)]

    p2 = voisins1[idx]

    max_points = 100000   # borne large mais sûre
    l = np.empty((max_points, N.shape[1]))
    l[0] = p
    l[1] = p1
    l[2] = p2
    k = 3

    angle_g = _angle_negatif_scalar(l[1], p_dep, p_arr)
    angle_oriente = _angle_negatif_scalar(l[1], p_dep, l[k-1])
    while angle_oriente < angle_g + ang or angle_oriente > 6:

        voisins = N[tree.query_ball_point(p2[:2], r)]
        angles = _angle_positif_batch(p1, p2, voisins)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles)
        candidats = np.where(np.abs(angles - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx = candidats[0]
        else:
            dists = np.linalg.norm(voisins[candidats, :2] - p2[:2], axis=1)
            idx = candidats[np.argmin(dists)]

        p3 = voisins[idx]

        if p3[3] < 0 and inactive_I_set and not bouclage(l[:k], p3):
            j = int(-p3[3]) - 1
            chain = [p3]
            jf = j + 1
            while jf in inactive_I_set:
                chain.append(inactive_lookup[jf])
                jf += 1
            if len(chain) > 100:
                for pt in chain[:-100]:
                    l[k] = pt
                    k += 1
                p1 = l[k - 2]
                p2 = l[k - 1]
                angle_oriente = _angle_negatif_scalar(l[1], p_dep, p2)
                continue

        if max_angle > np.pi:

            while bouclage(l[:k], p3):
                mask = np.ones(len(voisins), dtype=bool)
                mask[idx] = False
                voisins = voisins[mask]
                angles = angles[mask]
                if voisins.size == 0:
                    raise ValueError(f"Aucun voisin trouvé dans le cercle r.")

                # gestion cas 2 angles égaux
                max_angle = np.max(angles)
                candidats = np.where(np.abs(angles - max_angle) < eps)[0]
                if len(candidats) == 1:
                    idx = candidats[0]
                else:
                    dists = np.linalg.norm(voisins[candidats, :2] - p2[:2], axis=1)
                    idx = candidats[np.argmin(dists)]

                p3 = voisins[idx]

        l[k] = p3
        k += 1
        p1, p2 = p2, p3
        angle_oriente = _angle_negatif_scalar(l[1], p_dep, l[k-1])

    return l[1:k]
