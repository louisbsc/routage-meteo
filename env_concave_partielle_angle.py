import numpy as np
from scipy.spatial import cKDTree

import fonctions_utiles as f

# def intersection(A, B, C, D):
#     """
#     Vérifie l'intersection entre plusieurs segments.
    
#     A, B, C, D : arrays de shape (N, 2) représentant les points des segments
#     Segment1 = (A[i], B[i])
#     Segment2 = (C[i], D[i])
    
#     Retour : array booléen de longueur N
#     """
#     # Vecteurs
#     AB = B - A  # shape (N,2)
#     CD = D - C
#     AC = C - A
#     AD = D - A
#     CA = A - C
#     CB = B - C
    
#     # Calcul des produits croisés
#     def cross(u, v):
#         return u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0]
    
#     o1 = cross(AB, AC)
#     o2 = cross(AB, AD)
#     o3 = cross(CD, CA)
#     o4 = cross(CD, CB)
    
#     # Intersection stricte
#     intersect = (o1*o2 < 0) & (o3*o4 < 0)
    
#     # Cas colinéaire
#     colinear = (o1 == 0) & (((np.minimum(A[:,0], B[:,0]) <= C[:,0]) & (C[:,0] <= np.maximum(A[:,0], B[:,0])) &
#                              (np.minimum(A[:,1], B[:,1]) <= C[:,1]) & (C[:,1] <= np.maximum(A[:,1], B[:,1])) ) |
#                             ((np.minimum(A[:,0], B[:,0]) <= D[:,0]) & (D[:,0] <= np.maximum(A[:,0], B[:,0])) &
#                              (np.minimum(A[:,1], B[:,1]) <= D[:,1]) & (D[:,1] <= np.maximum(A[:,1], B[:,1])) ) )
    
#     intersect |= colinear
#     return intersect

# def intersection(A, B, C, D):
#     AB = B - A
#     AC = C - A
#     AD = D - A

#     CA = A - C
#     CB = B - C
#     CD = D - C  # C et D sont constants

#     o1 = AB[:,0]*AC[:,1] - AB[:,1]*AC[:,0]
#     o2 = AB[:,0]*AD[:,1] - AB[:,1]*AD[:,0]
#     o3 = CD[0]*CA[:,1] - CD[1]*CA[:,0]
#     o4 = CD[0]*CB[:,1] - CD[1]*CB[:,0]

#     return (o1*o2 < 0) & (o3*o4 < 0)

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

# def bouclage(l, p, atol=1e-8):
#     if len(l) < 3:
#         return False
#     l = l[:, :2]
#     p = p[:2]
#     if np.any(np.all(np.isclose(l[1:], p, atol=atol), axis=1)):
#         return True
#     A = l[:-2]
#     B = l[1:-1]
#     C = np.array([l[-1]] * len(A))
#     D = np.array([p] * len(A))
#     return intersection(A, B, C, D).any()

# def bouclage(l, p, atol=1e-8):
#     if len(l) < 3:
#         return False

#     l = np.asarray(l)[:, :2]   # conversion unique
#     p = p[:2]

#     if np.any(np.all(np.isclose(l[1:], p, atol=atol), axis=1)):
#         return True

#     A = l[:-2]
#     B = l[1:-1]

#     C = np.broadcast_to(l[-1], (len(A), 2))
#     D = np.broadcast_to(p, (len(A), 2))

#     return intersection(A, B, C, D).any()

def bouclage(l, p, atol=1e-8):
    if l.shape[0] < 3:
        return False

    l2 = l[:, :2]
    p2 = p[:2]

    # remplace np.isclose
    if np.any(np.all(np.abs(l2[1:] - p2) < atol, axis=1)):
        return True

    A = l2[:-2]
    B = l2[1:-1]

    # broadcasting naturel numpy


    # n = len(A)
    # C = np.broadcast_to(l2[-1], (n, 2))
    # D = np.broadcast_to(p2, (n, 2))

    return intersection(A, B, l2[-1], p2).any()

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


def enveloppe(N, r, p_dep, p_arr, ang, delta):

    p_dep = np.array(p_dep, dtype=float)
    p_arr = np.array(p_arr, dtype=float)
    dir = f.angle_direction(p_dep, p_arr)
    p1, p = points_les_plus_a(N, dir + ang, p_dep, delta)

    #voisins1 = pointsautour_np(p1, P, r)
    tree = cKDTree(N[:, :2])
    voisins1 = N[tree.query_ball_point(p1[:2], r)]

    #angles1 = np.abs(np.array([f.angleoriente3_np(p, p1, x) for x in voisins1]))
    angles1 = f.angle_oriente_positif(p, p1, voisins1)
    
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

    max_points = 10000   # borne large mais sûre
    l = np.empty((max_points, N.shape[1]))
    l[0] = p
    l[1] = p1
    l[2] = p2
    k = 3
    
    #l = [p, p1, p2]

    #while not np.array_equal(p2, l[0]): 
    angle_g = f.angle_oriente_negatif(l[1], p_dep, p_arr)
    angle_oriente = f.angle_oriente_negatif(l[1], p_dep, l[k-1])
    while angle_oriente < angle_g + ang or angle_oriente > 6:

        #voisins = pointsautour_np(p2, P, r)
        voisins = N[tree.query_ball_point(p2[:2], r)]
        #angles = np.abs(np.array([f.angleoriente3_np(p1, p2, x) for x in voisins]))
        angles = f.angle_oriente_positif(p1, p2, voisins)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles)
        candidats = np.where(np.abs(angles - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx = candidats[0]
        else:
            dists = np.linalg.norm(voisins[candidats, :2] - p2[:2], axis=1)
            idx = candidats[np.argmin(dists)]

        p3 = voisins[idx]

        if max_angle > np.pi:

            while bouclage(l[:k], p3):
            # while bouclage(np.array(l), p3):
                mask = np.ones(len(voisins), dtype=bool)
                mask[idx] = False
                voisins = voisins[mask]
                angles = angles[mask]
                #voisins = np.delete(voisins, idx, axis=0)
                #angles = np.delete(angles, idx)
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
        angle_oriente = f.angle_oriente_negatif(l[1], p_dep, l[k-1])

    return l[1:k]









