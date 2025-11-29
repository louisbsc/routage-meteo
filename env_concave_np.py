import numpy as np
from scipy.spatial import cKDTree

import fonctions_utiles_np as f

def pointsautour_np(p0, P, r):
    diff = P - p0
    distances = np.hypot(diff[:, 0], diff[:, 1])
    mask = (distances < r) & (distances > 0)
    return P[mask]

def pointsautour_np_2(p0, P, r):
    tree = cKDTree(P[:, :2])
    idx_voisins = tree.query_ball_point(p0[:2], r)
    return P[idx_voisins]

def intersec_np(seg1, seg2):
    p1, p2 = seg1
    p3, p4 = seg2

    if (p1 == p4).all():
        return False

    def orientation(a, b, c):
        return (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])

    o1 = orientation(p1, p2, p3)
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)

    if o1 * o2 < 0 and o3 * o4 < 0:
        return True

    # # Fonction pour tester chevauchement strict
    # def overlap_strict(a1, a2, b1, b2):
    #     """Retourne True si les intervalles [a1,a2] et [b1,b2]
    #        se chevauchent en plus qu’un seul point."""
    #     lo = max(min(a1, a2), min(b1, b2))
    #     hi = min(max(a1, a2), max(b1, b2))
    #     return hi > lo  # strictement > → exclut un point seul

    # # --- CAS 2 : Colinéaires
    # if o1 == 0 and o2 == 0 and o3 == 0 and o4 == 0:
    #     # Colinéaires : tester chevauchement strict en x et y
    #     ox = overlap_strict(p1[0], p2[0], p3[0], p4[0])
    #     oy = overlap_strict(p1[1], p2[1], p3[1], p4[1])
    #     return ox or oy


    def on_segment(p, q, r):
        return (min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and min(p[1], r[1]) <= q[1] <= max(p[1], r[1]))

    if o1 == 0 and on_segment(p1, p3, p2): return True
    if o2 == 0 and on_segment(p1, p4, p2): return True
    if o3 == 0 and on_segment(p3, p1, p4): return True
    if o4 == 0 and on_segment(p3, p2, p4): return True

    return False


def boucle_np(l, p):
    n = len(l)
    seg  = [l[-1], p]
    for k in range(n - 2):
        if intersec_np([l[k], l[k+1]], seg):
            return True
    return False






def segments_intersect_numpy(A, B, C, D):
    """
    Vérifie l'intersection entre plusieurs segments.
    
    A, B, C, D : arrays de shape (N, 2) représentant les points des segments
    Segment1 = (A[i], B[i])
    Segment2 = (C[i], D[i])
    
    Retour : array booléen de longueur N
    """
    # Vecteurs
    AB = B - A  # shape (N,2)
    CD = D - C
    AC = C - A
    AD = D - A
    CA = A - C
    CB = B - C
    
    # Calcul des produits croisés
    def cross(u, v):
        return u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0]
    
    o1 = cross(AB, AC)
    o2 = cross(AB, AD)
    o3 = cross(CD, CA)
    o4 = cross(CD, CB)
    
    # Intersection stricte
    intersect = (o1*o2 < 0) & (o3*o4 < 0)
    
    # Cas colinéaire
    colinear = (o1 == 0) & (((np.minimum(A[:,0], B[:,0]) <= C[:,0]) & (C[:,0] <= np.maximum(A[:,0], B[:,0])) &
                             (np.minimum(A[:,1], B[:,1]) <= C[:,1]) & (C[:,1] <= np.maximum(A[:,1], B[:,1])) ) |
                            ((np.minimum(A[:,0], B[:,0]) <= D[:,0]) & (D[:,0] <= np.maximum(A[:,0], B[:,0])) &
                             (np.minimum(A[:,1], B[:,1]) <= D[:,1]) & (D[:,1] <= np.maximum(A[:,1], B[:,1])) ) )
    
    intersect |= colinear
    return intersect

def boucle_np_vect(l, p):
    if len(l) < 3:
        return False
    l = l[:, :2]
    p = p[:2]
    A = l[:-2]
    B = l[1:-1]
    C = np.array([l[-1]] * len(A))
    D = np.array([p] * len(A))
    return segments_intersect_numpy(A, B, C, D).any()





def enveloppe_np(P, r):
    p1 = P[np.argmin(P[:, 0])]

    p = np.array([p1[0], p1[1] - 1])

    tree = cKDTree(P[:, :2])

    #voisins1 = pointsautour_np(p1, P, r)
    voisins1 = P[tree.query_ball_point(p1[:2], r)]

    #angles1 = np.abs(np.array([f.angleoriente3_np(p, p1, x) for x in voisins1]))
    angles1 = f.angleoriente3_np_vect(p, p1, voisins1)
    p2 = voisins1[np.argmax(angles1)]

    l = [p1, p2]

    #while not np.array_equal(p2, l[0]):
    while not np.array_equal(p2[:2], l[0][:2]):

        #voisins = pointsautour_np(p2, P, r)
        voisins = P[tree.query_ball_point(p2[:2], r)]
        #angles = np.abs(np.array([f.angleoriente3_np(p1, p2, x) for x in voisins]))
        angles = f.angleoriente3_np_vect(p1, p2, voisins)
        idx = np.argmax(angles)
        p3 = voisins[idx]

        #while boucle_np_vect(np.array(l), p3):
        while boucle_np_vect(np.array(l), p3):
            mask = np.ones(len(voisins), dtype=bool)
            mask[idx] = False
            voisins = voisins[mask]
            angles = angles[mask]
            #voisins = np.delete(voisins, idx, axis=0)
            #angles = np.delete(angles, idx)
            if voisins.size == 0:
                break
            idx = np.argmax(angles)
            p3 = voisins[idx]

        l.append(p3)
        p1, p2 = p2, p3

    return np.array(l)









