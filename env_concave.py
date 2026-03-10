import numpy as np
from scipy.spatial import cKDTree

import fonctions_utiles as f

# verifie l'intersection entre plusieurs segments
def intersection(A, B, C, D):
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

def bouclage(l, p, atol=1e-8):
    if len(l) < 3:
        return False
    l = l[:, :2]
    p = p[:2]
    if np.any(np.all(np.isclose(l[1:], p, atol=atol), axis=1)):
        return True
    A = l[:-2]
    B = l[1:-1]
    C = np.array([l[-1]] * len(A))
    D = np.array([p] * len(A))
    return intersection(A, B, C, D).any()

def enveloppe(P, r):
    p1 = P[np.argmin(P[:, 0])]

    p = np.array([p1[0], p1[1] - 1])

    tree = cKDTree(P[:, :2])

    #voisins1 = pointsautour_np(p1, P, r)
    voisins1 = P[tree.query_ball_point(p1[:2], r)]

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

    l = [p1, p2]

    #while not np.array_equal(p2, l[0]):
    while not np.array_equal(p2[:2], l[0][:2]):

        #voisins = pointsautour_np(p2, P, r)
        voisins = P[tree.query_ball_point(p2[:2], r)]
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

        #while bouclage(np.array(l), p3):
        while bouclage(np.array(l), p3):
            mask = np.ones(len(voisins), dtype=bool)
            mask[idx] = False
            voisins = voisins[mask]
            angles = angles[mask]
            #voisins = np.delete(voisins, idx, axis=0)
            #angles = np.delete(angles, idx)
            if voisins.size == 0:
                print("voisins vide")

            max_angle = np.max(angles)
            candidats = np.where(np.abs(angles - max_angle) < eps)[0]
            if len(candidats) == 1:
                idx = candidats[0]
            else:
                dists = np.linalg.norm(voisins[candidats, :2] - p2[:2], axis=1)
                idx = candidats[np.argmin(dists)]

            p3 = voisins[idx]

        l.append(p3)
        p1, p2 = p2, p3

    return np.array(l)







