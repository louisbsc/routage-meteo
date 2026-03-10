import numpy as np
from scipy.spatial import cKDTree

import fonctions_utiles as f

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

def distance_points_droite(points, p0, theta):
    points = np.asarray(points, dtype=float)
    p0 = np.asarray(p0, dtype=float)

    # Normale unitaire à la droite
    n = np.array([-np.sin(theta), np.cos(theta)])

    d = (points[..., :2] - p0[:2]) @ n
    return np.abs(d)

def points_les_plus_a(N, ang, p_dep, delta):
    dist = distance_points_droite(N, p_dep, ang)
    mask = dist <= delta
    N_select = N[mask]

    projections = N_select[:, 0] * np.cos(ang) + N_select[:, 1] * np.sin(ang)
    idx_max = np.argmax(projections)
    p = N_select[idx_max]
    p_ref = p[:2] + np.array([np.cos(ang), np.sin(ang)])
    p_ref = np.concatenate([p_ref, [0.0, 0.0]])
    return p, p_ref


def cut(P, A, B, eps=1e-12, return_on_line=False):
    """
    Coupe un nuage de points selon la droite passant par A et B.

    P : array (N,2) ou (N,>=2)
    A : array (2,) ou (>=2,)
    B : array (2,) ou (>=2,)
    
    eps : tolérance numérique
    return_on_line : si True, renvoie aussi les points sur la droite

    Retour :
      - (P_gauche, P_droite) 
      ou (P_gauche, P_droite, P_sur_droite) si return_on_line=True
    """
    P = np.asarray(P, dtype=float)
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    # On travaille en 2D uniquement
    P2 = P[:, :2]
    A2 = A[:2]
    B2 = B[:2]

    AB = B2 - A2
    AP = P2 - A2

    # déterminant 2D = cross(AB, AP)
    det = AB[0] * AP[:, 1] - AB[1] * AP[:, 0]

    mask_gauche = det > eps
    mask_droite = det < -eps
    mask_ligne  = np.abs(det) <= eps

    P_gauche = P[mask_gauche]
    P_droite = P[mask_droite]

    if return_on_line:
        P_ligne = P[mask_ligne]
        return P_gauche, P_droite, P_ligne

    return P_gauche, P_droite

def cut_direction(P, A, theta, eps=1e-12, return_on_line=False):
    """
    Coupe un nuage de points selon la droite passant par A
    et orientée selon l'angle theta (radians).

    P : array (N,2) ou (N,>=2)
    A : array (2,) ou (>=2,)
    theta : float (radians)

    eps : tolérance numérique
    return_on_line : si True, renvoie aussi les points sur la droite

    Retour :
      - (P_gauche, P_droite)
      ou (P_gauche, P_droite, P_sur_droite)
    """
    P = np.asarray(P, dtype=float)
    A = np.asarray(A, dtype=float)

    # On travaille en 2D uniquement
    P2 = P[:, :2]
    A2 = A[:2]

    # vecteur directeur unitaire
    d = np.array([np.cos(theta), np.sin(theta)])

    AP = P2 - A2

    # déterminant 2D = cross(d, AP)
    det = d[0] * AP[:, 1] - d[1] * AP[:, 0]

    mask_gauche = det > eps
    mask_droite = det < -eps
    mask_ligne  = np.abs(det) <= eps

    P_gauche = P[mask_gauche]
    P_droite = P[mask_droite]

    if return_on_line:
        P_ligne = P[mask_ligne]
        return P_gauche, P_droite, P_ligne

    return P_gauche, P_droite

def stop(l, p_dep, ang_lim = np.pi/2):
    i = len(l) // 2
    return f.angle_oriente_negatif(l[i], p_dep, l[-2]) < ang_lim and f.angle_oriente_positif(l[i], p_dep, l[0]) < ang_lim


def enveloppe(N, r, p_dep, dir, n_i, delta):
    p_dep = np.array(p_dep, dtype=float)

    p1, p = points_les_plus_a(N, dir, p_dep, delta)

    #voisins1 = pointsautour_np(p1, P, r)
    tree = cKDTree(N[:, :2])
    voisins1 = N[tree.query_ball_point(p1[:2], r)]

    #angles1 = np.abs(np.array([f.angleoriente3_np(p, p1, x) for x in voisins1]))
    angles1_d = f.angle_oriente_positif(p, p1, voisins1)

    # gestion cas 2 angles égaux
    max_angle = np.max(angles1_d)
    eps = 1e-12
    candidats = np.where(np.abs(angles1_d - max_angle) < eps)[0]
    if len(candidats) == 1:
        idx1_d = candidats[0]
    else:
        dists = np.linalg.norm(voisins1[candidats, :2] - p1[:2], axis=1)
        idx1_d = candidats[np.argmin(dists)]

    p2_d = voisins1[idx1_d]


    angles1_g = f.angle_oriente_negatif(p, p1, voisins1)

    # gestion cas 2 angles égaux
    max_angle = np.max(angles1_g)
    eps = 1e-12
    candidats = np.where(np.abs(angles1_g - max_angle) < eps)[0]
    if len(candidats) == 1:
        idx1_g = candidats[0]
    else:
        dists = np.linalg.norm(voisins1[candidats, :2] - p1[:2], axis=1)
        idx1_g = candidats[np.argmin(dists)]

    p2_g = voisins1[idx1_g]

    l = [p2_g, p1, p2_d]

    p1_g, p1_d = p1, p1

    N_d, N_g = cut(N, p, p1)
    tree_d = cKDTree(N_d[:, :2])
    tree_g = cKDTree(N_g[:, :2])

    #while not np.array_equal(p2, l[0]):   
    while len(l) < n_i and stop(l, p_dep) and (f.distance_np(np.array([l[0]]), np.array([p1])) < r or f.distance_np(np.array([l[-1]]), np.array([p1])) < r):

        #voisins = pointsautour_np(p2, P, r)
        voisins_d = N_d[tree_d.query_ball_point(p2_d[:2], r)]
        #angles = np.abs(np.array([f.angleoriente3_np(p1, p2, x) for x in voisins]))
        angles_d = f.angle_oriente_positif(p1_d, p2_d, voisins_d)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles_d)
        candidats = np.where(np.abs(angles_d - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx_d = candidats[0]
        else:
            dists = np.linalg.norm(voisins_d[candidats, :2] - p2_d[:2], axis=1)
            idx_d = candidats[np.argmin(dists)]

        p3_d = voisins_d[idx_d]

        #while bouclage(np.array(l), p3):
        while bouclage(np.array(l), p3_d):
            mask = np.ones(len(voisins_d), dtype=bool)
            mask[idx_d] = False
            voisins_d = voisins_d[mask]
            angles_d = angles_d[mask]
            #voisins_d = np.delete(voisins_d, idx, axis=0)
            #angles = np.delete(angles, idx)
            if voisins_d.size == 0:
                print("voisins_d vide")

            # gestion cas 2 angles égaux
            max_angle = np.max(angles_d)
            candidats = np.where(np.abs(angles_d - max_angle) < eps)[0]
            if len(candidats) == 1:
                idx_d = candidats[0]
            else:
                dists = np.linalg.norm(voisins_d[candidats, :2] - p2_d[:2], axis=1)
                idx_d = candidats[np.argmin(dists)]

            p3_d = voisins_d[idx_d]
        l.append(p3_d)



        voisins_g = N_g[tree_g.query_ball_point(p2_g[:2], r)]
        angles_g = f.angle_oriente_negatif(p1_g, p2_g, voisins_g)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles_g)
        candidats = np.where(np.abs(angles_g - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx_g = candidats[0]
        else:
            dists = np.linalg.norm(voisins_g[candidats, :2] - p2_g[:2], axis=1)
            idx_g = candidats[np.argmin(dists)]

        p3_g = voisins_g[idx_g]

        while bouclage(np.array(l)[::-1], p3_g):
            mask = np.ones(len(voisins_g), dtype=bool)
            mask[idx_g] = False
            voisins_g = voisins_g[mask]
            angles_g = angles_g[mask]
            if voisins_g.size == 0:
                print("voisins_g vide")

            # gestion cas 2 angles égaux
            max_angle = np.max(angles_g)
            candidats = np.where(np.abs(angles_g - max_angle) < eps)[0]
            if len(candidats) == 1:
                idx_g = candidats[0]
            else:
                dists = np.linalg.norm(voisins_g[candidats, :2] - p2_g[:2], axis=1)
                idx_g = candidats[np.argmin(dists)]

            p3_g = voisins_g[idx_g]

        l = [p3_g] + l

        p1_g, p1_d = p2_g, p2_d
        p2_g, p2_d = p3_g, p3_d


    while len(l) < n_i and stop(l, p_dep):

        #voisins = pointsautour_np(p2, P, r)
        voisins_d = N[tree.query_ball_point(p2_d[:2], r)]
        #angles = np.abs(np.array([f.angleoriente3_np(p1, p2, x) for x in voisins]))
        angles_d = f.angle_oriente_positif(p1_d, p2_d, voisins_d)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles_d)
        candidats = np.where(np.abs(angles_d - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx_d = candidats[0]
        else:
            dists = np.linalg.norm(voisins_d[candidats, :2] - p2_d[:2], axis=1)
            idx_d = candidats[np.argmin(dists)]

        p3_d = voisins_d[idx_d]

        #while bouclage(np.array(l), p3):
        while bouclage(np.array(l), p3_d):
            mask = np.ones(len(voisins_d), dtype=bool)
            mask[idx_d] = False
            voisins_d = voisins_d[mask]
            angles_d = angles_d[mask]
            #voisins_d = np.delete(voisins_d, idx, axis=0)
            #angles = np.delete(angles, idx)
            if voisins_d.size == 0:
                print("voisins_d vide")

            # gestion cas 2 angles égaux
            max_angle = np.max(angles_d)
            candidats = np.where(np.abs(angles_d - max_angle) < eps)[0]
            if len(candidats) == 1:
                idx_d = candidats[0]
            else:
                dists = np.linalg.norm(voisins_d[candidats, :2] - p2_d[:2], axis=1)
                idx_d = candidats[np.argmin(dists)]

            p3_d = voisins_d[idx_d]
        l.append(p3_d)



        voisins_g = N[tree.query_ball_point(p2_g[:2], r)]
        angles_g = f.angle_oriente_negatif(p1_g, p2_g, voisins_g)

        # gestion cas 2 angles égaux
        max_angle = np.max(angles_g)
        candidats = np.where(np.abs(angles_g - max_angle) < eps)[0]
        if len(candidats) == 1:
            idx_g = candidats[0]
        else:
            dists = np.linalg.norm(voisins_g[candidats, :2] - p2_g[:2], axis=1)
            idx_g = candidats[np.argmin(dists)]

        p3_g = voisins_g[idx_g]

        while bouclage(np.array(l)[::-1], p3_g):
            mask = np.ones(len(voisins_g), dtype=bool)
            mask[idx_g] = False
            voisins_g = voisins_g[mask]
            angles_g = angles_g[mask]
            if voisins_g.size == 0:
                print("voisins_g vide")

            # gestion cas 2 angles égaux
            max_angle = np.max(angles_g)
            candidats = np.where(np.abs(angles_g - max_angle) < eps)[0]
            if len(candidats) == 1:
                idx_g = candidats[0]
            else:
                dists = np.linalg.norm(voisins_g[candidats, :2] - p2_g[:2], axis=1)
                idx_g = candidats[np.argmin(dists)]

            p3_g = voisins_g[idx_g]

        l = [p3_g] + l

        p1_g, p1_d = p2_g, p2_d
        p2_g, p2_d = p3_g, p3_d

    return np.array(l)









