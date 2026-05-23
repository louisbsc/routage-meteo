import numpy as np

def modulo_np(ang):
    return np.mod(ang, 2 * np.pi)

def distance_np(p0, p1):
    return np.linalg.norm(p1[..., :2] - p0[..., :2], axis=1)

def vecteur_np(p0, p1):
    return np.array(p1) - np.array(p0)

def multiplication_np(a, u):
    return np.array(u) * a

def produitscalaire_np(u, v):
    return np.dot(u, v)

def angle_direction(u, v):
    x1, y1 = u[:2]
    x2, y2 = v[:2]
    dx = x2 - x1
    dy = y2 - y1
    return np.atan2(dy, dx) % (2 * np.pi)

def angle_oriente_positif(a, b, c):
    ba = a[..., :2] - b[..., :2]
    bc = c[..., :2] - b[..., :2]

    zero_mask = (np.linalg.norm(ba, axis=-1) < 1e-15) | \
                (np.linalg.norm(bc, axis=-1) < 1e-15)

    dot = np.sum(ba * bc, axis=-1)
    det = ba[..., 0] * bc[..., 1] - ba[..., 1] * bc[..., 0]

    ang = np.arctan2(det, dot)
    
    ang = np.where(np.abs(ang) < 1e-10, 0.0, ang)

    ang = np.where(ang < 0, ang + 2*np.pi, ang)

    ang = np.where(zero_mask, 0.0, ang)

    return ang

def angle_oriente_negatif(a, b, c):
    ba = a[..., :2] - b[..., :2]
    bc = c[..., :2] - b[..., :2]

    zero_mask = (np.linalg.norm(ba, axis=-1) < 1e-15) | \
                (np.linalg.norm(bc, axis=-1) < 1e-15)

    dot = np.sum(ba * bc, axis=-1)
    det = ba[..., 0] * bc[..., 1] - ba[..., 1] * bc[..., 0]

    ang = np.arctan2(-det, dot)

    ang = np.where(np.abs(ang) < 1e-10, 0.0, ang)

    ang = np.where(ang < 0, ang + 2*np.pi, ang)

    ang = np.where(zero_mask, 0.0, ang)

    return ang

def point_proche_route(route, P, seuil):
    """
    Retourne True/tableau de bool selon que la distance de P à la polyligne
    route est inférieure à seuil.

    route : (K, ≥2) — séquence de waypoints
    P     : (≥2,)       → retourne un bool scalaire
            (M, ≥2)     → retourne un tableau bool (M,)
    seuil : float
    """
    route  = np.asarray(route, dtype=float)
    P      = np.asarray(P,     dtype=float)
    scalar = P.ndim == 1
    if scalar:
        P = P[np.newaxis, :]       # (1, ≥2)

    P   = P[:, :2]                 # (M, 2)
    A   = route[:-1, :2]           # (K-1, 2)
    B   = route[1:,  :2]
    AB  = B - A                    # (K-1, 2)
    ab2 = np.sum(AB ** 2, axis=1)  # (K-1,)

    # AP : (M, K-1, 2)  via broadcasting
    AP  = P[:, np.newaxis, :] - A[np.newaxis, :, :]
    t   = np.where(ab2 > 1e-15,
                   np.sum(AP * AB, axis=2) / ab2,
                   0.0)                           # (M, K-1)
    t   = np.clip(t, 0.0, 1.0)

    closest = A + t[:, :, np.newaxis] * AB        # (M, K-1, 2)
    d2      = np.sum((P[:, np.newaxis, :] - closest) ** 2, axis=2)  # (M, K-1)
    result  = np.sqrt(d2.min(axis=1)) < seuil     # (M,)

    return bool(result[0]) if scalar else result