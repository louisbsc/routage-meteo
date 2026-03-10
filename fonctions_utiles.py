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



