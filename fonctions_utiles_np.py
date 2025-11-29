import numpy as np
import math


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

def angleoriente3_np_vect(a, b, c):
    # Vecteurs BA et BC
    ba = a[..., :2] - b[..., :2]
    bc = c[..., :2] - b[..., :2]

    # Normes nulles → angle = 0
    zero_mask = (np.linalg.norm(ba, axis=-1) < 1e-15) | \
                (np.linalg.norm(bc, axis=-1) < 1e-15)

    # Produit scalaire et déterminant 2D
    dot = np.sum(ba * bc, axis=-1)
    det = ba[..., 0] * bc[..., 1] - ba[..., 1] * bc[..., 0]

    # Angle orienté
    ang = np.arctan2(det, dot)
    
	# Très petits angles → 0
    ang = np.where(np.abs(ang) < 1e-10, 0.0, ang)

    # Normalisation dans [0, 2π[
    ang = np.where(ang < 0, ang + 2*np.pi, ang)

    # Appliquer 0 où les vecteurs sont nuls
    ang = np.where(zero_mask, 0.0, ang)

    return ang

