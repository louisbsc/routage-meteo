import matplotlib.pyplot as plt
from math import *

def V1(p, t):
	return 0, 5

def vent(V, t, intx, inty, n):
	x = (n + 1) * [intx[0] + k * (intx[1] - intx[0]) / n for k in range(n + 1)]
	y = []
	for k in range(n + 1):
		y = y + (n + 1) * [inty[0] + k * (inty[1] - inty[0]) / n]
	u = [- V([x[k], y[k]], t)[1] * sin(radians(V([x[k], y[k]], t)[0])) for k in range((n + 1)**2)]
	v = [- V([x[k], y[k]], t)[1] * cos(radians(V([x[k], y[k]], t)[0])) for k in range((n + 1)**2)]
	c = []
	for k in range((n + 1)**2):
		coul = V([x[k], y[k]], t)[1] / 30
		c.append((1 - coul, 1 - coul, 1 - coul))
	plt.quiver(x, y, u, v, color = c)
	plt.axis('equal')
	

def vent_circulaire_2(p, t):
    x, y = p[0], p[1]

    r = hypot(x, y)
    if r == 0:
        return 0.0, 0.0  # Direction indéfinie au centre

    # Force du vent : 10 kts proche du centre, asymptote à 20 kts
    # r exprimé en NM (supposé)
    force = 20 * (1 - exp(-r / 40))

    # Vent tangent enroulement horaire : (y, -x)
    vx = force * (y / r)
    vy = -force * (x / r)

    # Angle mathématique du vecteur (0 rad = +x, CCW)
    angle_rad = atan2(vy, vx)

    # Conversion en angle "rose des vents" (0° = Nord venant du haut, sens horaire)
    direction_deg = (90 - degrees(angle_rad)) % 360

    return direction_deg, force