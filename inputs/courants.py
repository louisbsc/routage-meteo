import numpy as np

def C0(p, t):
    return 0, 0

def C1(p, t):
    x, y = p[:2]
    if x < 0:
        return 180, 1
    else:
        return 180, 0
    
def C2(p, t):
    if t < 10:
        return 180, 1
    else:
        return 180, 0
