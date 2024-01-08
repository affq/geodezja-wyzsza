import numpy as np

e2 = 0.00669438002290

def Np(phi):
    a = 6378137.0
    e2 = 0.00669438002290
    N =  a / (np.sqrt(1 - e2 * np.sin(np.deg2rad(phi))**2))
    return N

def hirv(x,y,z):
    p = np.sqrt(x**2+y**2)
    
    phi = np.arctan(z/(p*(1-e2)))
    
    while True:
        N = Np(phi)
        h = p/np.cos(phi) - N
        phi_poprzednie = phi
        phi = np.arctan((z/p) * (1-(N*e2)/(N+h))**(-1))
        if abs(phi_poprzednie-phi)<(0.000001/60/60/360):
            break
    lam = np.arctan2(y,x)
    return phi, lam, h

def Rneu(phi, lamb):
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R