from pyproj import Geod
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vincenty import vincenty

s = 40000
az = np.deg2rad(45)

a = 6378137 # wielka półoś elipsoidy GRS80 w metrach
e2 = 0.00669438002290

def Np(phi):
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    return N

def Mp(phi):
    M = a*(1-e2)/((1-e2*np.sin(phi)**2)**(3/2))
    return M

def rad2dms(rad):
    dd = np.rad2deg(rad)
    dd = dd
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def kivioj(phi, lam, s, az):
    n = round(s / 1000)
    ds = s/n

    for i in range(n):
        #2. obliczenie M i N dla phi_i
        N_i = Np(phi)
        M_i = Mp(phi)

        #3. przyrost do szerokości i długości geograficznej
        dphi_i = ds*np.cos(az)/M_i
        dA_i = ds*np.sin(az)*np.tan(phi)/N_i

        # 4.
        phi_im = phi + dphi_i/2
        az_im = az + dA_i/2

        #5.
        N_i = Np(phi_im)
        M_i = Mp(phi_im)

        #6.
        dphi_i = ds*np.cos(az_im)/M_i
        dA_i = ds*np.sin(az_im)*np.tan(phi_im)/N_i
        dlam_i = ds*np.sin(az_im)/N_i/np.cos(phi_im)

        #7.
        phi = phi + dphi_i
        lam = lam + dlam_i
        az = az + dA_i

    az_odw = az + np.pi
    return phi, lam, az_odw

# g = Geod(ellps='WGS84')

# phi_k, lam_k, az_odw_k = kivioj(phi, lam, s, az)

# lam_p, phi_p, az_odw_p = g.fwd(lam, phi, az, s, radians=True)

# print(f'phi Kivioj {rad2dms(phi_k)}, phi Pyproj {rad2dms(phi_p)}')
nr = 4

punkt_1_phi = np.deg2rad(50 + nr * 15/60) 
punkt_1_lam = np.deg2rad(18 + nr * 15/60)
s_1 = 40000
az_1 = np.deg2rad(0)

punkt_1_phi_k, punkt_1_lam_k, az_1_odw_k = kivioj(punkt_1_phi, punkt_1_lam, s_1, az_1)
print (f'punkt_1_phi_k = {rad2dms(punkt_1_phi_k)}, punkt_1_lam_k = {rad2dms(punkt_1_lam_k)}, az_1_odw_k = {rad2dms(az_1_odw_k)}')
