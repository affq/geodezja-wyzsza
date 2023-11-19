from pyproj import Geod
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vincenty import vincenty

a = 6378137 # wielka półoś elipsoidy GRS80 w metrach
e2 = 0.00669438002290 # pierwsza mimośród elipsoidy GRS80

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

    #8. 
    az_odw = (az + np.pi) % (2*np.pi)
    #9.
    return phi, lam, az_odw

# g = Geod(ellps='WGS84')

# phi_k, lam_k, az_odw_k = kivioj(phi, lam, s, az)

# lam_p, phi_p, az_odw_p = g.fwd(lam, phi, az, s, radians=True)

# print(f'phi Kivioj {rad2dms(phi_k)}, phi Pyproj {rad2dms(phi_p)}')
nr = 4
point_1_phi = np.deg2rad(50 + nr * 15/60) 
point_1_lam = np.deg2rad(18 + nr * 15/60)

s = [40000, 100000, 40000, 100000] #odległości w metrach kolejno 1-2, 2-3, 3-4, 4-1*
az = np.deg2rad([0, 90, 180, 270]) #azymuty w radianach kolejno 1-2, 2-3, 3-4, 4-1*

current_phi = point_1_phi
current_lam = point_1_lam

new_points = []

for p in range(len(s)):
    current_phi, current_lam, az_odw = kivioj(current_phi, current_lam, s[p], az[p])
    new_points.append([current_phi, current_lam, az_odw])

new_points = np.array(new_points)
new_points = np.rad2deg(new_points)
print (new_points)

#Czy po obliczeniu kolejnych wierzchołków ‘trapezu’, na podstawie podanych obserwacji,
# zamkniemy otrzymamy figurę zamkniętą? Jaka będzie różnica położenia punktów 1 i 1*?
# Czy spowodowana będzie otrzymana różnica?

s_1_1g, az_1_1g, az_1g_1 = vincenty(point_1_phi, point_1_lam, np.deg2rad(new_points[3][0]), np.deg2rad(new_points[3][1]))

print(f'Odległość pomiędzy punktami 1 i 1* wynosi {s_1_1g/1000} km')

# Wyznacz właściwe obserwacje: odległość oraz azymut z punktu 4 do punktu 1 (zadanie
# odwrotne – algorytm Vincentego lub inny)
s_4_1, az_4_1, az_1_4 = vincenty(np.deg2rad(new_points[2][0]), np.deg2rad(new_points[2][1]), point_1_phi, point_1_lam)

print(f'Odległość z punktu 4 do punktu 1 wynosi {s_4_1/1000} km')
print(f'Azymut z punktu 4 do punktu 1 wynosi {np.rad2deg(az_4_1 % np.pi)} stopni')

# przedstawienie na mapie punktow - polaczenie ich liniami
import folium
m = folium.Map(location=[50, 18], zoom_start=6)
for i in range(len(new_points)):
    folium.Marker([new_points[i][0], new_points[i][1]], popup=f'Punkt {i+1}').add_to(m)
for i in range(len(new_points)):
    folium.PolyLine([[new_points[i][0], new_points[i][1]], [new_points[(i+1)%len(new_points)][0], new_points[(i+1)%len(new_points)][1]]], color='red').add_to(m)

m.save('mapa.html')

geod = Geod(ellps='WGS84')
area = geod.polygon_area_perimeter(new_points[:,0], new_points[:,1])
print(f'Pole powierzchni wynosi {area[0]/1000000} km^2')




