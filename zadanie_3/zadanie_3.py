from pyproj import Geod
import numpy as np

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

def vincenty(BA,LA,BB,LB):
    '''
    Parameters
    ----------
    BA : szerokosc geodezyjna punktu A [RADIAN]
    LA : dlugosc geodezyjna punktu A [RADIAN]
    BB : szerokosc geodezyjna punktu B [RADIAN]
    LB : dlugosc geodezyjna punktu B [RADIAN]

    Returns
    -------
    sAB : dlugosc linii geodezyjnej AB [METR]
    A_AB : azymut linii geodezyjnej AB [RADIAN]
    A_BA : azymut odwrotny linii geodezyjne [RADIAN]
    '''
    b = a * np.sqrt(1-e2)
    f = 1-b/a
    dL = LB - LA
    UA = np.arctan((1-f)*np.tan(BA))
    UB = np.arctan((1-f)*np.tan(BB))
    L = dL
    while True:
        sin_sig = np.sqrt((np.cos(UB)*np.sin(L))**2 +\
                          (np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L))**2)
        cos_sig = np.sin(UA)*np.sin(UB) + np.cos(UA) * np.cos(UB) * np.cos(L)
        sig = np.arctan2(sin_sig,cos_sig)
        sin_al = (np.cos(UA)*np.cos(UB)*np.sin(L))/sin_sig
        cos2_al = 1 - sin_al**2
        cos2_sigm = cos_sig - (2 * np.sin(UA) * np.sin(UB))/cos2_al
        C = (f/16) * cos2_al * (4 + f*(4 - 3 * cos2_al))
        Lst = L
        L = dL + (1-C)*f*sin_al*(sig+C*sin_sig*(cos2_sigm+C*cos_sig*(-1 + 2*cos2_sigm**2)))
        if abs(L-Lst)<(0.000001/206265):
            break
    
    u2 = (a**2 - b**2)/(b**2) * cos2_al
    A = 1 + (u2/16384) * (4096 + u2*(-768 + u2 * (320 - 175 * u2)))
    B = u2/1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    d_sig = B*sin_sig * (cos2_sigm + 1/4*B*(cos_sig*(-1+2*cos2_sigm**2)\
            - 1/6 *B*cos2_sigm * (-3 + 4*sin_sig**2)*(-3+4*cos2_sigm**2)))
    sAB = b*A*(sig-d_sig)
    A_AB = np.arctan2((np.cos(UB) * np.sin(L)),(np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L)))
    A_BA = np.arctan2((np.cos(UA) * np.sin(L)),(-np.sin(UA)*np.cos(UB) + np.cos(UA)*np.sin(UB)*np.cos(L))) + np.pi
    return sAB, A_AB, A_BA

def kivioj(phi, lam, s, az):
    n = round(s / 1000)
    ds = s/n

    for i in range(n):
        N_i = Np(phi)
        M_i = Mp(phi)

        dphi_i = ds*np.cos(az)/M_i
        dA_i = ds*np.sin(az)*np.tan(phi)/N_i

        phi_im = phi + dphi_i/2
        az_im = az + dA_i/2

        N_i = Np(phi_im)
        M_i = Mp(phi_im)

        dphi_i = ds*np.cos(az_im)/M_i
        dA_i = ds*np.sin(az_im)*np.tan(phi_im)/N_i
        dlam_i = ds*np.sin(az_im)/N_i/np.cos(phi_im)

        phi = phi + dphi_i
        lam = lam + dlam_i
        az = az + dA_i

    az_odw = (az + np.pi) % (2*np.pi)
    return phi, lam, az_odw

nr = 4

point_1_phi_rad = np.deg2rad(50 + nr * 15/60) 
point_1_lam_rad = np.deg2rad(18 + nr * 15/60)

s = [40000, 100000, 40000, 100000]
az = np.deg2rad([0, 90, 180, 270])

current_phi = point_1_phi_rad
current_lam = point_1_lam_rad

new_points_rad = []

for p in range(len(s)):
    current_phi, current_lam, az_odw = kivioj(current_phi, current_lam, s[p], az[p])
    new_points_rad.append([current_phi, current_lam, az_odw])

new_points_rad = np.array(new_points_rad)
new_points_deg = np.rad2deg(new_points_rad)

def angle_formatter(degrees, minutes, seconds):
    return f'{degrees}°{minutes:02d}\'{seconds:0.5f}\"'

point_2_phi_d, point_2_phi_m, point_2_phi_s = rad2dms(new_points_rad[0][0])
point_2_lam_d, point_2_lam_m, point_2_lam_s = rad2dms(new_points_rad[0][1])

print(f'Współrzędne punktu 2 wynoszą {angle_formatter(point_2_phi_d, point_2_phi_m, point_2_phi_s)} {angle_formatter(point_2_lam_d, point_2_lam_m, point_2_lam_s)}.')

point_3_phi_d, point_3_phi_m, point_3_phi_s = rad2dms(new_points_rad[1][0])
point_3_lam_d, point_3_lam_m, point_3_lam_s = rad2dms(new_points_rad[1][1])

print(f'Współrzędne punktu 3 wynoszą {angle_formatter(point_3_phi_d, point_3_phi_m, point_3_phi_s)} {angle_formatter(point_3_lam_d, point_3_lam_m, point_3_lam_s)}.')

point_4_phi_d, point_4_phi_m, point_4_phi_s = rad2dms(new_points_rad[2][0])
point_4_lam_d, point_4_lam_m, point_4_lam_s = rad2dms(new_points_rad[2][1])

print(f'Współrzędne punktu 4 wynoszą {angle_formatter(point_4_phi_d, point_4_phi_m, point_4_phi_s)} {angle_formatter(point_4_lam_d, point_4_lam_m, point_4_lam_s)}.')

point_1g_phi_d, point_1g_phi_m, point_1g_phi_s = rad2dms(new_points_rad[3][0])
point_1g_lam_d, point_1g_lam_m, point_1g_lam_s = rad2dms(new_points_rad[3][1])

print(f'Współrzędne punktu 1* wynoszą {angle_formatter(point_1g_phi_d, point_1g_phi_m, point_1g_phi_s)} {angle_formatter(point_1g_lam_d, point_1g_lam_m, point_1g_lam_s)}.')

s_1_1g, az_1_1g, az_1g_1 = vincenty(point_1_phi_rad, point_1_lam_rad, new_points_rad[3][0], new_points_rad[3][1])

print(f'Odległość pomiędzy punktami 1 i 1* wynosi {s_1_1g:.3f}m.')

s_4_1, az_4_1, az_1_4 = vincenty(new_points_rad[2][0], new_points_rad[2][1], point_1_phi_rad, point_1_lam_rad)
az_4_1 = (az_4_1) % (2*np.pi)
d, m, s = rad2dms(az_4_1)

print(f'Odległość z punktu 4 do punktu 1 wynosi {s_4_1:0.3f} m, czyli {s_4_1/1000:0.6f} km.')
print(f'Azymut z punktu 4 do punktu 1 wynosi {d}°{m:02d}\'{s:0.5f}\".')

#-------------------WIZUALIZACJA-------------------
import folium

m = folium.Map(location=[51, 19], zoom_start=9)

points = []

for i in range(len(new_points_deg)):
    point_label = (i + 1) % len(new_points_deg) + 1
    points.append([new_points_deg[i][0], new_points_deg[i][1]])
    folium.Marker(points[i], tooltip=f'Punkt {point_label}').add_to(m)

points.append([np.rad2deg(point_1_phi_rad), np.rad2deg(point_1_lam_rad)])
folium.Marker(points[-1], tooltip=f'Punkt 1*').add_to(m)

folium.Polygon(points, color='red', weight=2.5).add_to(m)

m.save('mapa.html')

#-------------------OBLICZANIE POLA-------------------
from pyproj import Geod
from shapely.geometry import Polygon

geod = Geod(ellps='GRS80')
poly = Polygon(points)
area = geod.geometry_area_perimeter(poly)[0]

print(f'Pole powierzchni wynosi {area:.6f} m^2, czyli {area/1000000:.12f} km^2.')



