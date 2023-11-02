# from geofunc import hirvonen, xyz2neu, neu2has
from geofunc import Rneu, flh2xyz
import numpy as np
from read_flightradar import read_flightradar
import folium
import matplotlib.pyplot as plt

map = folium.Map(location=[52.164318, 20.981787], zoom_start=6)

file = 'lot4.csv'
dane = read_flightradar(file)

flh = dane[:,[7,8,9]]
print (flh)
flh[:,-1] = flh[:,-1]*0.3048 + 135.40 # (wysokość elipsoidalna lotniska) przeliczenie na metry nad poziomem lotniska

flh_lotnisko = flh[0,:]
xyz_lotnisko = flh2xyz(np.deg2rad(flh_lotnisko[0]), np.deg2rad(flh_lotnisko[1]), flh_lotnisko[2])

flh = flh[66:,:] #- w zależności od wersji usunąć kiedy samolot stoi

R = Rneu(np.deg2rad(flh_lotnisko[0]), np.deg2rad(flh_lotnisko[1]))
azymuty = []


for fi, lam, h in flh:
    xyz = flh2xyz(fi, lam, h)

    xsl = np.array(xyz) - xyz_lotnisko
    neu = R.T @ xsl # 3-elementowy wektor [N, E, U]

    az = np.rad2deg(np.arctan2(neu[1], neu[0]))
    azymuty.append(az)
    print(az)

    folium.Marker([fi, lam], popup='Azymut: ' + str(az)).add_to(map)

map.save('mapa.html')

