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

spd = dane[:,-2]

speeds = []
azimuths = []
altitudes = []
distances = []
locations = []

for i in spd:
    speeds.append(i)
    print (i)

for fi, lam, h in flh:
    xyz = flh2xyz(np.deg2rad(fi), np.deg2rad(lam), h)

    xsl = np.array(xyz) - xyz_lotnisko
    neu = R.T @ xsl # 3-elementowy wektor [N, E, U]

    az = np.rad2deg(np.arctan2(neu[1], neu[0]))
    azimuths.append(az)

    s = np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2) / 1000
    distances.append(s)

    altitudes.append(h)

    locations.append([fi, lam])

folium.PolyLine(locations).add_to(map)

map.save('mapa.html')

#skyplot
plt.polar(np.deg2rad(azimuths), distances, 'ro')

# wykres wysokości samolotu w zależności od czasu
plt.figure()
plt.plot(altitudes)
plt.ylim(0)
plt.xlim(0)

# wykres zmian prędkości lotu samolotu
plt.figure()
plt.plot(speeds)
plt.ylim(0)
plt.xlim(0)

# wykres zmian odległości samolotu od lotniska
plt.figure()
plt.plot(distances)
plt.ylim(0)
plt.xlim(0)

plt.show()