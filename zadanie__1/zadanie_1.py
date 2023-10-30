# Zadanie 1. Transformacja współrzędnych gwiazdy z układu
# równikowego do horyzontalnego
# Treść zadania
# Wyznaczyć położenie danej gwiazdy, w układzie współrzędnych lokalnych (horyzontalnych),
# dla dwóch miejsc na powierzchni Ziemi. Obliczenia mają zostać wykonane dla całej doby, 1 lipca
# 2023, w godzinnych interwałach. Położenie gwiazdy należy zwizualizować na wykresach, np.
# wykresie typu skyplot, wykres 3D sfery niebieskiej, wykres liniowy zależności wysokości (i/lub
# azymutu) od czasu, lub wykres "panorama".
# 4.1 Dane do zadania
# Danymi do zadania są współrzędne gwiazdy (tab. 1), w układzie równikowym ekwinokcjal-
# nym, na epokę 2023.5 (czyli 1 lipca 2023). Obliczenia należy wykonać, poczynając od danego
# czasu urzędowego dla Polski (UTC+2).
# Zadanie proszę wykonać dla całej doby 1 lipca 2023, dla dwóch położeń obserwatora:
# •okolice Warszawy: φ = 52◦; λ = 21◦
# •równik, dla tej samej długości geograficznej: φ = 0◦; λ = 21◦

# nr gwiazdy z RA FK5: 168
# rektascensja: 4h 37min 16.316s
# deklinacja: 16° 33min 16.570s

# czas urzędowy dla Polski (UTC+2)
# 1 lipca 2023

# okolice Warszawy: φ = 52◦; λ = 21◦
# równik: φ = 0◦; λ = 21◦

# Kolejność wykonywania zadania
# 1. Uwzględnienie strefy czasowej Warszawy i obliczenie czasu UTC (w naszym przypadku,
# ignorujemy różnicę pomiędzy czasem UTC i UT1, nie uwzględniamy odpowiedniej po-
# prawki, dlatego obliczony przez nas czas UTC będzie równy czasowi UT1) – strefę cza-
# sową dla Polski (UTC+2) można uwzględnić na samym początku obliczeń: wówczas wy-
# niki otrzymamy dla całej doby 1 lipca, od godziny 0 do 24 w czasie UTC+2; lub możemy
# uwzględnić strefę czasową na samym końcu obliczeń, podczas etykietowania wartości na
# wykresach: wówczas otrzymamy wyniki dla doby 1 lipca, ale czasu UTC, czyli od 200
# (UTC+2) dnia 1 lipca do 200 (UTC+2) 2 lipca. Oba rozwiązania będą poprawne,
# 2. Obliczamy lokalny czas gwiazdowy (LST). Możemy to zrobić na dwa sposoby:
# 2.1. obliczamy czas średni gwiazdowy Greenwich, na epokę 0 UT1 (GMST0), na podstawie
# daty juliańskiej danej epoki na godzinę 0 (kod do obliczeń daty juliańskiej znajduje
# się w materiałach),
# 2.2. konkretną epokę obliczamy poprzez dodanie do GMST0 interesującej nas godziny,
# pomnożonej przez wartość 1.002737909 350795 (zamiana jednostek z godzin słonecz-
# nych, na gwiazdowe)
# 2.3. obliczamy LST uwzględniając długość geograficzna miejsca obserwacji.
# lub
# 2.1. obliczamy czas średni gwiazdowy Greenwich (GMST), na dowolną epokę danego dnia,
# na podstawie daty juliańskiej obliczonej dla konkretnej godziny,
# 2.2. obliczamy LST uwzględniając długość geograficzna miejsca obserwacji;
# 3. Obliczenie kąta godzinnego, znając lokalny czas gwiazdowy oraz rektascenzję gwiazdy
# (równanie 2);
# 4. Rozwiązanie trójkąta sferycznego i obliczenie wysokości oraz azymutu gwiazdy (wykonanie
# transformacji współrzędnych równikowych do współrzędnych horyzontalnych)
# 5. Wykonanie odpowiednich wizualizacji.

import datetime
import math
import time_transformations as tt #funkcje zamieniajace dzien na julianski, julianski na czas sredni gwiazdowy greenwich
import numpy as np
import matplotlib.pyplot as plt

# nr gwiazdy z RA FK5: 168
# rektascensja: 4h 37min 16.316s
# deklinacja: 16° 33min 16.570s

alfa = [4, 37, 16.316]
delta = [16, 33, 16.570]

# czas urzędowy dla Polski (UTC+2)
# 1 lipca 2023

# okolice Warszawy: φ = 52◦; λ = 21◦
# równik: φ = 0◦; λ = 21◦
warsaw = [52, 21]
equator = [0, 21]

#WYKRES 1
# fig = plt.figure(figsize = (10,10))
# ax = fig.add_subplot(projection = '3d')
# r = 1
# u, v = np.mgrid[0:(2 * np.pi+0.1):0.1, 0:np.pi:0.1]
# x = np.cos(u) * np.sin(v)
# y = np.sin(u) * np.sin(v)
# z = np.cos(v)
# z[z<0] = 0
# ax.plot_surface(x,y,z, alpha = 0.1)
# # label the axes
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# ax.grid(True)
# ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.8)

# WYKRES 2
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(polar = True)
# ax.set_theta_zero_location('N') # ustawienie kierunku północy na górze wykresu
# ax.set_theta_direction(-1)

# ax.set_yticks(range(0, 90+10, 10))                   # Define the yticks

# yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
# ax.set_yticklabels(yLabel)
# ax.set_rlim(0,90)

# # title for the plot
# ax.set_title('wykres skyplot - FK5 168', pad = 30, fontsize = 15)

#WYKRES 3
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot()
ax.set_xlabel('godzina')
ax.set_ylabel('wysokość [deg]')
ax.set_ylim(-90,90)
ax.set_xlim(0,24)

# ustaw przerywaną linię poziomą na wysokości 0
ax.axhline(0, linestyle = '-.', color = 'black')

#add title
ax.set_title('wykres wysokości gwiazdy FK5 168 w zależności od czasu- 1 lipca 2023', pad = 30, fontsize = 15)

#add grid
ax.grid(True)

#make the grid lighter
ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.5)


# WYKRES 3 - wykres "panorama"

# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot()
# ax.set_xlabel('azymut [deg]')
# ax.set_ylabel('wysokość [deg]')
# ax.set_ylim(0,90)
# ax.set_xlim(-180,180)

# ax.set_title('wykres panorama FK5 168', pad = 30, fontsize = 15)

# ax.grid(True)
# ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.5)

#do not show text for hour that are not visible

# zamiana stopni, minut, sekund na radiany
def dms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg)
    return rad

# zamiana godzin, minut, sekund na radiany
def hms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg*15)
    return rad

def hms2h(dms):
    return dms[0] + dms[1]/60 + dms[2]/3600

# right ascension to hours
alfa_h = hms2h(alfa)
# right ascension to radians
alfa_r = hms2rad(alfa)

# declination to hours, minutes, seconds
delta_r = dms2rad(delta)

# warsaw latitude to radians
warsaw_r = math.radians(warsaw[0])

# equator latitude to radians
equator_r = math.radians(equator[0])

# function to calculate height
def calculate_h(delta_r, fi_r, t_r):
    h = math.asin(math.sin(delta_r)*math.sin(fi_r) + math.cos(delta_r)*math.cos(fi_r)*math.cos(t_r))
    return h

# function to calculate azimuth
def calculate_a(delta_r, fi_r, t_r):
    a = -math.cos(delta_r)*math.sin(t_r)
    b = math.sin(delta_r)*math.cos(fi_r) - math.cos(delta_r)*math.sin(fi_r)*math.cos(t_r)
    az = math.atan2(a, b)
    return az

def add_point(a, h, r):
    x = r * np.sin(a) * np.cos(h)
    y = r * np.cos(a) * np.cos(h)
    z = r * np.sin(h)

    ax.scatter(x, y, z, marker='*', color="orange")

    # add label for every hour
    ax.text(x, y, z, hr)

# obliczamy czas średni gwiazdowy Greenwich, na epokę 0 UT1 (GMST0) dla 1 lipca 2023 co 1 godzine
for hr in range(2,26):
    # georgia date
    if hr >= 24:
        hr = hr - 24
        gd = datetime.datetime(2023, 7, 2, hr)
    else:
        gd = datetime.datetime(2023, 7, 1, hr)
    # julian date
    jd = tt.julday(gd.year, gd.month, gd.day, gd.hour)
    # greenwich mean sidereal time
    gmst = tt.GMST(jd)

    # local sidereal time for warsaw
    lst_warsaw = gmst + warsaw[1]/15
    # local sidereal time for equator
    lst_equator = gmst + equator[1]/15

    # warsaw hour angle
    t_warsaw = lst_warsaw - alfa_h
    t_warsaw_r = math.radians(t_warsaw * 15) 

    # equator hour angle
    t_equator = lst_equator - alfa_h
    t_equator_r = math.radians(t_equator * 15)
    

#############################CALCULATING H AND A FOR WARSAW####################################

    # calculate h for warsaw (in degrees)
    h_warsaw_r = calculate_h(delta_r, warsaw_r, t_warsaw_r)
    h_warsaw_d = math.degrees(h_warsaw_r)

    # calculate azimuth for warsaw (in degrees)
    a_warsaw_r = calculate_a(delta_r, warsaw_r, t_warsaw_r)
    a_warsaw_d = math.degrees(a_warsaw_r)

#############################CALCULATING H AND A FOR EQUATOR####################################

    # calculate azimuth for equator (in degrees)
    a_equator_r = calculate_a(delta_r, equator_r, t_equator_r)

    # calculate h for equator (in degrees)
    h_equator_r = calculate_h(delta_r, equator_r, t_equator_r)

#############################PLOTTING####################################
    # print height and azimuth for warsaw
    print("Warsaw:")
    print("Hour: ", hr)
    print("Height: ", h_warsaw_r)
    print("Azimuth: ", a_warsaw_r)


# #############################SPHERE####################################
    # ax.set_title('wykres 3D sfery niebieskiej - FK5 168', pad = 30, fontsize = 15)
    # add_point(a_equator_r, h_equator_r, r)

    # ax.set_title('wykres 3D sfery niebieskiej - FK5 168', pad = 30, fontsize = 15)
    # add_point(a_warsaw_r, h_warsaw_r, r)

#############################SKYPLOT####################################
    # narysowanie punktu na wykresie 
    # make the stars bigger
    # ax.scatter(a_warsaw_r, 90-np.rad2deg(h_warsaw_r), marker = "*", color = 'orange', s = 50)

    # # # #show text only within the plot
    # if (np.rad2deg(h_warsaw_r) >= 0):
    #     ax.text(a_warsaw_r, 90-np.rad2deg(h_warsaw_r), hr, color = 'black')
    #     
    # make poitns look like stars
    # if (np.rad2deg(h_equator_r) >= 0):
    #     ax.scatter(a_equator_r, 90-np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    #     ax.text(a_equator_r, 90-np.rad2deg(h_equator_r), hr)

#############################WYKRES 3####################################
    ax.scatter(hr, np.rad2deg(h_warsaw_r), marker = "*", color = 'orange')
    ax.text(hr, np.rad2deg(h_warsaw_r), hr)

    # # make poitns look like stars
    # ax.scatter(hr, np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    # ax.text(hr, np.rad2deg(h_equator_r), hr)

    #############################WYKRES 4####################################
    # if (np.rad2deg(h_warsaw_r) >= 0):
    #     ax.scatter(np.rad2deg(a_warsaw_r), np.rad2deg(h_warsaw_r), marker = "*", color = 'orange')
    #     ax.text(np.rad2deg(a_warsaw_r), np.rad2deg(h_warsaw_r), hr)

    # if (np.rad2deg(h_equator_r) >= 0):
    #     ax.scatter(np.rad2deg(a_equator_r), np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    #     ax.text(np.rad2deg(a_equator_r), np.rad2deg(h_equator_r), hr)
plt.show()

