import datetime
import math
import numpy as np
import matplotlib.pyplot as plt

def julday(y,m,d,h):
    if m <= 2:
        y = y - 1
        m = m + 12
    jd = np.floor(365.25*(y+4716))+np.floor(30.6001*(m+1))+d+h/24-1537.5;
    return jd

def GMST(jd):
    T = (jd - 2451545) / 36525
    Tu = jd - 2451545
    g = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2-T**3/38710000
    g = (g%360) / 15
    return g

def dms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg)
    return rad

def hms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg*15)
    return rad

def hms2h(dms):
    return dms[0] + dms[1]/60 + dms[2]/3600

def calculate_h(delta_r, fi_r, t_r):
    h = math.asin(math.sin(delta_r)*math.sin(fi_r) + math.cos(delta_r)*math.cos(fi_r)*math.cos(t_r))
    return h

def calculate_a(delta_r, fi_r, t_r):
    a = -math.cos(delta_r)*math.sin(t_r)
    b = math.sin(delta_r)*math.cos(fi_r) - math.cos(delta_r)*math.sin(fi_r)*math.cos(t_r)
    az = math.atan2(a, b)
    return az

def add_point(a, h, r, ax):
    x = r * np.sin(a) * np.cos(h)
    y = r * np.cos(a) * np.cos(h)
    z = r * np.sin(h)

    ax.scatter(x, y, z, marker='*', color="orange")
    ax.text(x, y, z, hr)

# right ascension
alfa = [4, 37, 16.316]
alfa_h = hms2h(alfa)
alfa_r = hms2rad(alfa)

# declination
delta = [16, 33, 16.570]
delta_r = dms2rad(delta)

# warsaw coordinates
warsaw = [52, 21]
warsaw_r = math.radians(warsaw[0])

# equator coordinates
equator = [0, 21]
equator_r = math.radians(equator[0])

#WYKRES 3D
# fig = plt.figure(figsize = (10,10))
# ax = fig.add_subplot(projection = '3d')
# r = 1
# u, v = np.mgrid[0:(2 * np.pi+0.1):0.1, 0:np.pi:0.1]
# x = np.cos(u) * np.sin(v)
# y = np.sin(u) * np.sin(v)
# z = np.cos(v)
# z[z<0] = 0
# ax.plot_surface(x,y,z, alpha = 0.1)
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('wykres 3D sfery niebieskiej - FK5 168', pad = 30, fontsize = 15)
# ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.8)

# WYKRES SKYPLOT
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(polar = True)
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1)
# ax.set_yticks(range(0, 90+10, 10))
# yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
# ax.set_yticklabels(yLabel)
# ax.set_rlim(0,90)
# ax.set_title('wykres skyplot - FK5 168', pad = 30, fontsize = 15)

#WYKRES ZALEZNOSCI WYSOKOSCI OD CZASU
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot()
# ax.set_xlabel('godzina')
# ax.set_ylabel('wysokość [deg]')
# ax.set_ylim(-90,90)
# ax.set_xlim(0,24)
# ax.axhline(0, linestyle = '-.', color = 'black')
# ax.set_title('wykres wysokości gwiazdy FK5 168 w zależności od czasu- 1 lipca 2023', pad = 30, fontsize = 15)
# ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.5)


# WYKRES PANORAMA
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot()
# ax.set_xlabel('azymut [deg]')
# ax.set_ylabel('wysokość [deg]')
# ax.set_ylim(0,90)
# ax.set_xlim(-180,180)
# ax.set_title('wykres panorama FK5 168', pad = 30, fontsize = 15)
# ax.grid(color = 'grey', linestyle = '-', linewidth = 0.25, alpha = 0.5)

for hr in range(2,26):
    if hr >= 24:
        hr = hr - 24
        gd = datetime.datetime(2023, 7, 2, hr)
    else:
        gd = datetime.datetime(2023, 7, 1, hr)
    # julian date
    jd = julday(gd.year, gd.month, gd.day, gd.hour)
    # greenwich mean sidereal time
    gmst = GMST(jd)

    # local sidereal time
    lst_warsaw = gmst + warsaw[1]/15

    # hour angle
    t_warsaw = lst_warsaw - alfa_h
    t_warsaw_r = math.radians(t_warsaw * 15) 
    t_equator_r = t_warsaw_r
    

    # calculate h and a for warsaw
    h_warsaw_r = calculate_h(delta_r, warsaw_r, t_warsaw_r)
    a_warsaw_r = calculate_a(delta_r, warsaw_r, t_warsaw_r)

    # calculate h and a for equator
    h_equator_r = calculate_h(delta_r, equator_r, t_equator_r)
    a_equator_r = calculate_a(delta_r, equator_r, t_equator_r)

# adding points to SPHERE
    # add_point(a_equator_r, h_equator_r, r, ax)
    # add_point(a_warsaw_r, h_warsaw_r, r, ax)

# adding points to SKYPLOT
    # if (np.rad2deg(h_warsaw_r) >= 0):
    #     ax.scatter(a_warsaw_r, 90-np.rad2deg(h_warsaw_r), marker = "*", color = 'orange', s = 50)
    #     ax.text(a_warsaw_r, 90-np.rad2deg(h_warsaw_r), hr, color = 'black')
    #     
    # if (np.rad2deg(h_equator_r) >= 0):
    #     ax.scatter(a_equator_r, 90-np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    #     ax.text(a_equator_r, 90-np.rad2deg(h_equator_r), hr)

# adding points to WYKRES ZALEZNOSCI WYSOKOSCI OD CZASU
    # ax.scatter(hr, np.rad2deg(h_warsaw_r), marker = "*", color = 'orange')
    # ax.text(hr, np.rad2deg(h_warsaw_r), hr)

    # ax.scatter(hr, np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    # ax.text(hr, np.rad2deg(h_equator_r), hr)

# adding points to PANORAMA
    # if (np.rad2deg(h_warsaw_r) >= 0):
    #     ax.scatter(np.rad2deg(a_warsaw_r), np.rad2deg(h_warsaw_r), marker = "*", color = 'orange')
    #     ax.text(np.rad2deg(a_warsaw_r), np.rad2deg(h_warsaw_r), hr)

    # if (np.rad2deg(h_equator_r) >= 0):
    #     ax.scatter(np.rad2deg(a_equator_r), np.rad2deg(h_equator_r), marker = '*', color = 'orange')
    #     ax.text(np.rad2deg(a_equator_r), np.rad2deg(h_equator_r), hr)

# plt.show()

