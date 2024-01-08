import numpy as np

dane = np.genfromtxt('0NOA.txyz2.txt')

t = dane[:, 2]
xyz = dane[:, 3:6]

dxyz = xyz - xyz[0, :]

import matplotlib.pyplot as plt
fig, ax = plt.subplots(3, 1, sharex=True)
for i in range(len(ax)):
    ax[i].plot(t, dxyz[:, i])

plt.show()

jedynki = np.ones(len(t))
A = np.column_stack((t, jedynki))
y_linia_lista = []

for i in range(3):
    y = dxyz[:, 0]
    xx = np.linalg.inv(A.T @ A) @ (A.T @ y)

    y_linia = A @ xx
    y_linia_lista.append(y_linia)
y_linia_lista = np.array(y_linia_lista)

fig, ax = plt.subplots()

for i in range(len(ax)):
    ax[i].plot(t, dxyz[:, i])

plt.show()

######################
from hirvonen import hirv, Rneu

fi, lam, h = hirv(*xyz[0, :])
R = Rneu(fi, lam)

dneu = []

for dx in dxyz:
    dneu.append(R.T @ dx)

dneu = np.array(dneu)

y_linia_lista = []

for i in range(3):
    y = dneu[:, 0]
    xx = np.linalg.inv(A.T @ A) @ (A.T @ y)

    y_linia = A @ xx
    y_linia_lista.append(y_linia)
y_linia_lista = np.array(y_linia_lista)

fig, ax = plt.subplots(3, 1, sharex=True)

for i in range(len(ax)):
    ax[i].plot(t, dxyz[:, i])
    ax[i].plot(t, y_linia_lista[:,i])

plt.show()

import geopandas as gpd
shpfile = gpd.read_file('CNTR_BN_03M_2020_4326.shp')

fig, ax = plt.subplots()
shpfile.plot(ax=ax)

vn = (y_linia_lista[-1, 0] - y_linia_lista[0, 0]) / (t[-1] - t[0])*100
ve = (y_linia_lista[-1, 1] - y_linia_lista[0, 1]) / (t[-1] - t[0])*100

ax.quiver(np.rad2deg(lam), np.rad2deg(fi), ve, vn)

plt.show()