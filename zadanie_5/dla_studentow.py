# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 16:32:41 2023

@author: mgrzy
"""

import numpy as np
from funkcje import gpx2df, read_route_df, save_df

import gpxpy
import plotly.graph_objects as go
import plotly.io as pio
'''
zaimportowanie funkcji, która przelicza współrzędne phi, lambda, h do XYZ - w moim przypadku funkcja blh2xyz
'''
# plik = 'szlak-rowerowy-hajnowka-topilo-hajnowka.gpx'

# gpx2df(plik, 'szlak-rowerowy-hajnowka-topilo-hajnowka.csv')

df = read_route_df('szlak-rowerowy-hajnowka-topilo-hajnowka.csv')

pio.renderers.default = 'browser'

fig = go.Figure(data=go.Scattermapbox(lat = df['latitude'],
    lon=df['longitude'], mode='lines'))
fig.update_layout(mapbox_style='open-street-map')
fig.show()

phis = df['latitude']
lams = df['longitude']

dists = []
for i, _ in enumerate(phis):
    if i < len(phis)-1:
        phi1 = phis[i]
        lam1 = lams[i]
        phi2 = phis[i+1]
        lam2 = lams[i+1]
        
        xyz1 = blh2xyz(np.deg2rad(phi1), np.deg2rad(lam1), 200)
        xyz2 = blh2xyz(np.deg2rad(phi2), np.deg2rad(lam2), 200)
        dxyz = np.array(xyz2) - np.array(xyz1)
        
        d = np.linalg.norm(dxyz)
        dists.append(d)
dists.append(np.nan)

df['dists'] = dists
df['distance'] = df['dists'].cumsum()

model = np.genfromtxt('siatka_gugik-geoid2011-PL-KRON86-NH.txt')

def zmniejszenieSiatki(model, phis, lams, grid_step=0.01):
    max_phi = np.max(phis)
    min_phi = np.min(phis)
    max_lam = np.max(lams)
    min_lam = np.min(lams)
    
    ind_phi = np.logical_and(model[:,0]<(max_phi+grid_step), model[:,0]>(min_phi-grid_step))
    ind_lam = np.logical_and(model[:,1]<(max_lam+grid_step), model[:,1]>(min_lam-grid_step))
    
    indeksy = np.logical_and(ind_phi, ind_lam)
    model2 = model[indeksy, :]
    return model2
model2 = zmniejszenieSiatki(model, phis, lams)
#%%


def interpolacja(model2, phi, lam, grid_step=0.01):
    ind_phi = np.logical_and(model2[:,0]<(phi+grid_step), 
                             model2[:,0]>(phi-grid_step))
    ind_lam = np.logical_and(model2[:,1]<(lam+grid_step), 
                             model2[:,1]>(lam-grid_step))
    
    indeksy = np.logical_and(ind_phi, ind_lam)
    model3 = model2[indeksy, :]

    (y1,x1,Q11),(_y1,x2,Q21),(y2,_x1,Q12),(_y2,_x2,Q22) = model3
    x = lam
    y = phi
    R1 = Q11 + (Q21-Q11)/(x2-x1) * (x-x1)
    R2 = Q12 + (Q22-Q12)/(x2-x1) * (x-x1)
    
    P = R1 + (R2-R1)/(y2-y1) * (y-y1)
    return P


dzety = []
for phi, lam in zip(phis, lams):
    dzeta = interpolacja(model2, phi, lam)
    dzety.append(dzeta)
    print(phi, lam)

df['hel'] = df['elevation'] + np.array(dzety)








