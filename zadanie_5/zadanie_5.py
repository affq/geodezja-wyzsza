from funkcje import gpx2df, read_route_df, save_df
import gpxpy.gpx
import pandas as pd
import scipy
import plotly.graph_objects as go   
import plotly.io as pio
import numpy as np

dists = []

plik = 'szlak-orla-bialego.gpx'

gpx2df(plik, new_file_name = 'route_df.csv')
df = read_route_df('route_df.csv')
pio.renderers.default = "browser"

fig = go.Figure(data = go.Scattermapbox(lat=df['latitude'], lon=df['longitude'], mode='lines', marker=go.scattermapbox.Marker(size=9)))
fig.update_layout(mapbox_style="open-street-map")
fig.show()

phis = df['latitude']
lams = df['longitude']

def flh2xyz(phi, lamb, h, a=6378137, e2=0.0066943800229007):
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    x = (N + h)*np.cos(phi)*np.cos(lamb)
    y = (N + h)*np.cos(phi)*np.sin(lamb)
    z = (N*(1-e2)+h)*np.sin(phi)
    return [x, y, z] 

for i, _ in enumerate(phis):
    if i < len(phis) - 1:
        phi1 = phis[i]  
        lam1 = lams[i]
        phi2 = phis[i+1]
        lam2 = lams[i+1]

        xyz1 = flh2xyz(np.deg2rad(phi1), np.deg2rad(lam1), 0)
        xyz2 = flh2xyz(np.deg2rad(phi2), np.deg2rad(lam2), 0)
        dist = np.linalg.norm((np.array(xyz2)-(xyz1)))
        dists.append(dist)
dists.append(np.nan)
df['dists'] = dists
df['distance'] = df['dists'].cumsum()

model = np.genfromtxt('siatka_gugik-geoid2011-PL-KRON86-NH.txt')

max_phi = np.max(phis)
max_lam = np.max(lams)
min_phi = np.min(phis)
min_lam = np.min(lams)

ind_phi = np.logical_and(model[:,0]<(max_phi+0.1), model[:,0]>(min_phi-0.1))
ind_lam = np.logical_and(model[:,1]<(max_lam+0.1), model[:,1]>(min_lam-0.1))

indeksy = np.logical_and(ind_phi, ind_lam)
model2 = model[indeksy,:]

for phi, lam in zip(phis, lams):
    print (phi, lam)

phi = phis[0]
lam = lams[0]

phi_u = np.unique(phis)
lam_u = np.unique(lams)

phi_grid = np.reshape(model2[:,0], (np.len(phi_u), -1))
lam_grid = np.reshape(model2[:,1], (np.len(phi_u), -1))
zeta_grid = np.reshape(model2[:,2], (np.len(phi_u), -1))

phi_0_ind = np.argwhere(phi_grid[:,0]<phi)[-1][0]
phi_1_ind = np.argwhere(phi_grid[:,0]>phi)[-1][0]
phi_0 = phi_grid[phi_0_ind,0]
phi_1 = phi_grid[phi_1_ind,0]

lam_0_ind = np.argwhere(lam_grid[0,:]<lam)[-1][0]
lam_1_ind = np.argwhere(lam_grid[0,:]>lam)[0][-1]
lam_0 = lam_grid[0,lam_0_ind]
lam_1 = lam_grid[0,lam_1_ind]

zeta_siatka = zeta_grid[[phi_0_ind, phi_1_ind], [lam_0_ind, lam_1_ind]]
phi_siatka = phi_grid[[phi_0_ind, phi_1_ind], [lam_0_ind, lam_1_ind]]
lam_siatka = lam_grid[[phi_0_ind, phi_1_ind], [lam_0_ind, lam_1_ind]]
