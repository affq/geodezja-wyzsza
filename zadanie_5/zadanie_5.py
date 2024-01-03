import numpy as np
from funkcje import gpx2df, read_route_df, flh2xyz, zmniejszenieSiatki, interpolacja
import plotly.graph_objects as go
import plotly.io as pio

szlak = 'obwodnica-rowerowa-po-narwianskim-parku-narodowym'

plik = f'{szlak}.gpx'
gpx2df(plik, f'{szlak}.csv')
df = read_route_df(f'{szlak}.csv')

pio.renderers.default = 'browser'

fig = go.Figure(data=go.Scattermapbox(lat = df['latitude'],
    lon=df['longitude'], mode='lines'))
fig.update_layout(mapbox_style='open-street-map')
fig.update_layout(mapbox_center_lon=22.859, mapbox_center_lat=53.073)
fig.update_layout(mapbox_zoom=9)
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
        
        xyz1 = flh2xyz(np.deg2rad(phi1), np.deg2rad(lam1), 200)
        xyz2 = flh2xyz(np.deg2rad(phi2), np.deg2rad(lam2), 200)
        dxyz = np.array(xyz2) - np.array(xyz1)
        
        d = np.linalg.norm(dxyz)
        dists.append(d)

dists.append(np.nan)

df['dists'] = dists
df['distance'] = df['dists'].cumsum()

model = np.genfromtxt('2011-KRON86-NH.txt')
model2 = zmniejszenieSiatki(model, phis, lams)

dzety = []

for phi, lam in zip(phis, lams):
    dzeta = interpolacja(model2, phi, lam)
    dzety.append(dzeta)

df['hel'] = df['elevation'] + np.array(dzety)

###################################33
dzety = []

model = np.genfromtxt('2021-EVRF2007-NH.txt')
model2 = zmniejszenieSiatki(model, phis, lams)

for phi, lam in zip(phis, lams):
    dzeta = interpolacja(model2, phi, lam)
    dzety.append(dzeta)

df['HN'] = df['hel'] - np.array(dzety)


# obliczenie różnic wysokości normalnych między kronsztadem a amsterdamem
df['HN_diff'] = df['HN'] - df['elevation']

fig = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['HN_diff'], mode='lines'))
fig.update_layout(title='Różnica wysokości normalnych między kronsztadem (2011-KRON86-NH) a amsterdamem(2021-EVRF2007-NH)',
                  xaxis_title='Dystans [km]',
                  yaxis_title='Różnica wysokości normalnych [m]')
fig.show()

#profil trasy - wysokość elipsoidalna
fig = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['hel'], mode='lines'))
fig.update_layout(title='Profil trasy - wysokość elipsoidalna',
                  xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość elipsoidalna [m]')
fig.show()

#profil trasy - 2021-EVRF2007-NH
fig = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['HN'], mode='lines'))
fig.update_layout(title='Profil trasy - 2021-EVRF2007-NH',
                  xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość normalna [m]')
fig.show()

# profil trasy - 2011-KRON86-NH
fig = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['elevation'], mode='lines'))
fig.update_layout(title='Profil trasy - 2011-KRON86-NH',
                  xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość normalna [m]')
fig.show()


