import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import gpxpy

def gpx2df(gpx_file_name, new_file_name = 'route_df.csv'):
    with open(gpx_file_name, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)
    route_info = []
    
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                route_info.append({
                    'latitude': point.latitude,
                    'longitude': point.longitude,
                    'elevation': point.elevation,
                    'time':point.time
                })
    
    route_df = pd.DataFrame(route_info)
    
    save_df(route_df, new_file_name)
    

def save_df(route_df, new_file_name= 'route_df.csv'):
    route_df.to_csv(new_file_name, index=False)
    print(f'file saved to {new_file_name}')
    
def read_route_df(file):
    route_df = pd.read_csv(file)
    return route_df

def flh2xyz(phi, lamb, h, a=6378137, e2=0.0066943800229007):
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    x = (N + h)*np.cos(phi)*np.cos(lamb)
    y = (N + h)*np.cos(phi)*np.sin(lamb)
    z = (N*(1-e2)+h)*np.sin(phi)
    return [x, y, z]

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

app = Dash(__name__)

szlak = 'obwodnica-rowerowa-po-narwianskim-parku-narodowym'

plik = f'{szlak}.gpx'
gpx2df(plik, f'{szlak}.csv')
df = read_route_df(f'{szlak}.csv')

fig = px.line_mapbox(df, lat='latitude', lon='longitude', hover_name='elevation', zoom=9, custom_data=['elevation'])
fig.update_layout(mapbox_style='open-street-map')
fig.update_layout(mapbox_center_lon=22.859, mapbox_center_lat=53.073)
fig.update_traces(line=dict(color='rgb(160, 75, 196)'), hovertemplate='%{customdata[0]:.2f}m')

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

###################################
dzety = []

model = np.genfromtxt('2021-EVRF2007-NH.txt')
model2 = zmniejszenieSiatki(model, phis, lams)

for phi, lam in zip(phis, lams):
    dzeta = interpolacja(model2, phi, lam)
    dzety.append(dzeta)

df['HN'] = df['hel'] - np.array(dzety)

df['HN_diff'] = df['HN'] - df['elevation']

dif = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['HN_diff'], mode='lines', hovertemplate='Dystans: %{x:.2f} km<br>Różnica wysokości normalnych: %{y:.4f} m<extra></extra>'))
dif.update_layout(title='Różnica wysokości między 2021-EVRF2007-NH a 2011-KRON86-NH',
                    xaxis_title='Dystans [km]',
                  yaxis_title='Różnica wysokości[m]',
                  title_x=0.5)
dif.update_traces(line=dict(color='rgb(160, 75, 196)'))

#profil trasy - wysokość elipsoidalna
eli = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['hel'], mode='lines', hovertemplate='Dystans: %{x:.2f} km<br>Wysokość elipsoidalna: %{y:.2f} m<extra></extra>', fill='tozeroy'))
eli.update_layout(xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość elipsoidalna [m]',
                  margin=dict(t=0))
eli.update_yaxes(range=[100, 200])
eli.update_traces(line=dict(color='rgb(160, 75, 196)'), fillcolor='rgb(201, 133, 230)')

#profil trasy - 2021-EVRF2007-NH
ams = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['HN'], mode='lines', hovertemplate='Dystans: %{x:.2f} km<br>Wysokość normalna: %{y:.2f} m<extra></extra>', fill='tozeroy'))
ams.update_layout(xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość normalna [m]',
                  margin=dict(t=0))
ams.update_yaxes(range=[100, 200])
ams.update_traces(line=dict(color='rgb(160, 75, 196)'), fillcolor='rgb(201, 133, 230)')

# profil trasy - 2011-KRON86-NH
kro = go.Figure(data=go.Scatter(x=df['distance']/1000, y=df['elevation'], mode='lines', hovertemplate='Dystans: %{x:.2f} km<br>Wysokość normalna: %{y:.2f} m<extra></extra>', fill='tozeroy'))
kro.update_layout(xaxis_title='Dystans [km]',
                  yaxis_title='Wysokość normalna [m]',
                  margin=dict(t=0))
kro.update_yaxes(range=[100, 200])
kro.update_traces(line=dict(color='rgb(160, 75, 196)'), fillcolor='rgb(201, 133, 230)')


dlugosctrasy = "{:.2f}".format(df['dists'].sum() / 1000)
najnizszy = "{:.2f}".format(df['elevation'].min())
najwyzszy = "{:.2f}".format(df['elevation'].max())


table = go.Figure(data=[go.Table(header=dict(values=['Długość trasy [km]', 'Najwyższy punkt trasy [m]', 'Najniższy punkt trasy [m]']),
                 cells=dict(values=[[dlugosctrasy], [najwyzszy], [najnizszy]]))])
table.update_layout(margin=dict(t=0))

app.layout = html.Div(children=[
    html.H1(children='Zadanie 5', style={'text-align': 'center'}),
    html.H2(children='Obwodnica rowerowa po Narwiańskim Parku Narodowym', style={'text-align': 'center'}),
    html.Div(
        style={'display': 'inline-block', 'width': '48%'},
        children=[
            html.H3(children='Mapa trasy', style={'text-align': 'center', 'margin-bottom': '0px'}),
            dcc.Graph(
                id='map',
                figure=fig,
                style={'height': '60vh'}
            ),
            dcc.Graph(
                id='table',
                figure=table,
                style={'height': '20vh', 'margin': '0px'}
            )
        ]
    ),

    html.Div(
        style={'display': 'inline-block', 'width': '48%', 'vertical-align': 'top'},
        children=[
            html.H3(children='Proszę wybrać wykres:', style={'text-align': 'center'}),
            dcc.Dropdown(
                id='chart-dropdown',
                options=[
                    {'label': 'Wysokość elipsoidalna', 'value': 'eli'},
                    {'label': 'Wysokość normalna - 2021-EVRF2007-NH', 'value': 'ams'},
                    {'label': 'Wysokość normalna - 2011-KRON86-NH', 'value': 'kro'}
                ],
                value='ams',
                style = {'margin-bottom': '3vh', 'width': '80%', 'margin-left': '4vh'}
            ),
            dcc.Graph(
                id='selected-chart',
                style={'margin': '0px', 'height': '35vh'}
            ),
            dcc.Graph(
                id='dif',
                figure=dif,
                style={'margin': '0px', 'height': '35vh'}
            )
        ]
    ),
    html.Footer(children=[html.P(children=["Źródło trasy: ", html.A(children="GreenVelo", href="https://greenvelo.pl/detal/1217-greenvelo-obwodnica-rowerowa-po-narwianskim-parku-narodowym")])], style={'text-align': 'center', 'margin-top': '3vh'})
],
style={'margin': '0px', 'font-family': 'Arial'})

@app.callback(
    Output('selected-chart', 'figure'),
    [Input('chart-dropdown', 'value')]
)

def update_chart(selected_chart):
    if selected_chart == 'eli':
        return eli
    elif selected_chart == 'ams':
        return ams
    elif selected_chart == 'kro':
        return kro

if __name__ == '__main__':
    app.run_server(debug=True)


