import numpy as np
import gpxpy
import gpxpy.gpx
import pandas as pd

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