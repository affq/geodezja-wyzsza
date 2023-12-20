import math
import numpy as np

e2 = 0.00669438002290
a = 6378137
b2 = a**2 * (1 - e2)
eprim2 = (a**2 - b2) / b2

A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
A2 = 3/8 * (e2 + e2**2/4 + 15*e2**3/128)
A4 = 15/256 * (e2**2 + 3*e2**3/4)
A6 = 35*e2**3/3072

def to_gk(fi, lam, lam0):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    lam0 = np.deg2rad(lam0)
    delta_lam = lam - lam0

    t = np.tan(fi)
    eta2 = eprim2 * np.cos(fi)**2
    N = a / np.sqrt(1-e2*np.sin(fi)**2)

    sigma = a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))
    
    xgk = sigma + (((delta_lam**2) / 2) * N * np.sin(fi) * np.cos(fi)) * (1 + ((delta_lam**2) / 12) * np.cos(fi)**2 * (5 - t**2 + 9*eta2 + 4*eta2**2) + ((delta_lam**4) / 360) * np.cos(fi)**4 * (61 - 58*t**2 + t**4 + 270*eta2 - 330*t**2*eta2))

    ygk = delta_lam * N * np.cos(fi) * (1 + ((delta_lam**2) / 6) * np.cos(fi)**2 * (1 - t**2 + eta2) + ((delta_lam**4) / 120) * np.cos(fi)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*t**2*eta2))

    return xgk, ygk


def from_gk(xgk, ygk, lam0):
    fi = xgk / (a * A0)
    sigma = a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))

    while True:
        fi1 = fi + (xgk - sigma) / (a * A0)

        N = a / math.sqrt(1-e2*np.sin(fi1)**2)
        M = a * (1 - e2) / math.sqrt(1-e2*np.sin(fi1)**2)**3
        t = math.tan(fi1)
        eta2 = eprim2 * np.cos(fi1)**2
        sigma = a * (A0 * fi1 - A2 * np.sin(2 * fi1) + A4 * np.sin(4 * fi1) - A6 * np.sin(6 * fi1))

        if abs(fi1 - fi) < (0.000001 / 3600):
            break

        fi = fi1

    fi = fi1 - ((ygk**2 * t) / (2 * M * N)) * (1 - ygk**2 / 12 * N**2 * (5 + 3*t**2 + eta2 - 9*t**2*eta2 - 4*eta2**2) + ygk**4 / 360 * N**4 * (61 + 90*t**2 + 45*t**4))

    lam = lam0 + (ygk / (N * np.cos(fi))) * (1 - ygk**2 / 6 * N**2 * (1 + 2*t**2 + eta2) + ygk**4 / 120 * N**4 * (5 + 28*t**2 + 24*t**4 + 6*eta2 + 8*t**2*eta2))

    return fi, lam

def to_2000(fi, lam, lam0):
    xgk, ygk = to_gk(fi, lam, lam0)
    m0 = 0.999923

    nr = 0
    
    if lam >= np.deg2rad(13.5) and lam < np.deg2rad(16.5):
        nr = 5
    elif lam >= np.deg2rad(16.5) and lam < np.deg2rad(19.5):
        nr = 6
    elif lam >= np.deg2rad(19.5) and lam < np.deg2rad(22.5):
        nr = 7
    elif lam >= np.deg2rad(22.5) and lam < np.deg2rad(25.5):
        nr = 8

    x2000 = m0 * xgk
    y2000 = m0 * ygk + nr * 1000000 + 500000

    return x2000, y2000


def to_1992(fi, lam, lam0):
    xgk, ygk = to_gk(fi, lam, lam0)
    m0 = 0.9993

    x1992 = m0 * xgk - 5300000
    y1992 = m0 * ygk + 500000

    return x1992, y1992


if __name__ == "__main__":
    points = [[51.0, 19.0], [51.35954501388889, 19.0], [51.350750175, 20.435489580555554], [50.991204616666664, 20.435489580555554]]

    gk = []
    pl2000 = []
    pl1992 = []

    for point in points:
        fi = point[0]
        lam = point[1]
        lam0 = 19

        xgk, ygk = to_gk(fi, lam, lam0)
        gk.append([xgk, ygk])
        print("GK: ", xgk, ygk)

        x2000, y2000 = to_2000(fi, lam, lam0)
        pl2000.append([x2000, y2000])
        print("2000: ", x2000, y2000)

        x1992, y1992 = to_1992(fi, lam, lam0)
        pl1992.append([x1992, y1992])
        print("1992: ", x1992, y1992)
        print("")
    
        
    print("Długości odcinków między punktami na płaszczyźnie PL2000:")
    print("1-2: ", math.sqrt((pl2000[0][0] - pl2000[1][0])**2 + (pl2000[0][1] - pl2000[1][1])**2))
    print("2-3: ", math.sqrt((pl2000[1][0] - pl2000[2][0])**2 + (pl2000[1][1] - pl2000[2][1])**2))
    print("3-4: ", math.sqrt((pl2000[2][0] - pl2000[3][0])**2 + (pl2000[2][1] - pl2000[3][1])**2))
    print("4-1: ", math.sqrt((pl2000[3][0] - pl2000[0][0])**2 + (pl2000[3][1] - pl2000[0][1])**2))
    print("")

    length_gk = []
    for i in range(0, len(gk)):
        if i == len(gk)-1:
            length_gk.append(math.sqrt((gk[i][0] - gk[0][0])**2 + (gk[i][1] - gk[0][1])**2))
        else:
            length_gk.append(math.sqrt((gk[i][0] - gk[i+1][0])**2 + (gk[i][1] - gk[i+1][1])**2))

    print("Długości odcinków między punktami na płaszczyźnie G-K:")
    print("1-2: ", math.sqrt((gk[0][0] - gk[1][0])**2 + (gk[0][1] - gk[1][1])**2))
    print("2-3: ", math.sqrt((gk[1][0] - gk[2][0])**2 + (gk[1][1] - gk[2][1])**2))
    print("3-4: ", math.sqrt((gk[2][0] - gk[3][0])**2 + (gk[2][1] - gk[3][1])**2))
    print("4-1: ", math.sqrt((gk[3][0] - gk[0][0])**2 + (gk[3][1] - gk[0][1])**2))
    print("")

    print("Długości odcinków między punktami na płaszczyźnie 1992:")
    print("1-2: ", math.sqrt((pl1992[0][0] - pl1992[1][0])**2 + (pl1992[0][1] - pl1992[1][1])**2))
    print("2-3: ", math.sqrt((pl1992[1][0] - pl1992[2][0])**2 + (pl1992[1][1] - pl1992[2][1])**2))
    print("3-4: ", math.sqrt((pl1992[2][0] - pl1992[3][0])**2 + (pl1992[2][1] - pl1992[3][1])**2))
    print("4-1: ", math.sqrt((pl1992[3][0] - pl1992[0][0])**2 + (pl1992[3][1] - pl1992[0][1])**2))
    print("")


    middle_phis = []
    middle_ms = []
    middle_ns = []
    
    for i in range(0, len(points)):
        if i == len(points)-1:
            middle_phis.append((points[i][0] + points[0][0]) / 2)
        else:
            middle_phis.append((points[i][0] + points[i+1][0]) / 2)

    middle_phis = np.deg2rad(middle_phis)
    
    for phi in middle_phis:
        middle_ms.append(a * (1 - e2) / np.sqrt(1-e2*np.sin(phi)**2)**3)
        middle_ns.append(a / np.sqrt(1-e2*np.sin(phi)**2))

    reductions = []
    for i in range(0, len(points)):
        if i == len(points)-1:
            rAB = length_gk[i] * (gk[i][1]**2 + gk[i][1]*gk[0][1] + gk[0][1]**2) / (6 * middle_ms[i] * middle_ns[i])
        else:
            rAB = length_gk[i] * (gk[i][1]**2 + gk[i][1]*gk[i+1][1] + gk[i+1][1]**2) / (6 * middle_ms[i] * middle_ns[i])
        reductions.append(rAB)
    
    print("Redukcje długości:")
    print("1-2: ", reductions[0])
    print("2-3: ", reductions[1])
    print("3-4: ", reductions[2])
    print("4-1: ", reductions[3])
    print("")

    lengths = []
    for i in range(0, len(points)):
        if i == len(points)-1:
            lengths.append(length_gk[i] - reductions[i])
        else:
            lengths.append(length_gk[i] - reductions[i])
    
    print("Długości odcinków na elipsoidzie:")
    print("1-2: ", lengths[0])
    print("2-3: ", lengths[1])
    print("3-4: ", lengths[2])
    print("4-1: ", lengths[3])
    print("")


#pl2000
polepl2000 = (pl2000[1][0] * (pl2000[2][1] - pl2000[0][1]) + pl2000[2][0] * (pl2000[3][1] - pl2000[1][1]) + pl2000[3][0] * (pl2000[0][1] - pl2000[2][1]) + pl2000[0][0] * (pl2000[1][1] - pl2000[3][1])) / 2

print("Pole powierzchni na płaszczyźnie PL-2000: ", polepl2000)

#pl1992
polepl1992 = (pl1992[1][0] * (pl1992[2][1] - pl1992[0][1]) + pl1992[2][0] * (pl1992[3][1] - pl1992[1][1]) + pl1992[3][0] * (pl1992[0][1] - pl1992[2][1]) + pl1992[0][0] * (pl1992[1][1] - pl1992[3][1])) / 2

print("Pole powierzchni na płaszczyźnie PL-1992: ", polepl1992)


from pyproj import CRS, Transformer

output_proj = CRS.from_epsg(3035)
input_proj = CRS.from_epsg(4326)

transformer = Transformer.from_crs(input_proj, output_proj)

plalea = []
for point in points:
    plalea.append(transformer.transform(point[0], point[1]))

polepllaea = (plalea[1][0] * (plalea[2][1] - plalea[0][1]) + plalea[2][0] * (plalea[3][1] - plalea[1][1]) + plalea[3][0] * (plalea[0][1] - plalea[2][1]) + plalea[0][0] * (plalea[1][1] - plalea[3][1])) / 2

print("Pole powierzchni na płaszczyźnie PL-LAEA: ", polepllaea)



