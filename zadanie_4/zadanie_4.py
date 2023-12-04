from pyproj import Proj, transform, CRS, Transformer
import math

phi = 52.2296756
lam = 21.0122287

output_proj = CRS.from_epsg(2180) # układ wyjściowy
input_proj = CRS.from_epsg(4326) # układ wejściowy

x_out, y_out = Transformer.from_proj(input_proj, output_proj).transform(phi, lam)

# fi, lambda -> do pl-92 i pl-2000

# jeżeli z bibliotek pythona to pl-92, pl-2000, utm, laea, lcc

#strefa odwzorowawcza odpowiednia dla punktu 1

#jeśli lam = 17 stopni, wyubieramy lambda 18 stopni (strefa 6)
# jeśli lam = 23,5 stopnia, wybieramy lambda 24 stopnie (strefa 8)
# jeśli lam = 19.5 stopnia, wybieramy lambda 21 stopni 

# w przypadku utm, wybieramy strefę 33/34

## 2 czesc zadania
# 1. obliczenie odleglosci miedzy punktami na plaszczyznie (pitagoras)
# sprawdzenie różnic między Sel-d (dflugosc miedzy punktamia na plaszczyznie a dlugosc miedzy punktami na elipsoidzie)
# 2. obliczenie redukcji dlugosci z ppodanego wzoru
# Sgk = d/m0
# r = Sgk * (...) Y W UKLADZIE gaussa-krugera
# 3. Sgk + r => porównanie z Sel -> ~0,00m

# zrobic to w ukladzie pl-92 i pl-2000

# redukcje azymutow dla chetnych - nie trzeba na 5

# mozna sprobowac policzyc pola powierzchni w laea - wzor 14

print(x_out, y_out)

def gk (fi, lam, lam0, a = 6378137, e2 = 0.00669438002290):
    b2 = a**2 * (1 - e2)
    eprim2 = (a**2 - b2) / b2

    delta_lam = lam - lam0
    t = math.tan(fi)
