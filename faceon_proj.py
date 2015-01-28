import matplotlib
matplotlib.use('Agg')
from lauren import *
import matplotlib.pyplot as plt

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 25.0/pf['kpc']
ds = pf.h.sphere(pos,rad)

dense_ad = ds.cut_region(['grid["Density"] > 1e-24'])
dense_cold= dense_ad.cut_region(['grid["Temperature"] < 1e4'])

L = dense_cold.quantities['AngularMomentumVector']()

proj = OffAxisProjectionPlot(pf, L, "Density", ds.center, (320, "kpc"))

frb = proj["Density"].image.get_array()

proj.save()



