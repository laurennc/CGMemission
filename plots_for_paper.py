import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *
import numpy as np
import sys

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val, pos = pf.h.find_max('Density')
#pos = [0.40328598,0.47176743,0.46131516]
#rad = 108.0/pf['kpc']
#rad = 320./pf['kpc']

pc = PlotCollection(pf,pos)
p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Temperature","CellMassMsun"],weight=None)
pc.set_cmap('spectral')
pc.save('paper1')

pc = PlotCollection(pf,pos)
p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Metallicity","CellMassMsun"],weight=None,cmap='spectral',x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(-5,10),axes=[])
pc.set_cmap('spectral')
pc.save('paper1')

pc = PlotCollection(pf,pos)
p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Density","CellMassMsun"],weight=None)
pc.set_cmap('gnuplot2')
pc.save('paper1')



