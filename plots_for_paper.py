import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *
import numpy as np
import sys
import matplotlib.colorbar as cb

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
orient = 'horizontal'

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

fig,axes,colorbars = get_multi_plot(3,1,colorbar=orient,bw=4)

pc = PlotCollection(pf,pos)

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Temperature","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**4.,10**7.),figure=fig,axes=[0,1])
p.set_cmap('spectral')

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Metallicity","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**-5.,10),figure=fig,axes=[0,2])
p.set_cmap('spectral')

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Density","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**-32.,10**-22.),figure=fig,axes=[0,0])
p.set_cmap('gnuplot2')

for p,cax in zip(pc.plots,colorbars):
	cbar = cb.Colorbar(cax,p.image,orientation=orient)
	p.colorbar = cbar
	p._autoset_label()

fig.savefig('paper1_sim_basic_profiles.png')



