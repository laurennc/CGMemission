import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *
import numpy as np
import sys
import matplotlib.colorbar as cb


##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"



orient = 'horizontal'

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

#fig,axes,colorbars = get_multi_plot(3,1,colorbar=orient,bw=4)

pc = PlotCollection(pf,pos)

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Temperature","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**4.,10**7.))
p.set_cmap('spectral')

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Metallicity","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**-5.,10))
p.set_cmap('spectral')

p = pc.add_phase_sphere(1000.,'kpc',["Radiuskpc","Density","CellMassMsun"],weight=None,x_log=True,x_bounds=(1.,1000.),y_log=True,y_bounds=(10**-32.,10**-22.))
p.set_cmap('gnuplot2')


pc.save('paper1_z1')



