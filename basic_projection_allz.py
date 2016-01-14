#from lauren import *
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib
import sys
from matplotlib import colors

##################### z = 0 #########################
fn_z0="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
fn_z02="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
fn_z05="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn_z1="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"

pf = load(fn_z0, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(600.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.save('Metallicity_z0')

pf = load(fn_z02, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(600.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.save('Metallicity_z02')

pf = load(fn_z05, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(600.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.save('Metallicity_z05')

pf = load(fn_z1, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(600.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.save('Metallicity_z1')

