#from lauren import *
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib
import sys
from matplotlib import colors
import seaborn as sns

##################### z = 0 #########################
fn_z0="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
fn_z02="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
fn_z05="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn_z1="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"

metal_cmap = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
rvir0,rvir02,rvir05,rvir1 = 316.86,294.55,226.32,162.66


pf = load(fn_z0, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(1000.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.set_cmap('all',metal_cmap)
pp.set_zlim('all',1e-4,1e-2)
pp.annotate_sphere(pos,radius=(rvir0,'kpc'),circle_args={'color':'white','ls':'--','lw':2})
pp.save('Metallicity_z0')

pf = load(fn_z02, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(1000.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.set_cmap('all',metal_cmap)
pp.set_zlim('all',1e-4,1e-2)
pp.annotate_sphere(pos,radius=(rvir02,'kpc'),circle_args={'color':'white','ls':'--','lw':2})
pp.save('Metallicity_z02')

pf = load(fn_z05, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(1000.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.set_cmap('all',metal_cmap)
pp.set_zlim('all',1e-4,1e-2)
pp.annotate_sphere(pos,radius=(rvir05,'kpc'),circle_args={'color':'white','ls':'--','lw':2})
pp.save('Metallicity_z05')

pf = load(fn_z1, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')
pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(1000.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
pp.set_cmap('all',metal_cmap)
pp.set_zlim('all',1e-4,1e-2)
pp.annotate_sphere(pos,radius=(rvir1,'kpc'),circle_args={'color':'white','ls':'--','lw':2})
pp.save('Metallicity_z1')

