#from lauren import *
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib
import sys
from matplotlib import colors

##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 0.75 #########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0043/redshift0043"
####################  z = 1.0 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"


################satellite position at z=1#################
#pos =  (0.38492584228516125, 0.46051788330078125, 0.49160003662109875)

##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')



pp = ProjectionPlot(pf,'x','Density',center=pos,width=(600.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,cmap=cmap,norm=norm)
cmap = 'gnuplot2'
pp.set_cmap('all',cmap)
pp.save('density_z075_wide')

#field = ['Temperature','Metallicity']
#pp = ProjectionPlot(pf,'x',field,center=pos,width=(320.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.,weight_field='Density')#,cmap=cmap,norm=norm)
#cmap = 'spectral'
#pp.set_cmap('all',cmap)
#pp.save('paper1_z1')


