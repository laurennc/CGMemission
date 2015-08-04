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
####################  z = 1.0 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"


################satellite position at z=1#################
#pos =  (0.38492584228516125, 0.46051788330078125, 0.49160003662109875)

##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')


#field = ['HI_NumberDensity']
#field = ['OVI_Density','CIII_Density','SiIV_Density']
field = ['Emission_OVI_Scaled','Emission_CIII_977_Scaled','Emission_CIV_Scaled','Emission_SiIV_Scaled']
#field = ['OVI_Density','Emission_OVI']
#field = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density']

#pp = ProjectionPlot(pf,'x',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')

####LYMAN ALPHA COLORS
#cmap = colors.ListedColormap(['Gray','DarkViolet','DarkTurquoise','Chartreuse','HotPink'])
#bounds = [np.power(10.,14.),np.power(10.,15.5),np.power(10.,17.2),np.power(10.,19),np.power(10.,20.3),np.power(10.,23.)]
#norm = colors.BoundaryNorm(bounds,cmap.N)

####EMISSION COLORS
cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
bounds = [np.power(10.,-10.),np.power(10.,1.),np.power(10.,2.),np.power(10.,3.),np.power(10.,10.)]
norm = colors.BoundaryNorm(bounds,cmap.N)

pp = ProjectionPlot(pf,'x',field,center=pos,width=(320.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,cmap=cmap,norm=norm)

#cmap = 'YlGnBu'
#cmap = 'RdPu'
pp.set_cmap('all',cmap)
#pp.set_zlim('all',10**14.0,10**21.)
#pp.set_zlim('all',10**12.0,10**24.0)
#pp.set_zlim('all',10**12.0,10**16.0)
pp.set_zlim('all',1.,10**7)
#pp.save('euvb_proj/ColDens/control/euvb_control_z02')
#pp.save('bertone_proj/ColDens/grid_galquas/g1q1/g1q1')
pp.save('paper1_z1_probcolor')


