#from lauren import *
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib
import sys

#field = str(sys.argv[1])
#print field

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

########THIS IS FOR Z=0 !!!!!!!!!!! ######################
pos = [0.40328598,0.47176743,0.46131516]
########THIS IS FOR Z=0.2 !!!!!!!!!!######################
#pos = [0.39871597,0.46913528,0.46808243]

#field = ['HI_NumberDensity']
field = ['OVI_Density','Emission_OVI']
#field = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density']

pp = ProjectionPlot(pf,'x',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
pp.set_cmap('all','spectral')
#pp.set_zlim('all',10**14.0,10**23.5)
#pp.set_zlim('all',10**12.0,10**24.0)
#pp.set_zlim('all',1.,10**7)
#pp.save('euvb_proj/ColDens/control/euvb_control_z02')
#pp.save('bertone_proj/ColDens/grid_galquas/g1q1/g1q1')
pp.save('check_Ximena')


