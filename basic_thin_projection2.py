#from lauren import *
from yt.mods import *
import numpy as np
import matplotlib
import sys

#field = str(sys.argv[1])
#print field

#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

########THIS IS FOR Z=0 !!!!!!!!!!! ######################
#pos = [0.40328598,0.47176743,0.46131516]
########THIS IS FOR Z=0.2 !!!!!!!!!!######################
pos = [0.39871597,0.46913528,0.46808243]
thick = 500./pf['kpc']

fields = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density']
for field in fields:
	pc = PlotCollection(pf,center=pos)
	pc.add_thin_projection(field,'x',thickness=thick)
	pc.set_cmap('spectral')
	pc.set_width(500.,'kpc')
	pc.set_zlim(10**12.0,10**24.0)
	pc.save('euvb_proj/ColDens/control/euvb_control_z02_THIN')



