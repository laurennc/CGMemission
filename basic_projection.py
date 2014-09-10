#from lauren import *
from yt.mods import *
import numpy
import numpy as np
import matplotlib
import sys


#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val, pos = pf.h.find_max('Density')
#pos = [0.40328598,0.47176743,0.46131516]
#rad = 108.0/pf['kpc']
rad = 320./pf['kpc']

#field = ['Emission_OVI','Emission_MgII','Emission_CIV','Emission_HAlpha']
field = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']

pp = ProjectionPlot(pf,'x',field,center=pos,width=(640.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
pp.set_cmap('all','spectral')
#pp.set_zlim('all',10**-5.0,10**2.0)
pp.set_zlim('all',10**-2,10**7)
pp.save('bertone_proj/grid_galquas/g1q01')



