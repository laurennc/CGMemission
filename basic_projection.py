#from lauren import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *
import numpy as np
import sys


#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val, pos = pf.h.find_max('Density')
#pos = [0.40328598,0.47176743,0.46131516]
#rad = 108.0/pf['kpc']
rad = 320./pf['kpc']

#field = ['Emission_OVI','Emission_MgII','Emission_CIV','Emission_HAlpha']
#field = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']
#field = ['Emission_HAlpha']
#field = ["Density"]

#for i in "xyz":
#	pp = ProjectionPlot(pf,i,field,center=pos,width=(216.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
#	pp.set_cmap('all','spectral')
#pp.set_zlim('all',10**-5.0,10**2.0)
#pp.set_zlim('all',10**-2,10**7)
#	pp.save('visualize_axes')
#pp.save('bertone_proj/grid_galquas/g1q01')


#For making radial profiles of different quanitites:
pc = PlotCollection(pf)
#Default of this type of profile is weight='CellMassMsun', x_bins=128, x_log=True
pc.add_profile_sphere(320.0,"kpc",["Radiuskpc","Temperature"],x_bins=320,x_log=False)
pc.save()

pc = PlotCollection(pf)
pc.add_profile_sphere(320.0,"kpc",["Radiuskpc","H_NumberDensity"],x_bins=320,x_log=False)
pc.save()



