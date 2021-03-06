import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib
import sys

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 108.0/pf['kpc']

#data_source = pf.h.sphere(pos,rad)

#field = ['Density','HI_NumberDensity']
#field = ['Emission_SiIII_1883']
#field = ['H_NumberDensity','Temperature','Metallicity']
#field = ['RadialVelocityKMSABS']
#field = ['H_NumberDensity']
field = ['z-velocity']

font = {'weight':'bold','size':22,'family':'sans-serif'}
matplotlib.rc('font',**font)

#sp = SlicePlot(pf,'y',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)

sp = SlicePlot(pf,'z',field,center=pos,width=(320.,'kpc'))

#sp.annotate_particles(0.05,p_size=1.0,dm_only=True,alpha=0.03) #also stars_only keyword
sp.annotate_grids()
sp.set_font(font)
sp.set_cmap('all','spectral')
#pp.set_zlime('all',low,high)

#sp.save('003')
#sp.save('slice-test')
#sp.save('AMRgrid_basic_fields_slice')
sp.save('investigate')

