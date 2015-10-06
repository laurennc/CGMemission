#from lauren import *
from yt.mods import *
import numpy as np
import matplotlib
import sys
import cPickle

def make_thin_frb(pf,axis,fields,weight,center,width,thickness,resolution):
        axis = fix_axis(axis)
        center = np.array(center)

        LE,RE = center.copy(),center.copy()
        LE[axis] -= thickness/2.0
        RE[axis] += thickness/2.0

        area_axes = [0,1,2]
        i = area_axes.index(axis)
        del area_axes[i]

        LE[area_axes] -= width/2.0
        RE[area_axes] += width/2.0

        region = pf.h.region(center, LE, RE)
        obj = pf.h.proj(axis,fields,weight_field=weight,source=region,center=center)
        frb = obj.to_frb(width,resolution,center=center)

        return frb



##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"


pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')



fields = ['H_NumberDensity','Temperature']
#weights = ['CIII_Density','OVI_Density','SiIV_Density']
#weights = ['Density']
weights = ['Emission_CIII_977','Emission_CIV','Emission_SiIV','Emission_OVI']
center = pos
width = 320./pf['kpc']
thickness = 500./pf['kpc']
resolution = [320,320]

file_begin = 'bertone_frbs/final/basic/z1/g1q1/'

for weight in weights:
	frbx = make_thin_frb(pf,'x',fields,weight,center,width,thickness,resolution)
	for field in fields:
		fileout = file_begin+'frbx_1kpc_500kpc_z1_'+field+'_'+weight+'.cpkl'
		cPickle.dump(frbx[field],open(fileout,'wb'),protocol=-1)




