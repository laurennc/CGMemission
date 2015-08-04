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


#field = str(sys.argv[1])
#print field

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

########THIS IS FOR Z=0 !!!!!!!!!!! ######################
pos = [0.40328598,0.47176743,0.46131516]
########THIS IS FOR Z=0.2 !!!!!!!!!!######################
#pos = [0.39871597,0.46913528,0.46808243]
#thick = 500./pf['kpc']

#fields = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density']
#fields = ['H_NumberDensity','HI_NumberDensity','OVI_Density','CIV_Density','Emission_OVI','Emission_CIV']

#for field in fields:
#	pc = PlotCollection(pf,center=pos)
#	pc.add_thin_projection(field,'z',thickness=thick)
#	pc.set_cmap('spectral')
#	pc.set_width(500.,'kpc')
	#pc.set_zlim(10**12.0,10**24.0)
#	pc.save('diagnose/z0')


fields = ['H_NumberDensity','Temperature']
weights = ['CIII_Density','OVI_Density','SiIV_Density']
center = pos
width = 320./pf['kpc']
thickness = 500./pf['kpc']
resolution = [320,320]

file_begin = 'bertone_frbs/final/basic/z0/g1q01/'

for weight in weights:
	frbx = make_thin_frb(pf,'x',fields,weight,center,width,thickness,resolution)
	for field in fields:
		fileout = file_begin+'frbx_1kpc_500kpc_z0_'+field+'_'+weight+'.cpkl'
		cPickle.dump(frbx[field],open(fileout,'wb'),protocol=-1)




