from yt.mods import *
import numpy as np
import cPickle
import re

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
pos =  (0.38492584228516125, 0.46051788330078125, 0.49160003662109875)

##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
#val, pos = pf.h.find_max('Density')


#fields = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']
fields = ['Emission_OVI','Emission_CIII_977']

#ORIGINAL
#width = 981.5/pf['kpc']
res = [151,151]
#thickness = 1000./pf['kpc']

#ZOOM IN
width = 500./pf['kpc']
thickness = 500./pf['kpc']

project_X = False
project_Y = False
project_Z = True
project_faceon = False

if project_X:
	frbx = make_thin_projection(pf,'x',fields,pos,width,thickness,res)

	for field in fields:
		key = re.split('Emission_|',field)[1]
		fileout = 'bertone_frbs/emis/grid_galquas/forIMO/frbx_500kpc_z1_'+key+'.cpkl'
		cPickle.dump(frbx[field],open(fileout,'wb'),protocol=-1)

if project_Y:
	projy = pf.h.proj(1,fields)
	frby = projy.to_frb(width,res,center=pos)

	for field in fields:
      	        key = re.split('Emission_|',field)[1]
                fileout = 'bertone_frbs/emis/grid_galquas/forIMO/frby_500kpc_z1_'+key+'.cpkl'
                cPickle.dump(frby[field],open(fileout,'wb'),protocol=-1)
	

if project_Z:
	projz = pf.h.proj(2,fields)
	frbz = projz.to_frb(width,res,center=pos)

        for field in fields:
                key = re.split('Emission_|',field)[1]
                fileout = 'bertone_frbs/emis/grid_galquas/forIMO/frbz_500kpc_z1_satellite_'+key+'.cpkl'
                cPickle.dump(frbz[field],open(fileout,'wb'),protocol=-1)


def make_thin_frb(pf,axis,fields,center,width,thickness,resolution):
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
        obj = pf.h.proj(axis,fields,source=region,center=center)
        frb = obj.to_frb(width,resolution,center=center)

        return frb

