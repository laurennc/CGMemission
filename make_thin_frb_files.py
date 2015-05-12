from yt.mods import *
import numpy as np
import cPickle
import re

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

##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"


pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

resolution = [320,320]
resolution = [64,64]
#resolution = [13,13]

center = pos
#thickness = 320./pf['kpc']
#thickness = 500./pf['kpc']
#width = 320./pf['kpc']
thickness = 1./pf['Mpc']
width = 1./pf['Mpc']


emission = True
coldens = False

if coldens:	
	#fields = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density']#,'HI_NumberDensity']
	fields = ['HI_NumberDensity']	

	file_begin = 'bertone_frbs/coldens/grid_galquas/g1q10/'
	key_string,keydx = '_Density',0
	file_end = 'dens.cpkl'

if emission:
	#fields = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']
	#file_begin = 'bertone_frbs/emis/grid_galquas/z02/g1q10/'	
	fields = ['Emission_OVI']
	file_begin = 'wide_'
	key_string,keydx = 'Emission_|',1
	file_end = '.cpkl'

project_X = True
project_Y = False
project_Z = False

if project_X:
	axis = 'x'
	frbx = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

	for field in fields:
                key = re.split(key_string,field)[keydx]
                fileout = file_begin + 'frbx_1kpc_500kpc_z02_'+key+file_end
                cPickle.dump(frbx[field],open(fileout,'wb'),protocol=-1)

if project_Y:
	axis = 'y'
        frby = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

	for field in fields:
                key = re.split(key_string,field)[keydx]
                fileout = file_begin + 'frby_1kpc_500kpc_z02_'+key+file_end
                cPickle.dump(frby[field],open(fileout,'wb'),protocol=-1)

if project_Z:
	axis = 'z'
        frbz = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

        for field in fields:
                key = re.split(key_string,field)[keydx]
                fileout = file_begin + 'frbz_1kpc_500kpc_z02_'+key+file_end
                cPickle.dump(frbz[field],open(fileout,'wb'),protocol=-1)




