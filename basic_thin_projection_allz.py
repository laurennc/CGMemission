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
fn_z0="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
pf = load(fn_z0, file_style="%s.grid.cpu%%04i")
file_begin = 'bertone_frbs/final/basic/z0/g1q1/'
rvir = 316.864810195 #kpc
#pos = [ 0.40320258,  0.4718025 ,  0.46123817]
val, pos = pf.h.find_max('Density')
width = 475./pf['kpc'] ##roughly 1.5*rvir
resolution = [475,475]

frbx = make_thin_frb(pf,'x','Temperature','Density',pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z0_Temperature_Densityweight.cpkl'
cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

frbx = make_thin_frb(pf,'x','Density',None,pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z0_Density.cpkl'
cPickle.dump(frbx['Density'],open(fileout,'wb'),protocol=-1)

##################### z = 0.2 ##########################
fn_z02="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn_z02, file_style="%s.grid.cpu%%04i")
file_begin = 'bertone_frbs/final/basic/z02/g1q1/'
rvir = 294.55286015 #kpc
#pos = [ 0.39898226,  0.46904468,  0.46830814]
val, pos = pf.h.find_max('Density')
width = 442./pf['kpc'] ##roughly 1.5*rvir
resolution = [442,442]

frbx = make_thin_frb(pf,'x','Temperature','Density',pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z02_Temperature_Densityweight.cpkl'
cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

frbx = make_thin_frb(pf,'x','Density',None,pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z02_Density.cpkl'
cPickle.dump(frbx['Density'],open(fileout,'wb'),protocol=-1)


####################  z = 0.5 ###########################
fn_z05="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
pf = load(fn_z05, file_style="%s.grid.cpu%%04i")
file_begin = 'bertone_frbs/final/basic/z05/g1q1/'
rvir = 226.321520623 #kpc
#pos = [ 0.3932969 ,  0.46559631,  0.47849673]
val, pos = pf.h.find_max('Density')
width = 340./pf['kpc'] ##roughly 1.5*rvir
resolution = [340,340]

frbx = make_thin_frb(pf,'x','Temperature','Density',pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z05_Temperature_Densityweight.cpkl'
cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

frbx = make_thin_frb(pf,'x','Density',None,pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z05_Density.cpkl'
cPickle.dump(frbx['Density'],open(fileout,'wb'),protocol=-1)



####################  z = 1.0 ##########################
fn_z1="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"
pf = load(fn_z1, file_style="%s.grid.cpu%%04i")
file_begin = 'bertone_frbs/final/basic/z1/g1q1/'
rvir = 162.656338662 #kpc
#pos = [ 0.38580719,  0.46065774,  0.49183819]
val, pos = pf.h.find_max('Density')
width = 245./pf['kpc'] ##roughly 1.5*rvir
resolution = [245,245]

frbx = make_thin_frb(pf,'x','Temperature','Density',pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z1_Temperature_Densityweight.cpkl'
cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

frbx = make_thin_frb(pf,'x','Density',None,pos,width,width,resolution)
fileout = file_begin+'frbx_1kpc_500kpc_z1_Density.cpkl'
cPickle.dump(frbx['Density'],open(fileout,'wb'),protocol=-1)





