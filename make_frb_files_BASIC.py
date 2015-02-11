from yt.mods import *
import numpy as np
import cPickle

##################### z = 0 #########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"
##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

#fields = ['']
weightfield = 'Density'

width = 320./pf['kpc']
#res = [13,13]
res = [64,64]
#res = [320,320]

project_X = True
project_Y = True
project_Z = True
project_faceon = False
velocity = False

if project_X:
	projx = pf.h.proj(0,'Temperature',weight_field=weightfield)
	frbx = projx.to_frb(width,res,center=pos)
	fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbx_temp_mass.cpkl'
        cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

        if velocity:
		projx = pf.h.proj(0,'x-velocity',weight_field=weightfield)
	        frbx = projx.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbx_xvel_mass.cpkl'
		cPickle.dump(frbx['x-velocity'],open(fileout,'wb'),protocol=-1)

if project_Y:
	projy = pf.h.proj(1,'Temperature',weight_field=weightfield)
        frby = projy.to_frb(width,res,center=pos)
        fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frby_temp_mass.cpkl'
        cPickle.dump(frby['Temperature'],open(fileout,'wb'),protocol=-1)

	if velocity:
     		projy = pf.h.proj(1,'y-velocity',weight_field=weightfield)
        	frby = projy.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frby_yvel_mass.cpkl'
        	cPickle.dump(frby['y-velocity'],open(fileout,'wb'),protocol=-1)

if project_Z:
	projz = pf.h.proj(2,'Temperature',weight_field=weightfield)
        frbz = projz.to_frb(width,res,center=pos)
        fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbz_temp_mass.cpkl'
        cPickle.dump(frbz['Temperature'],open(fileout,'wb'),protocol=-1)

	if velocity:
		projz = pf.h.proj(2,'z-velocity',weight_field=weightfield)
		frbz = projz.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbz_zvel_mass.cpkl'
		cPickle.dump(frbz['z-velocity'],open(fileout,'wb'),protocol=-1)


