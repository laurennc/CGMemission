from yt.mods import *
import numpy as np
import cPickle

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
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"
##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

#fields = ['']
weightfield = 'Density'

width = 320./pf['kpc']
#res = [13,13]
#res = [64,64]
res = [320,320]
thickness = 500./pf['kpc']

project_X = True
project_Y = True
project_Z = True
project_faceon = False
velocity = False
temperature = False
HIdens = True

if project_X:
	if temperature:
		frbx = make_thin_frb(pf,'x',fields,width,thickness,res,weight_field=weightfield)	
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/''.cpkl'
        	cPickle.dump(frbx['Temperature'],open(fileout,'wb'),protocol=-1)

	if HIdens:
		projx = pf.h.proj(0,'HI_NumberDensity')
		frbx = projx.to_frb(width,res,center=pos)
		fileout = 'bertone_frbs/coldens/grid_galquas/g1q1/frbx_1kpc_z1_HIdens.cpkl'
		cPickle.dump(frbx['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

        if velocity:
		projx = pf.h.proj(0,'x-velocity',weight_field=weightfield)
	        frbx = projx.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbx_xvel_mass.cpkl'
		cPickle.dump(frbx['x-velocity'],open(fileout,'wb'),protocol=-1)

if project_Y:
	if temperature:
		projy = pf.h.proj(1,'Temperature',weight_field=weightfield)
        	frby = projy.to_frb(width,res,center=pos)
        	fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frby_temp_mass.cpkl'
        	cPickle.dump(frby['Temperature'],open(fileout,'wb'),protocol=-1)

        if HIdens:
                projy = pf.h.proj(1,'HI_NumberDensity')
                frby = projy.to_frb(width,res,center=pos)
                fileout = 'bertone_frbs/coldens/grid_galquas/g1q1/frby_1kpc_z1_HIdens.cpkl'
                cPickle.dump(frby['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

	if velocity:
     		projy = pf.h.proj(1,'y-velocity',weight_field=weightfield)
        	frby = projy.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frby_yvel_mass.cpkl'
        	cPickle.dump(frby['y-velocity'],open(fileout,'wb'),protocol=-1)

if project_Z:
	if temperature:
		projz = pf.h.proj(2,'Temperature',weight_field=weightfield)
        	frbz = projz.to_frb(width,res,center=pos)
        	fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbz_temp_mass.cpkl'
        	cPickle.dump(frbz['Temperature'],open(fileout,'wb'),protocol=-1)

        if HIdens:
                projz = pf.h.proj(2,'HI_NumberDensity')
                frbz = projz.to_frb(width,res,center=pos)
                fileout = 'bertone_frbs/coldens/grid_galquas/g1q1/frbz_1kpc_z1_HIdens.cpkl'
                cPickle.dump(frbz['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

	if velocity:
		projz = pf.h.proj(2,'z-velocity',weight_field=weightfield)
		frbz = projz.to_frb(width,res,center=pos)
		fileout = '/u/10/l/lnc2115/vega/data/Ryan/pickles/frbz_zvel_mass.cpkl'
		cPickle.dump(frbz['z-velocity'],open(fileout,'wb'),protocol=-1)


