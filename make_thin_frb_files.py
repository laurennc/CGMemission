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

fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

resolution = [320,320]
center = pos
thickness = 500./pf['kpc']
width = 320./pf['kpc']

fields = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density','HI_NumberDensity']

project_X = True
project_Y = True
project_Z = True

if project_X:
	axis = 'x'
	frbx = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_CIIIdens.cpkl'
	cPickle.dump(frbx['CIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_CIVdens.cpkl'
	cPickle.dump(frbx['CIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_OVIdens.cpkl'
	cPickle.dump(frbx['OVI_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_MgIIdens.cpkl'
	cPickle.dump(frbx['MgII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_SiIIdens.cpkl'
	cPickle.dump(frbx['SiII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_SiIIIdens.cpkl'
	cPickle.dump(frbx['SiIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_SiIVdens.cpkl'
	cPickle.dump(frbx['SiIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbx_1kpc_500kpc_z02_HIdens.cpkl'
	cPickle.dump(frbx['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

if project_Y:
	axis = 'y'
        frbx = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_CIIIdens.cpkl'
        cPickle.dump(frbx['CIII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_CIVdens.cpkl'
        cPickle.dump(frbx['CIV_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_OVIdens.cpkl'
        cPickle.dump(frbx['OVI_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_MgIIdens.cpkl'
        cPickle.dump(frbx['MgII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_SiIIdens.cpkl'
        cPickle.dump(frbx['SiII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_SiIIIdens.cpkl'
        cPickle.dump(frbx['SiIII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_SiIVdens.cpkl'
        cPickle.dump(frbx['SiIV_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frby_1kpc_500kpc_z02_HIdens.cpkl'
        cPickle.dump(frbx['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

if project_Z:
	axis = 'z'
        frbx = make_thin_frb(pf,axis,fields,center,width,thickness,resolution)

        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_CIIIdens.cpkl'
        cPickle.dump(frbx['CIII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_CIVdens.cpkl'
        cPickle.dump(frbx['CIV_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_OVIdens.cpkl'
        cPickle.dump(frbx['OVI_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_MgIIdens.cpkl'
        cPickle.dump(frbx['MgII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_SiIIdens.cpkl'
        cPickle.dump(frbx['SiII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_SiIIIdens.cpkl'
        cPickle.dump(frbx['SiIII_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_SiIVdens.cpkl'
        cPickle.dump(frbx['SiIV_Density'],open(fileout,'wb'),protocol=-1)
        fileout = 'bertone_frbs/coldens/grid_galquas/g1q01/frbz_1kpc_500kpc_z02_HIdens.cpkl'
        cPickle.dump(frbx['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)


