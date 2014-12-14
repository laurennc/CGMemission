from yt.mods import *
import numpy as np
import cPickle

#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i")

fields = ['CIII_Density','CIV_Density','OVI_Density','MgII_Density','SiII_Density','SiIII_Density','SiIV_Density','HI_NumberDensity']

#width = 200./pf['kpc']
width = 320./pf['kpc']
res = [320,320]
#res = [800,800]
#res = [40,40]
#pos = [0.40328598,0.47176743,0.46131516]
pos = [0.39871597,0.46913528,0.46808243]
#dims = "xyz"

project_X = False
project_Y = False
project_Z = True

if project_X:
	projx = pf.h.proj(0,fields)
	frbx = projx.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_CIIIdens.cpkl'
	cPickle.dump(frbx['CIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_CIVdens.cpkl'
	cPickle.dump(frbx['CIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_OVIdens.cpkl'
	cPickle.dump(frbx['OVI_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_MgIIdens.cpkl'
	cPickle.dump(frbx['MgII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_SiIIdens.cpkl'
	cPickle.dump(frbx['SiII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_SiIIIdens.cpkl'
	cPickle.dump(frbx['SiIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_SiIVdens.cpkl'
	cPickle.dump(frbx['SiIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbx_1kpc_z02_HIdens.cpkl'
	cPickle.dump(frbx['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

if project_Y:
	projy = pf.h.proj(1,fields)
	frby = projy.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_CIIIdens.cpkl'
	cPickle.dump(frby['CIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_CIVdens.cpkl'
	cPickle.dump(frby['CIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_OVIdens.cpkl'
	cPickle.dump(frby['OVI_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_MgIIdens.cpkl'
	cPickle.dump(frby['MgII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_SiIIdens.cpkl'
	cPickle.dump(frby['SiII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_SiIIIdens.cpkl'
	cPickle.dump(frby['SiIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_SiIVdens.cpkl'
	cPickle.dump(frby['SiIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frby_1kpc_z02_HIdens.cpkl'
	cPickle.dump(frby['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

if project_Z:
	projz = pf.h.proj(2,fields)
	frbz = projz.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_CIIIdens.cpkl'
	cPickle.dump(frbz['CIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_CIVdens.cpkl'
	cPickle.dump(frbz['CIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_OVIdens.cpkl'
	cPickle.dump(frbz['OVI_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_MgIIdens.cpkl'
	cPickle.dump(frbz['MgII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_SiIIdens.cpkl'
	cPickle.dump(frbz['SiII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_SiIIIdens.cpkl'
	cPickle.dump(frbz['SiIII_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_SiIVdens.cpkl'
	cPickle.dump(frbz['SiIV_Density'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/coldens/grid_galquas/g1q0.05/frbz_1kpc_z02_HIdens.cpkl'
	cPickle.dump(frbz['HI_NumberDensity'],open(fileout,'wb'),protocol=-1)

