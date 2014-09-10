from yt.mods import *
import numpy as np
import cPickle

#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')

fields = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']

#width = 200./pf['kpc']
width = 320./pf['kpc']
res = [320,320]
#res = [800,800]
#res = [40,40]
#pos = [0.40328598,0.47176743,0.46131516]
pos = [0.39871597,0.46913528,0.46808243]

project_X = True
project_Y = True
project_Z = True

if project_X:
	projx = pf.h.proj(0,fields)
	frbx = projx.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_CIII.cpkl'
	cPickle.dump(frbx['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_CIII_977.cpkl'
	cPickle.dump(frbx['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_CIV.cpkl'
	cPickle.dump(frbx['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_OVI.cpkl'
	cPickle.dump(frbx['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_MgII.cpkl'
	cPickle.dump(frbx['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_SiII.cpkl'
	cPickle.dump(frbx['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_SiIII_1207.cpkl'
	cPickle.dump(frbx['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_SiIII_1883.cpkl'
	cPickle.dump(frbx['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_SiIV.cpkl'
	cPickle.dump(frbx['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbx_1kpc_z02_HAlpha.cpkl'
	cPickle.dump(frbx['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

if project_Y:
	projy = pf.h.proj(1,fields)
	frby = projy.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_CIII.cpkl'
	cPickle.dump(frby['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_CIII_977.cpkl'
	cPickle.dump(frby['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_CIV.cpkl'
	cPickle.dump(frby['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_OVI.cpkl'
	cPickle.dump(frby['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_MgII.cpkl'
	cPickle.dump(frby['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_SiII.cpkl'
	cPickle.dump(frby['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_SiIII_1207.cpkl'
	cPickle.dump(frby['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_SiIII_1883.cpkl'
	cPickle.dump(frby['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_SiIV.cpkl'
	cPickle.dump(frby['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frby_1kpc_z02_HAlpha.cpkl'
	cPickle.dump(frby['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

if project_Z:
	projz = pf.h.proj(2,fields)
	frbz = projz.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_CIII.cpkl'
	cPickle.dump(frbz['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_CIII_977.cpkl'
	cPickle.dump(frbz['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_CIV.cpkl'
	cPickle.dump(frbz['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_OVI.cpkl'
	cPickle.dump(frbz['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_MgII.cpkl'
	cPickle.dump(frbz['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_SiII.cpkl'
	cPickle.dump(frbz['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_SiIII_1207.cpkl'
	cPickle.dump(frbz['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_SiIII_1883.cpkl'
	cPickle.dump(frbz['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_SiIV.cpkl'
	cPickle.dump(frbz['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/control/frbz_1kpc_z02_HAlpha.cpkl'
	cPickle.dump(frbz['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

