from yt.mods import *
import numpy as np
import cPickle

##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"


##############################################################3
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
val, pos = pf.h.find_max('Density')


fields = ['Emission_HAlpha','Emission_CIV','Emission_OVI','Emission_CIII_977','Emission_CIII','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']
#fields = ['Emission_HAlpha','Emission_OVI','Emission_SiII','Emission_SiIII_1207','Emission_SiIII_1883','Emission_SiIV','Emission_MgII']


#width = 200./pf['kpc']
width = 320./pf['kpc']
#res = [13,13]
#res = [64,64]
res = [320,320]
#res = [800,800]
#res = [40,40]

project_X = True
project_Y = True
project_Z = True
project_faceon = False

if project_X:
	projx = pf.h.proj(0,fields)
	frbx = projx.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_CIII.cpkl'
	cPickle.dump(frbx['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_CIII_977.cpkl'
	cPickle.dump(frbx['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_CIV.cpkl'
	cPickle.dump(frbx['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_OVI.cpkl'
	cPickle.dump(frbx['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_MgII.cpkl'
	cPickle.dump(frbx['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_SiII.cpkl'
	cPickle.dump(frbx['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_SiIII_1207.cpkl'
	cPickle.dump(frbx['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_SiIII_1883.cpkl'
	cPickle.dump(frbx['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_SiIV.cpkl'
	cPickle.dump(frbx['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbx_1kpc_z05_HAlpha.cpkl'
	cPickle.dump(frbx['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

if project_Y:
	projy = pf.h.proj(1,fields)
	frby = projy.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_CIII.cpkl'
	cPickle.dump(frby['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_CIII_977.cpkl'
	cPickle.dump(frby['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_CIV.cpkl'
	cPickle.dump(frby['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_OVI.cpkl'
	cPickle.dump(frby['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_MgII.cpkl'
	cPickle.dump(frby['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_SiII.cpkl'
	cPickle.dump(frby['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_SiIII_1207.cpkl'
	cPickle.dump(frby['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_SiIII_1883.cpkl'
	cPickle.dump(frby['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_SiIV.cpkl'
	cPickle.dump(frby['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frby_1kpc_z05_HAlpha.cpkl'
	cPickle.dump(frby['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

if project_Z:
	projz = pf.h.proj(2,fields)
	frbz = projz.to_frb(width,res,center=pos)
	
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_CIII.cpkl'
	cPickle.dump(frbz['Emission_CIII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_CIII_977.cpkl'
	cPickle.dump(frbz['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_CIV.cpkl'
	cPickle.dump(frbz['Emission_CIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_OVI.cpkl'
	cPickle.dump(frbz['Emission_OVI'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_MgII.cpkl'
	cPickle.dump(frbz['Emission_MgII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_SiII.cpkl'
	cPickle.dump(frbz['Emission_SiII'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_SiIII_1207.cpkl'
	cPickle.dump(frbz['Emission_SiIII_1207'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_SiIII_1883.cpkl'
	cPickle.dump(frbz['Emission_SiIII_1883'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_SiIV.cpkl'
	cPickle.dump(frbz['Emission_SiIV'],open(fileout,'wb'),protocol=-1)
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbz_1kpc_z05_HAlpha.cpkl'
	cPickle.dump(frbz['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

if project_faceon:
	rad = 25.0/pf['kpc']
	ds = pf.h.sphere(pos,rad)
	
	dense_ad = ds.cut_region(['grid["Density"] > 1e-24'])
	dense_cold= dense_ad.cut_region(['grid["Temperature"] < 1e4'])
	
	L = dense_cold.quantities['AngularMomentumVector']()
	
	projFace = OffAxisProjectionPlot(pf, L, fields, ds.center, (320, "kpc"))
	
	#frbFace = projFace['Emission_CIII'].image.get_array()
        #fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_CIII.cpkl'
        #cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	
	#frbFace = projFace['Emission_CIII_977'].image.get_array()
	#fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_CIII_977.cpkl'
        #cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
	
	#frbFace = projFace['Emission_CIV'].image.get_array()
        #fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_CIV.cpkl'
        #cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_OVI'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_OVI.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_MgII'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_MgII.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_SiII'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_SiII.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_SiIII_1207'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_SiIII_1207.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_SiIII_1883'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_SiIII_1883.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_SiIV'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_SiIV.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)
        
	frbFace = projFace['Emission_HAlpha'].image.get_array()
	fileout = 'bertone_frbs/emis/grid_galquas/z05/g1q1/frbFace_1kpc_z05_HAlpha.cpkl'
        cPickle.dump(frbFace,open(fileout,'wb'),protocol=-1)



