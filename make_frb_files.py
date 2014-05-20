from lauren import *
import cPickle

fn="/media/caldisk/lauren/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#fields = ["HAlpha_Emissivity_R","Emission_CIV","Emission_OVI"]
fields = ["Emission_HAlpha","Emission_CIII","Emission_CIV","Emission_CIV","Emission_CIII_977"]

width = 200./pf['kpc']
#res = [800,800]
res = [40,40]
pos = [0.39871597,0.46913528,0.46808243]
dims = "xyz"

projx = pf.h.proj(0,fields)
frbx = projx.to_frb(width,res,center=pos)

fileout = 'frbx_5kpc_CIV.cpkl'
cPickle.dump(frbx['Emission_CIV'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_OVI.cpkl'
cPickle.dump(frbx['Emission_OVI'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_CIII.cpkl'
cPickle.dump(frbx['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_CIII_977.cpkl'
cPickle.dump(frbx['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_HAlpha.cpkl'
cPickle.dump(frbx['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

projy = pf.h.proj(1,fields)
frby = projy.to_frb(width,res,center=pos)

fileout = 'frby_5kpc_CIV.cpkl'
cPickle.dump(frby['Emission_CIV'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_OVI.cpkl'
cPickle.dump(frby['Emission_OVI'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_CIII.cpkl'
cPickle.dump(frby['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_CIII_977.cpkl'
cPickle.dump(frby['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_HAlpha.cpkl'
cPickle.dump(frby['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

projz = pf.h.proj(2,fields)
frbz = projz.to_frb(width,res,center=pos)

fileout = 'frbz_5kpc_CIV.cpkl'
cPickle.dump(frbz['Emission_CIV'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_OVI.cpkl'
cPickle.dump(frbz['Emission_OVI'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_CIII.cpkl'
cPickle.dump(frbz['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_CIII_977.cpkl'
cPickle.dump(frbz['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_HAlpha.cpkl'
cPickle.dump(frbz['Emission_HAlpha'],open(fileout,'wb'),protocol=-1)

