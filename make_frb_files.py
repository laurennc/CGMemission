from lauren import *
import cPickle

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#fields = ["HAlpha_Emissivity_R","Emission_CIV","Emission_OVI"]
#fields = ["Emission_CIV","Emission_OVI","Emission_CIV_Scaled","Emission_OVI_Scaled","Emission_CIV_Scaled_ncut","Emission_OVI_Scaled_ncut"]
#fields = ["Emission_CIII_977","Emission_CIII_977_Scaled","Emission_CIII_977_Scaled_ncut"]
fields = ["Emission_CIII","Emission_CIII_Scaled","Emission_CIII_Scaled_ncut"]


width = 200./pf['kpc']
#res = [800,800]
res = [40,40]
pos = [0.40328598,0.47176743,0.46131516]
dims = "xyz"

projx = pf.h.proj(0,fields)
frbx = projx.to_frb(width,res,center=pos)

#fileout = 'frbx_5kpc_CIV.cpkl'
#cPickle.dump(frbx['Emission_CIV'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_CIV_Scaled.cpkl'
#cPickle.dump(frbx['Emission_CIV_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_CIV_Scaled_ncut.cpkl'
#cPickle.dump(frbx['Emission_CIV_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_OVI.cpkl'
#cPickle.dump(frbx['Emission_OVI'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_OVI_Scaled.cpkl'
#cPickle.dump(frbx['Emission_OVI_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_OVI_Scaled_ncut.cpkl'
#cPickle.dump(frbx['Emission_OVI_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_CIII.cpkl'
cPickle.dump(frbx['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_CIII_Scaled.cpkl'
cPickle.dump(frbx['Emission_CIII_Scaled'],open(fileout,'wb'),protocol=-1)
fileout = 'frbx_5kpc_CIII_Scaled_ncut.cpkl'
cPickle.dump(frbx['Emission_CIII_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_CIII_977.cpkl'
#cPickle.dump(frbx['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_CIII_977_Scaled.cpkl'
#cPickle.dump(frbx['Emission_CIII_977_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbx_5kpc_CIII_977_Scaled_ncut.cpkl'
#cPickle.dump(frbx['Emission_CIII_977_Scaled_ncut'],open(fileout,'wb'),protocol=-1)

projy = pf.h.proj(1,fields)
frby = projy.to_frb(width,res,center=pos)

#fileout = 'frby_5kpc_CIV.cpkl'
#cPickle.dump(frby['Emission_CIV'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_CIV_Scaled.cpkl'
#cPickle.dump(frby['Emission_CIV_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_CIV_Scaled_ncut.cpkl'
#cPickle.dump(frby['Emission_CIV_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_OVI.cpkl'
#cPickle.dump(frby['Emission_OVI'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_OVI_Scaled.cpkl'
#cPickle.dump(frby['Emission_OVI_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_OVI_Scaled_ncut.cpkl'
#cPickle.dump(frby['Emission_OVI_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_CIII_977.cpkl'
#cPickle.dump(frby['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_CIII_977_Scaled.cpkl'
#cPickle.dump(frby['Emission_CIII_977_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frby_5kpc_CIII_977_Scaled_ncut.cpkl'
#cPickle.dump(frby['Emission_CIII_977_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_CIII.cpkl'
cPickle.dump(frby['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_CIII_Scaled.cpkl'
cPickle.dump(frby['Emission_CIII_Scaled'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_5kpc_CIII_Scaled_ncut.cpkl'
cPickle.dump(frby['Emission_CIII_Scaled_ncut'],open(fileout,'wb'),protocol=-1)

projz = pf.h.proj(2,fields)
frbz = projz.to_frb(width,res,center=pos)

#fileout = 'frbz_5kpc_CIV.cpkl'
#cPickle.dump(frbz['Emission_CIV'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_CIV_Scaled.cpkl'
#cPickle.dump(frbz['Emission_CIV_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_CIV_Scaled_ncut.cpkl'
#cPickle.dump(frbz['Emission_CIV_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_OVI.cpkl'
#cPickle.dump(frbz['Emission_OVI'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_OVI_Scaled.cpkl'
#cPickle.dump(frbz['Emission_OVI_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_OVI_Scaled_ncut.cpkl'
#cPickle.dump(frbz['Emission_OVI_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_CIII_977.cpkl'
#cPickle.dump(frbz['Emission_CIII_977'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_CIII_977_Scaled.cpkl'
#cPickle.dump(frbz['Emission_CIII_977_Scaled'],open(fileout,'wb'),protocol=-1)
#fileout = 'frbz_5kpc_CIII_977_Scaled_ncut.cpkl'
#cPickle.dump(frbz['Emission_CIII_977_Scaled_ncut'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_CIII.cpkl'
cPickle.dump(frbz['Emission_CIII'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_CIII_Scaled.cpkl'
cPickle.dump(frbz['Emission_CIII_Scaled'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_5kpc_CIII_Scaled_ncut.cpkl'
cPickle.dump(frbz['Emission_CIII_Scaled_ncut'],open(fileout,'wb'),protocol=-1)


#patt = 'frbs/frb%c_250pc_Ha.cpkl'
#for i in range(3):
#	fileout = patt%dim[i]
#	cPickle.dump(frbs[i]['HAlpha_Emissivity_R'],open(fileout,'wb'),protocol=-1)

