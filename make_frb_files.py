from lauren import *
import cPickle

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

fields = ["HAlpha_Emissivity_R","Emission_CIV","Emission_OVI"]
width = 200./pf['kpc']
res = [800,800]
pos = [0.40328598,0.47176743,0.46131516]
dims = "xyz"

projx = pf.h.proj(0,fields)
frbx = projx.to_frb(width,res,center=pos)

projy = pf.h.proj(1,fields)
frby = projy.to_frb(width,res,center=pos)

projz = pf.h.proj(2,fields)
frbz = projz.to_frb(width,res,center=pos)

fileout = 'frbx_250pc_Ha.cpkl'
cPickle.dump(frbx['HAlpha_Emissivity_R'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_250pc_Ha.cpkl'
cPickle.dump(frby['HAlpha_Emissivity_R'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_250pc_Ha.cpkl'
cPickle.dump(frbz['HAlpha_Emissivity_R'],open(fileout,'wb'),protocol=-1)

fileout = 'frbx_250pc_CIV.cpkl'
cPickle.dump(frbx['Emission_CIV'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_250pc_CIV.cpkl'
cPickle.dump(frby['Emission_CIV'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_250pc_CIV.cpkl'
cPickle.dump(frbz['Emission_CIV'],open(fileout,'wb'),protocol=-1)

fileout = 'frbx_250pc_OVI.cpkl'
cPickle.dump(frbx['Emission_OVI'],open(fileout,'wb'),protocol=-1)
fileout = 'frby_250pc_OVI.cpkl'
cPickle.dump(frby['Emission_OVI'],open(fileout,'wb'),protocol=-1)
fileout = 'frbz_250pc_OVI.cpkl'
cPickle.dump(frbz['Emission_OVI'],open(fileout,'wb'),protocol=-1)










#patt = 'frbs/frb%c_250pc_Ha.cpkl'
#for i in range(3):
#	fileout = patt%dim[i]
#	cPickle.dump(frbs[i]['HAlpha_Emissivity_R'],open(fileout,'wb'),protocol=-1)

