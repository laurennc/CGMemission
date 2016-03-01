## For a single pixel I want:
## NHI, NCIII, NCIV, NOVI, NSiIV
##      ECIII, ECIV, EOVI, ESiIV

import matplotlib.pyplot as plt
import numpy as np
import cPickle

count = 1
model_gqs = ['g1q01','g1q1','g1q10']

ions_coldens = ['CIII','CIV','OVI']
ions_emis = ['CIII_977','CIV','OVI']

model_beg_coldens = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/coldens/z02/'
model_beg_emis = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/z02/'
model_mid = '/frbx_1kpc_500kpc_z02_'


data = {}
fileout = 'absemis_test.cpkl'

xL = np.linspace(-160,160,320)
maxnow = 160.
xL, yL = np.meshgrid(xL,xL)
r = abs(xL+1j*yL)
dr = np.abs([xL[0,0] - xL[0,1]])
radial = np.arange(maxnow/dr)*dr + dr/2
nrad = len(radial)



## just for now I'm going to make a cut such that the CIII emission is above our detection limit
CIIIemis = model_beg_emis+model_gqs[count]+model_mid+'CIII_977.cpkl'
minrad = 25.
maxrad = 100.
thisindex = (r>=minrad) * (r<maxrad)
CIIIemis_rselect = CIIIemis[thisindex]
r_rselect = r[thisindex]
idx = np.where(np.log10(CIIIemis_rselect) > 1.2)[0]

## then just for this I'll pick index 100!

for ion in ions_emis:
	modelname_emis = model_beg_emis+model_gqs[count]+model_mid+ion+'.cpkl'
	emisdata = cPickle.load(open(modelname_emis,'rb'))
	emisdata = emisdata[thisindex]
	data['emis_'+ion] = emisdata[idx[100]]

for ion in ions_coldens:
	modelname_coldens = model_beg_coldens+model_gqs[count]+model_mid+ion+'dens.cpkl'
	coldensdata = cPickle.load(open(modelname_coldens,'rb'))
	coldensdata = coldensdata[thisindex]
	data['coldens_'+ion] = coldensdata[idx[100]]

modelname_HI = model_beg_coldens+model_gqs[count]+model_mid+'HI_NumberDensitydens.cpkl'
coldensdata = cPickle.load(open(modelname_HI,'rb'))
coldensdata = coldensdata[thisindex]
data['coldens_HI'] = coldensdata[idx[100]]

cPickle.dump(data,open(fileout,'wb'),protocol=-1)



