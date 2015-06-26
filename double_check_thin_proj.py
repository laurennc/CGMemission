import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import itertools

###OK WHAT I WANT THIS TO DO IS PLOT THE FRBS SIDE BY SIDE WITH THE SAME COLOR BAR SO I CAN GET BETTER IDEA AS TO WHATS HAPPENING

fileout = 'frb_comparison_g1q1_Zfixed.png'
xlen,ylen = 3,4
fig,ax = plt.subplots(ylen,xlen)
ax = ax.flat
fig.set_size_inches(12,16)#,16)

i = 0

model_gqs = ['g1q1']
ions = ['HI','SiIV','CIII','OVI']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
model_mid3 = '/frbz_1kpc_320kpc_z02_'
model_mid2 = '/frbz_1kpc_500kpc_z02_'
model_mid1 = '/frbz_1kpc_z02_'
ions = ['HI','SiIV','CIII','OVI']

i = 0
for ion in ions:
	m1 = model_beg+model_gqs[0]+model_mid1+ion+'dens.cpkl'
	frb1 = cPickle.load(open(m1,'rb'))
	frb1 = np.log10(frb1)
	m2 = model_beg+model_gqs[0]+model_mid2+ion+'dens.cpkl'
	frb2 = cPickle.load(open(m2,'rb'))
	frb2 = np.log10(frb2)
	m3 = model_beg+model_gqs[0]+model_mid3+ion+'dens.cpkl'
	frb3 = cPickle.load(open(m3,'rb'))
	frb3 = np.log10(frb3)

	im = ax[i].imshow(frb1,vmin=10.,vmax=18.)
	ax[i].set_axis_off()
	ax[i+1].imshow(frb2,vmin=10.,vmax=18.)
	ax[i+1].set_axis_off()
	ax[i+2].imshow(frb3,vmin=10.,vmax=18.)
	ax[i+2].set_axis_off()
	
	if i == 0:
		ax[i].set_title('full proj')
		ax[i+1].set_title('500kpc ')
		ax[i+2].set_title('320kpc '+ion)
	else:
		ax[i+2].set_title(ion)

	i = i + 3

#fig.subplots_adjust(right=0.75)
#cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
#fig.colorbar(im, cax=cbar_ax)

plt.tight_layout()
plt.savefig(fileout)
plt.close()




