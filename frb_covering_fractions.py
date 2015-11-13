#Code to plot the frbs and radial profiles or my emission predictions
#Created by: Lauren
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import re
from lauren import find_covering_fraction

#########################################################################

#ions = ['SiIV','CIV','OVI']
#ions = ['CIII_977','CIV','OVI']

no_axes_labels = True
obs_colors = True

#ASTROFEST PARAMETERS
#model_gqs = ['g1q1']
#res_keys = ['1kpc']
#redshift_keys = ['z02']
#znow = [0.]
#ions = ['SiIV','CIII_977','CIV','OVI']
#xlen,ylen = 1,4
#figxlen,figylen = 3,9
#fileout = 'frb_emis_z02_obs_nozscaling.png'

##REDSHIFT EVOLUTION PARAMETERS
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'
temp_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'
model_gqs = ['g1q1','g1q1','g1q1','g1q1']
res_keys = ['1kpc','1kpc','1kpc','1kpc']
redshift_keys = ['z0','z02','z05','z1']
znow = [0.,0.2,0.5,1.0]
#znow = [0.,0.,0.,0.]
zplot = [0.,0.2,0.5,1.0]
ions = ['SiIV','CIII_977','CIV','OVI']	
xlen,ylen = 3,1
figxlen,figylen = 12,4
fileout = 'frb_covering_fraction_nodisk.png'#_nozscaling.png'
colors = ['Red','Magenta','Green','Blue']
SB_lims = [1.,2.,3.]

fig,ax = plt.subplots(ylen,xlen,sharey=True)
fig.set_size_inches(figxlen,figylen)
ax = ax.flat
i = 0

for ion in ions:
	print ion
	count = 0
	lim1 = []
	lim2 = []
	lim3 = []
	while count < len(model_gqs):
		print ion, count, model_gqs[count],i
		modelnames = [model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'.cpkl'] 
		tempname = temp_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_Temperature_Density.cpkl'

		#print tempname

		fractions = find_covering_fraction(modelnames[0],tempname,SB_lims,znow[count])
		lim1 = np.append(lim1,fractions[0])
		lim2 = np.append(lim2,fractions[1])
		lim3 = np.append(lim3,fractions[2])	
		count = count + 1

	ax[0].plot(zplot,lim1,'-o',color=colors[i],label=ion)
	ax[1].plot(zplot,lim2,'-o',color=colors[i],label=ion)
	ax[2].plot(zplot,lim3,'-o',color=colors[i],label=ion)
	i = i + 1

ax[0].set_title('SB > 1')
ax[1].set_title('SB > 2')
ax[2].set_title('SB > 3')
ax[0].set_ylabel('Fraction of Pixels')
ax[1].set_xlabel('Redshift')
ax[2].legend()
plt.tight_layout()
plt.savefig(fileout)#,transparent=True)
plt.close()



