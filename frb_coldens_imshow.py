#Code to plot the frbs and radial profiles or my emission predictions
#Created by: Lauren
#Created on: 12/11/2014
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.integrate as integrate
import itertools
import matplotlib.gridspec as gridspec
import re
from lauren_holding import make_SB_profile
from radial_data_lauren import *
from plotting_routines import *
from matplotlib import colors
import brewer2mpl as brew

############FUNCTIONS FOR THE PLOTTING####################################

def plot_frb(modelname,ax,vmin=12,vmax=16,z=0.,include_colorbar=False,obs_colors=True):
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	if obs_colors:
		cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
		bounds = [-5,1,2,3,5]	
		norm = colors.BoundaryNorm(bounds,cmap.N)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=vmin,vmax=max,interpolation='none',cmap=cmap,norm=norm,origin='lower')
	else:
		#bmap = brew.get_map('GnBu','Sequential',9)
		#cmap = bmap.get_mpl_colormap(N=1000, gamma=2.0)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=vmin,vmax=vmax,interpolation='none',cmap='YlGnBu',origin='lower')

	return im

#########################################################################

no_axes_labels = True
obs_colors = False

#model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/' ##CHANGED FORM Z02
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'

#ASTROFEST PARAMETERS
model_gqs = ['g1q1']
res_keys = ['1kpc']
redshift_keys = ['z02']
znow = [0.]
ions = ['HI','SiIV','CIII','OVI']
xlen,ylen = 1,4
figxlen,figylen = 3,9
#fileout = 'frb_coldens_z02_imshow.png'
fileout = 'frb_coldens_z02_imshow.pdf'

##REDSHIFT EVOLUTION PARAMETERS
#model_gqs = ['g1q1','g1q1','g1q1']
#res_keys = ['1kpc','1kpc','1kpc']
#redshift_keys = ['z0','z05','z1']
#znow = [0.,0.5,1.0]
#znow = [0.,0.,0.]
#ions = ['CIII_977','CIV','OVI']	
#xlen,ylen = 3,3
#figxlen,figylen = 12,12
#fileout = 'frb_emis_theory_nozscaling.png'

max_r = 160.
percentile = 0.01
ncontours = 4

fig,ax = plt.subplots(ylen,xlen,sharey=True)
fig.set_size_inches(figxlen,figylen)
#fig.set_size_inches(24,6)
#gs1 = gridspec.GridSpec(1, 4)
#gs1.update(wspace=0.025, hspace=0.05)
ax = ax.flat
i = 0
print fileout

for ion in ions:
	print ion
	count = 0
	while count < len(model_gqs):
		print ion, count, model_gqs[count],i
		modelnames = [model_beg+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'dens.cpkl']

		if ion=='HI':
			im_out_HI = plot_frb(modelnames[0],ax[i],vmin=12,vmax=24,z=znow[count],obs_colors=obs_colors)
			#ax[i].colorbar(im_out_HI)
		else:	
			im_out = plot_frb(modelnames[0],ax[i],vmin=12,vmax=16,z=znow[count],obs_colors=obs_colors)

		if no_axes_labels:
			#ax[i].set_title(model_gqs[count])
			plt.axis('on')
			ax[i].set_xticklabels([])
			ax[i].set_yticklabels([])
			ax[i].set_yticklabels([])

		ax[i].set_adjustable('box-forced')		
		#ax[i].set_title(ion)		
		ax[i].text(75,-100,ion)
		i = i + 1
		count = count + 1	

fig.subplots_adjust(right=0.85)
#divider1 = make_axes_locatable(ax[0])
#cax1 = divider1.append_axes("right",size="10%",pad=0.25)
cax1 = fig.add_axes([0.85,0.725,0.05,0.17])
cbar = plt.colorbar(im_out_HI,cax=cax1,ticks=[12,16,20,24])

##                     left, bottom, width, height
cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.59])
fig.colorbar(im_out, cax=cbar_ax,ticks=[12,13,14,15,16])

#for j in range(len(model_gqs)):
        #ax[j].set_title(model_gqs[j])
#	ax[j].set_title(redshift_keys[j])
	
#for k in range(len(ions)):
#        ax[k*xlen+2].text(75,-100,ions[k])

#plt.tight_layout()
plt.savefig(fileout)#,transparent=True)
plt.close()



