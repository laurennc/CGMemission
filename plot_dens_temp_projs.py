#Code to plot the frb projections of density and temperature for paper figure
#Created by: Lauren
#Created on: 10/26/2014
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import brewer2mpl as brew

no_axes_labels = True
model_gqs = ['g1q1','g1q1','g1q1','g1q1']
res_keys = ['1kpc','1kpc','1kpc','1kpc']
redshift_keys = ['z0','z02','z05','z1']
xlen,ylen = 2,4
figxlen,figylen = 4,8
fileout = 'dens_temp_projs.pdf'
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'

fig,ax = plt.subplots(ylen,xlen,sharey=True,sharex=True)
fig.set_size_inches(figxlen,figylen)
ax = ax.flat
count = 0
i = 0

bmap = brew.get_map('PuBu','Sequential',9)
dens_cmap = bmap.get_mpl_colormap(N=1000,gamma=2.0)
bmap = brew.get_map('OrRd','Sequential',9)
temp_cmap = bmap.get_mpl_colormap(N=1000,gamma=2.0)

dens_cmap = 'gnuplot2'
#temp_cmap = 'jet'
plt.axis('on')

while count < len(res_keys):
	tempname = model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_Temperature_Densityweight.cpkl'
	densname = model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_Density.cpkl'

	temp_frb = np.log10(cPickle.load(open(tempname,'rb')))
	dens_frb = np.log10(cPickle.load(open(densname,'rb')))

	###PLOT THE DENSITY FIRST### 
	im = ax[i].imshow(dens_frb,extent=(-160,160,160,-160),interpolation='none',cmap=dens_cmap,origin='lower',vmin=-5,vmax=-1)
	###PLOT TEMPERATURE###
	im2 = ax[i+1].imshow(temp_frb,extent=(-160,160,160,-160),interpolation='none',cmap=temp_cmap,origin='lower',vmin=4,vmax=6)

	###TURN AXIS LABELS OFF###
	ax[i].set_xticklabels([])
	ax[i+1].set_xticklabels([])
	ax[i].set_yticklabels([])
	ax[i+1].set_yticklabels([])

	i = i + 2
	count = count + 1

##rect = l,b,w,h
fig.subplots_adjust(bottom=0.125)
cbar_ax = fig.add_axes([0.125,0.08,0.353,0.03])
fig.colorbar(im,cax=cbar_ax,orientation='horizontal',ticks=[-4, -3, -2,-1])
cbar_ax = fig.add_axes([0.548,0.08,0.353,0.030])
fig.colorbar(im2,cax=cbar_ax,orientation='horizontal',ticks=[4,5,6])

plt.savefig(fileout)
plt.close()


