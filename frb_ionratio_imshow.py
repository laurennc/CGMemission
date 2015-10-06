#Code to plot the frbs and radial profiles or my emission predictions
#Created by: Lauren
#Created on: 12/11/2014
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from plotting_routines import *
from matplotlib import colors
import brewer2mpl as brew

############FUNCTIONS FOR THE PLOTTING####################################

def plot_frb(modelname,ax,z=0.,obs_colors=True):
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	if obs_colors:
		cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
		bounds = [-5,1,2,3,5]	
		norm = colors.BoundaryNorm(bounds,cmap.N)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',cmap=cmap,norm=norm,origin='lower')
	else:
		#bmap = brew.get_map('PuRd','Sequential',9)
		bmap = brew.get_map('PRGn','Diverging',9)
		cmap = bmap.get_mpl_colormap(N=1000, gamma=2.0)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',cmap=cmap,origin='lower')

	return im


#########################################################################

no_axes_labels = True
obs_colors = True

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/' ##CHANGED FORM Z02

#ASTROFEST PARAMETERS
model_gqs = ['g1q10']
res_keys = ['1kpc']
redshift_keys = ['z02']
znow = [0.]
ions = ['SiIV','CIII_977','CIV','OVI']
xlen,ylen = 1,3
figxlen,figylen = 3,9
fileout = 'frb_ionfracs_z02_nozscaling.png'


max_r = 160.
percentile = 0.01
ncontours = 4

#fig,ax = plt.subplots(ylen,xlen,sharey=True)
#fig.set_size_inches(figxlen,figylen)
#ax = ax.flat

count = 0

frb_SiIV = np.array(cPickle.load(open(model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_SiIV.cpkl','rb')))

frb_CIII = np.array(cPickle.load(open(model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_CIII_977.cpkl','rb')))

frb_CIV = np.array(cPickle.load(open(model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_CIV.cpkl','rb')))

frb_OVI = np.array(cPickle.load(open(model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_OVI.cpkl','rb')))


#bmap = brew.get_map('PuRd','Sequential',9)
bmap = brew.get_map('PRGn','Diverging',9)
cmap = bmap.get_mpl_colormap(N=1000, gamma=2.0)



plt.imshow(np.log10(frb_SiIV/frb_CIII),extent=(-160,160,160,-160),interpolation='none',cmap=cmap,origin='lower')
plt.colorbar()
plt.savefig('SiIV_CIII_'+redshift_keys[count]+'_'+model_gqs[count]+'.png')
plt.close()

plt.imshow(np.log10(frb_CIII/frb_CIV),extent=(-160,160,160,-160),interpolation='none',cmap=cmap,origin='lower',vmax=1.5,vmin=-2.5)
plt.colorbar()
plt.savefig('CIII_CIV'+redshift_keys[count]+'_'+model_gqs[count]+'.png')
plt.close()

plt.imshow(np.log10(frb_CIII/frb_OVI),extent=(-160,160,160,-160),interpolation='none',cmap=cmap,origin='lower',vmax=1.5,vmin=-4.0)
plt.colorbar()
plt.savefig('CIII_OVI'+redshift_keys[count]+'_'+model_gqs[count]+'.png')
plt.close()



