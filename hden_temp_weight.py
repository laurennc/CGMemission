import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cPickle
#from yt.mods import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MultipleLocator
import sys


def hden_temp_hist(weight_basic,weight_col,z_key,HMkey):	
	hden_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'+z_key+'/'+HMkey+'/frbx_1kpc_500kpc_'+z_key+'_H_NumberDensity_'+weight_basic+'.cpkl'
	temp_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'+z_key+'/'+HMkey+'/frbx_1kpc_500kpc_'+z_key+'_Temperature_'+weight_basic+'.cpkl'
	column_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'+HMkey+'/frbx_1kpc_500kpc_z02_'+weight_col+'dens.cpkl'
	
	hden = cPickle.load(open(hden_file,'rb'))
	temp = cPickle.load(open(temp_file,'rb'))
	coldens = cPickle.load(open(column_file,'rb'))
	
	hden,temp,coldens = np.log10(hden.flat),np.log10(temp.flat),np.log10(coldens.flat)
	
	hist, xedges, yedges  = np.histogram2d(hden,temp,range=[[-6,1],[4, 6]],bins=100)#,weights=coldens)
	
	xbins = np.digitize(hden,xedges[1:-1])
	ybins = np.digitize(temp,yedges[1:-1])
	
	totvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))
	numvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))
	maxvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))
	
	for i in range(len(coldens)):
	    totvals[xbins[i],ybins[i]] += 10.**coldens[i]
	    numvals[xbins[i],ybins[i]] += 1
	    if 10.**coldens[i] > maxvals[xbins[i],ybins[i]]:
        	maxvals[xbins[i],ybins[i]] = 10.**coldens[i]

	avgvals = totvals/numvals
	
	idx = np.where(totvals == 0)
	avgvals[idx] = 0

	return hist,avgvals,maxvals


try:
	plot_key = sys.argv[1]
except:
	raise Exception('hey! you need a key so it knows what plots to make! use dev or paper!')

emis = cPickle.load(open('cloudywerk.cpkl','rb'))

weight_col = 'SiIII'
weight_basic = 'SiIII_Density'
redshift_key = 'z02'
HMkey = 'g1q1'


if plot_key == 'paper':
	fileout = 'hden_temp_paper_'+redshift_key+'_'+HMkey+'.png'
	
	weights_col = ['SiIV','CIII','OVI']
	weights_basic = ['SiIV_Density','CIII_Density','OVI_Density']

	fig,ax = plt.subplots(3,2,sharex=True,sharey=True)
	fig.set_size_inches(8,12)
	ax= ax.flat	
	i,j = 0,0

	while i < len(weights_col):
		hist,avgvals,maxvals = hden_temp_hist(weights_basic[i],weights_col[i],redshift_key,HMkey)
		im1 = ax[j].imshow(np.log10(hist.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',cmap='Blues')
		im2 = ax[j+1].imshow(np.log10(avgvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',vmin=9.,vmax=16.,cmap='jet')
		x0,x1 = ax[j].get_xlim()
        	y0,y1 = ax[j].get_ylim()

		ax[j].set_aspect((x1-x0)/(y1-y0))
		ax[j].set_xlim(-6,1)
           	ax[j].set_ylim(3.5,6)

		ax[j+1].set_aspect((x1-x0)/(y1-y0))
                ax[j+1].set_xlim(-6,1)
                ax[j+1].set_ylim(3.5,6)
                ax[j+1].set_title(weights_col[i])

		ax[j].plot(emis['hden'],emis['T'],'s',color='HotPink',alpha=0.75,ms=3.5)
		i = i + 1
		j = j + 2

	fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im2, cax=cbar_ax)#,shrink=0.6,ticks=MultipleLocator(0.5), format="%.2f")

        plt.savefig(fileout)
        plt.close()

if plot_key == 'dev':
	fileout = 'histdens_z02_g1q01_'+weight_basic+'_nohistweight.png'
	
	fig,ax = plt.subplots(1,3)
	fig.set_size_inches(12,4)
	ax = ax.flat

	hist,avgvals,maxvals = hden_temp_hist(weight_basic,weight_col,redshift_key,HMkey)	
	
	im1 = ax[0].imshow(np.log10(hist.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',cmap='Blues')
	im2 = ax[1].imshow(np.log10(avgvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower')#,vmin=9.,vmax=16.)
	im3 = ax[2].imshow(np.log10(maxvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower')#,vmin=9.,vmax=16.)
	
	x0,x1 = ax[2].get_xlim()
	y0,y1 = ax[2].get_ylim()
	titles = ['Number','Avg Column Density','Max Column Density']

	for j in range(len(ax)):
	   ax[j].set_aspect((x1-x0)/(y1-y0))
	   ax[j].set_xlim(-6,1)
	   ax[j].set_ylim(3.5,6)
	   ax[j].set_title(titles[j])
  
	ax[0].plot(emis['hden'],emis['T'],'s',color='HotPink',alpha=0.75,ms=3.5)

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(im3, cax=cbar_ax)#,shrink=0.6,ticks=MultipleLocator(0.5), format="%.2f")

	plt.savefig(fileout)
	plt.close()

