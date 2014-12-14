#Code to plot the frbs and radial profiles or my emission predictions
#Created by: Lauren
#Created on: 12/11/2014
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import itertools
from Line_ORG import *
#from Galaxy import *
import re
from lauren import make_SB_profile
from radial_data_lauren import *
from plotting_routines import *


############FUNCTIONS FOR THE PLOTTING####################################

def basic_profile_vals(modelnames,xL):
        rpr,rpmean,rpmedian,rpmax,rpmin,rp_std = make_SB_profile(modelnames[0],modelnames[1],modelnames[2],xL)
        return rpr,rpmean,rpmedian,rpmax,rpmin

def make_radius_array(key):
	if key=='1kpc':
        	xL = np.linspace(-160,160,320)
        elif key=='5kpc':
		xL = np.linspace(-160,160,64)
	elif key=='25kpc':
		xL = np.linspace(-160,160,13)
	else:
		print 'Failure in the key'
	maxnow = 160.
        xL2, yL = np.meshgrid(xL,xL)
        r = abs(xL2+1j*yL)
        dr = np.abs([xL2[0,0] - xL2[0,1]])
        radial = np.arange(maxnow/dr)*dr + dr/2
        nrad = len(radial)
        return xL,r, dr, nrad

def plot_scatter_points_percentile(ax,frbarr,r,dr,nrad,percentile,key):
        mslope = (0.0075-0.07)/nrad
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
        	alphahere = mslope*irad + 0.07
		if key=='1kpc':
        		ax.plot(r[thisindex],np.log10(frbarr[thisindex]),'k.',alpha=alphahere)
		else:
			ax.plot(r[thisindex],np.log10(frbarr[thisindex]),'k.',alpha=0.37)
		rhere = r[thisindex].flatten()
                datahere = frbarr[thisindex].flatten()
                meanhere = frbarr[thisindex].mean()
                lenperc = int(len(datahere)*percentile)
                idSort = np.argsort(datahere)[::-1]
                wanted = idSort[0:lenperc]
                if len(wanted) == 0:
                        wanted = np.where(datahere == datahere.max())
                ax.plot(rhere[wanted],np.log10(datahere[wanted]),'g.',alpha=0.5)
        return

def full_scatter_plot(modelnames,ion,ax,res_key,max_r,percentile):
	xL,r,dr,nrad = make_radius_array(res_key)
	frbarr = np.array(cPickle.load(open(modelnames[0],'rb')))
	plot_scatter_points_percentile(ax,frbarr,r,dr,nrad,percentile,res_key)
	rpr,rpmean,rpmedian,rpmax,rpmin = basic_profile_vals(modelnames,xL)
	
	idr = np.where(rpr <= 160.0)
	ax.plot(rpr[idr],np.log10(rpmedian[idr]),color='Goldenrod',linewidth=2.2)	

	ax.set_xticks(range(0,160,30))
	ax.set_xlim(0,160)
	#ax.axis([0.,160.,-2.,5.])
	return

def plot_frb(modelname,ax,include_colorbar=False):
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	im = ax.imshow(np.log10(frbarr),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none')
	#if include_colorbar==True:
	#	plt.colorbar(im)
	return
#########################################################################

#ions = ['HI','MgII','SiII','SiIII']#,'SiIV','CIII','OVI']
ions = ['SiIV','CIII','OVI']
#ions = ['HAlpha','CIII_977','CIV','MgII','SiII','SiIII_1207','SiIII_1883']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/'
model_gqs = ['g1q01','g1q1','g1q10','g1q01','g1q1','g1q10','g1q01','g1q1','g1q10']
model_mid = '/frbz_1kpc_z02_'
res_keys = ['1kpc','1kpc','1kpc','5kpc','5kpc','5kpc','25kpc','25kpc','25kpc']

werk_data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
l,u = 10.,20.
model_width = 0.0
nbins,ndraws = 500,10
max_r,percentile = 160.,0.01


for ion in ions:
	fileout = ion+'_frb_profile_matrix_nointerp.png'
	xlen,ylen = 3,6
	first = [k*xlen for k in range(ylen)]
	last = np.array(first)+xlen-1
	fig,ax = plt.subplots(ylen,xlen,sharey='row')
	fig.set_size_inches(8,16)
	plt.subplots_adjust(.1,.1,.9,.9,0,0.1)
	ax = ax.flat
	i = 0
	yes,no,count = 0,0,0

	while i < len(ax):
		if yes < xlen:
			modelnames = [model_beg+model_gqs[count]+'/frbx_'+res_keys[count]+'_z02_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frby_'+res_keys[count]+'_z02_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frbz_'+res_keys[count]+'_z02_'+ion+'.cpkl'] 

			if i in last:
				plot_frb(modelnames[0],ax[i],include_colorbar=True)
			else:
				plot_frb(modelnames[0],ax[i])

			full_scatter_plot(modelnames,ion,ax[i+xlen],res_keys[count],max_r,percentile)
			yes  = yes + 1
			count = count + 1
		else:
			no = no + 1
			if no == xlen:
				no = 0
				yes = 0
		i = i + 1	


	for j in range(xlen):
		ax[j].set_title(model_gqs[j])
	
	k = 1 
	while k  < len(first):
		ax[first[k]].set_ylabel('Emission [photon s^-1 cm^-2 sr^-1]')
		k = k + 2

	
	plt.tight_layout()
	plt.savefig(fileout)
	plt.close()



