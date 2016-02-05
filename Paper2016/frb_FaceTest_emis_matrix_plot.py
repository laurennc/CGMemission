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
from matplotlib import colors

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
		rnow,frbnow = r[thisindex],np.log10(frbarr[thisindex])	
	
		lims = [10,3,2,1,-10]
                colors = ['Chartreuse','DarkTurquoise','HotPink','Gray']
                i = 0
		while i < len(lims)-1:
			wanted = np.where((frbnow < lims[i]) & (frbnow >= lims[i+1]))[0]
			if key=='1kpc':
				ax.plot(rnow[wanted],frbnow[wanted],'.',color=colors[i],alpha=alphahere)
			else:
				ax.plot(rnow[wanted],frbnow[wanted],'.',color=colors[i],alpha=0.37)
			i = i + 1

		#rhere = r[thisindex].flatten()
                #datahere = frbarr[thisindex].flatten()
                #meanhere = frbarr[thisindex].mean()
                #lenperc = int(len(datahere)*percentile)
                #idSort = np.argsort(datahere)[::-1]
                #wanted = idSort[0:lenperc]
                #if len(wanted) == 0:
                #        wanted = np.where(datahere == datahere.max())
                #ax.plot(rhere[wanted],np.log10(datahere[wanted]),'g.',alpha=0.5)
        return

def full_scatter_plot(modelnames,ion,ax,res_key,max_r,percentile):
	xL,r,dr,nrad = make_radius_array(res_key)
	frbarr = np.array(cPickle.load(open(modelnames[0],'rb')))
	plot_scatter_points_percentile(ax,frbarr,r,dr,nrad,percentile,res_key)
	rpr,rpmean,rpmedian,rpmax,rpmin = basic_profile_vals(modelnames,xL)
	
	idr = np.where(rpr <= 160.0)
	ax.plot(rpr[idr],np.log10(rpmedian[idr]),color='k',linewidth=2.2)#previously Goldenrod	

	ax.set_xticks(range(0,160,30))
	#ax.set_xlim(0,160)
	ax.axis([0.,160.,-6.,6.5])
	return

def plot_frb(modelname,ax,include_colorbar=False):
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
	bounds = [-5,1,2,3,5]
	norm = colors.BoundaryNorm(bounds,cmap.N)
	im = ax.imshow(np.log10(frbarr),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',cmap=cmap,norm=norm)
	#if include_colorbar==True:
	#	plt.colorbar(im)
	return

def add_HI_contour(ax,HIfrb,ncontours):
	xL,r,dr,nrad = make_radius_array('1kpc')
	X,Y = np.meshgrid(xL,xL)
	CS = ax.contour(X,Y,HIfrb,[15,18],color='k')
	#plt.clabel(CS, fontsize=9, inline=1)
	return
#########################################################################

#ions = ['HI','MgII','SiII','SiIII']#,'SiIV','CIII','OVI']
ions = ['SiIV','CIII','OVI']
#ions = ['HAlpha','CIII_977','CIV','MgII','SiII','SiIII_1207','SiIII_1883']
redshift_key = 'z0'

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/'+redshift_key'/' ##CHANGED FORM Z02
HI_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
model_gqs = ['g1q01','g1q1','g1q10','g1q01','g1q1','g1q10','g1q01','g1q1','g1q10']
#model_mid = '/frbz_1kpc_z02_'
res_keys = ['1kpc','1kpc','1kpc','5kpc','5kpc','5kpc','25kpc','25kpc','25kpc']

max_r = 160.
percentile = 0.01

ncontours = 4

for ion in ions:
	fileout = ion+'_frbFace_profile_matrix_nointerp.png'
	xlen,ylen = 3,2
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
			modelnames = [model_beg+model_gqs[count]+'/frbx_'+res_keys[count]+'_'+redshift_key+'_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frby_'+res_keys[count]+'_'+redshift_key+'_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frbz_'+res_keys[count]+'_'+redshift_key+'_'+ion+'.cpkl'] 

			HIfile = HI_beg+model_gqs[count]+'/frbx_1kpc_'+redshift_key+'_HIdens.cpkl'
			HIfrb = cPickle.load(open(HIfile,'rb'))
			HIfrb = np.log10(HIfrb)

			if i in last:
				modelFace = model_beg+model_gqs[count]+'/frbFace_default_z0_'+ion+'.cpkl'
				plot_frb(modelnamesFace,ax[i],include_colorbar=True)
				add_HI_contour(ax[i],HIfrb,ncontours)
			else:
				modelFace = model_beg+model_gqs[count]+'/frbFace_default_z0_'+ion+'.cpkl'
				plot_frb(modelFace,ax[i])
				add_HI_contour(ax[i],HIfrb,ncontours)

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
	plt.savefig(fileout, dpi=1000)
	plt.close()



