#Code to assess my statistical model
#Created by: Lauren
#Created on: 12/4/2014
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

def basic_profile_vals(patt,pattend):
        rpr,rpmean,rpmedian,rpmax,rpmin = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)
        return rpr,rpmean,rpmedian,rpmax,rpmin

def make_radius_array():
        xL = np.linspace(-160,160,320)
        maxnow = 160.
        xL, yL = np.meshgrid(xL,xL)
        r = abs(xL+1j*yL)
        dr = np.abs([xL[0,0] - xL[0,1]])
        radial = np.arange(maxnow/dr)*dr + dr/2
        nrad = len(radial)
        return r, dr, nrad

def plot_scatter_points(ax,frbarr,r,dr,nrad):
        mslope = (0.0075-0.07)/nrad
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
        	alphahere = mslope*irad + 0.07
        	ax.plot(r[thisindex],np.log10(frbarr[thisindex]),'k.',alpha=alphahere)
        return

def full_scatter_plot(modelname,ion,ax,werk_data,max_r):
	r,dr,nrad = make_radius_array()
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	plot_scatter_points(ax,frbarr,r,dr,nrad)
	###FIX THIS###
	#rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)
	#idr = np.where(rpr <= max_r)[0]
	#ax.plot(rpr[idr],np.log10(rpmedian[idr]),'c-',linewidth=2.2)
	plot_Werk_ColDens(ax,werk_data,ion,'Rperp',xmax=max_r)
	ax.text(100,17.05,re.split('galquas/|/frb',modelname)[1])
	ax.set_xticks(range(0,160,30))
	ax.axis([0.,160.,8,18])
	return
#########################################################################

#ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
ions = ['SiIV','CIII','OVI']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
#removed g01q0.05 to compare and fit all on one line
#model_gqs = ['g01q01', 'g01q0.1','g01q0.5','g01q1','g01q2','g01q10','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
model_gqs = ['g1q01','g1q1','g1q10']
model_mid = '/frbz_1kpc_z02_'

werk_data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
l,u = 10.,20.
model_width = 0.0
nbins,ndraws = 500,10
max_r = 160.

for ion in ions:
	fileout = ion+'_evaluate.png'
	print fileout
	xlen,ylen = 6,6
	fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey='row')
	fig.set_size_inches(16,16)
	plt.subplots_adjust(.1,.1,.9,.9,0,0.1)
	ax = ax.flat
	i = 0
	yes,no,count = 0,0,0
	#Here, I'm looping through the subplots by sets of two rows.
	#The top row gets the scatter plot of the projected image data with the Werk data overplotted
	#The bottom plot is the Likelihood(Rperp) for each of the data points
	#The symbol of the data in the top and the bottom plots match to help guide the eye
	while i < len(ax):
		if count == len(model_gqs):
			if i+6 < xlen*ylen:
				ax[i+6].axis([0.,160.,-12,2])
			print 'done'	
		elif yes < xlen:
			print i, count,model_gqs[count]
			modelname = model_beg+model_gqs[count]+model_mid+ion+'dens.cpkl'
			full_scatter_plot(modelname,ion,ax[i],werk_data,max_r)
			
			lines = Line(ion,modelname,l,u,model_width,nbins,ndraws)
			lines.total_probability()	
			idx = np.where(lines.label=='u')[0]
			if len(idx) > 0:
				ax[i+6].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.4,uplims=True,fmt=None,ecolor='m',capsize=5,elinewidth=2,mew=0)

			idx = np.where(lines.label=='l')[0]
                        if len(idx) > 0:
                                ax[i+6].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.4,lolims=True,fmt=None,ecolor='Crimson',capsize=5,elinewidth=2,mew=0)

			idx = np.where(lines.label=='n')[0]
                        if len(idx) > 0:
                                ax[i+6].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.001,fmt='o',color='DarkOrange',ecolor='m')
			
			ax[i+6].axis([0.,160.,-12,2])
	
			yes = yes + 1
			count = count + 1
		else:
			no = no + 1
			if no == xlen:
				no = 0
				yes = 0
		i = i + 1

	plt.savefig(fileout)
	plt.close()



