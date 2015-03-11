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

def plot_total_fractions(ion):
	count = 0
        fileout = ion+'_pixelFractions.png'
        fractions = np.zeros((len(model_gqs),4))
        while count < len(model_gqs):

                modelname = model_beg+model_gqs[count]+'/frbx_'+res_keys[count]+redshifts[count]+ion+'.cpkl'
                frbarr = np.array(cPickle.load(open(modelname,'rb')))

                totnum = float(len(frbarr.flatten()))
                lims = [10,3,2,1,-10]
                colors = ['Chartreuse','DarkTurquoise','HotPink','Gray']
                i = 0
                while i < len(lims)-1:
                        wanted = np.where((frbarr < lims[i]) & (frbarr >= lims[i+1]))[0]
                        fractions[count,i] = float(len(wanted.flatten()))/totnum
                        i = i + 1

                count = count + 1

        x = np.array(range(len(model_gqs)))+0.5

        for j in range(3):
                plt.plot(x,fractions[:,j],'o',color=colors[j],label=labels[j])

        plt.axvline(x=3,linestyle='--',color='DimGray',linewidth=1.5)
        plt.axvline(x=6,linestyle='--',color='DimGray',linewidth=1.5)
        plt.xticks(x, ticklabels, rotation=45)
        plt.ylabel('Fraction of Pixels')
        plt.legend(loc=1,scatterpoints=1)
        #plt.ylim(0.,1.)
        plt.tight_layout()
        plt.savefig(fileout)
        plt.close()

	return

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

def plot_fraction_vs_radius(ax,ion,frbarr,res_key,title_now):
	xL,r,dr,nrad = make_radius_array(res_key)
	fractions = np.zeros((nrad,4))	
	rplot = []
	for irad in range(nrad):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (r>=minrad) * (r<maxrad)
		frbnow = np.log10(frbarr[thisindex])	
		rplot = np.append(rplot,minrad+(dr/2.))

		lims = [10,3,2,1,-10]
		i = 0
                while i < len(lims)-1:
                        wanted = np.where((frbnow < lims[i]) & (frbnow >= lims[i+1]))[0]
                        fractions[irad,i] = float(len(wanted.flatten()))/float(len(frbnow))
                        i = i + 1
	colors = ['Chartreuse','DarkTurquoise','HotPink','Gray']
	labels = ['x>3','2<x<3','1<x<2','x<1']
	for j in range(4):
		#print 'len(r): ',len(rplot)
		#print 'len(fractions): ',len(fractions[:,j])
		ax.plot(rplot,fractions[:,j],'-',linewidth=2.0,color=colors[j],label=labels[j])
		ax.set_xticks(range(0,160,30))
		ax.set_title(title_now)

	return	

#ions = ['SiIV','CIII','OVI']
ions = ['HAlpha','CIII_977','CIV','MgII','SiII','SiIII_1207','SiIII_1883']
ions = ['OVI']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/z02/'
model_gqs = ['g1q01','g1q01','g1q01','g1q1','g1q1','g1q1','g1q10','g1q10','g1q10']
res_keys = ['1kpc','5kpc','25kpc','1kpc','5kpc','25kpc','1kpc','5kpc','25kpc']
redshifts = ['_z02_','_z02_','_z02_','_z02_','_z02_','_z02_','_z02_','_z02_','_z02_']

max_r,percentile = 160.,0.01
labels = ['x>3','2<x<3','1<x<2','x<1']
ticklabels = ['g1q01 (1kpc)','g1q01 (5kpc)','g1q01 (25kpc)','g1q1 (1kpc)','g1q1 (5kpc)','g1q1 (25kpc)','g1q10 (1kpc)','g1q10 (5kpc)','g1q10 (25kpc)']

for ion in ions:
	fileout = ion+'_pixel_fraction_radius.png'
	xlen,ylen = 3,3
        fig,ax = plt.subplots(ylen,xlen,sharey=True,sharex=True)
        fig.set_size_inches(8,8)
        plt.subplots_adjust(.1,.1,.9,.9,0,0.1)
        ax = ax.flat
        i = 0
	
	while i < len(ax):
		modelname = model_beg+model_gqs[i]+'/frbx_'+res_keys[i]+redshifts[i]+ion+'.cpkl'
                frbarr = np.array(cPickle.load(open(modelname,'rb')))

		plot_fraction_vs_radius(ax[i],ion,frbarr,res_keys[i],model_gqs[i])
		i = i + 1

	plt.tight_layout()
	plt.savefig(fileout)
	plt.close()

	#plot_total_fractions(ion)


