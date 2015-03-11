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

def plot_scatter_points(ax,frbarr,r,dr,nrad,rpmean):
        mslope = (0.0075-0.07)/nrad
	total_pts, num_above = 0.,0.
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
        	#For the percentage calculation
		total_pts = total_pts + len(frbarr[thisindex])
		num_above = num_above + len(np.where(np.log10(frbarr[thisindex])>np.log10(rpmean[irad]))[0])
		#FOr the scatter plot
		alphahere = mslope*irad + 0.07
        	ax.plot(r[thisindex],np.log10(frbarr[thisindex]),'k.',alpha=alphahere)
	#print float(num_above),float(total_pts)
	percent_above = float(num_above)/float(total_pts)
        return percent_above

def full_scatter_plot(modelname,profile_names,ion,ax,werk_data,max_r):
	r,dr,nrad = make_radius_array()
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	###FIX THIS###
	xL = np.linspace(-160,160,320)
        rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(profile_names[0],profile_names[1],profile_names[2],xL)
	percent_above = plot_scatter_points(ax,frbarr,r,dr,nrad,rpmean)
	idr = np.where(rpr <= max_r)[0]
	print rpr[50]
	ax.plot(rpr[idr],np.log10(rpmedian[idr]),'c-',linewidth=2.2)
	plot_Werk_ColDens(ax,werk_data,ion,'Rperp',xmax=max_r)
	#ax.text(100,17.05,re.split('galquas/|/frb',modelname)[1])
	ax.set_xticks(range(0,160,30))
	ax.axis([0.,160.,8,18])
	return percent_above
#########################################################################

#ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
ions = ['SiIV','CIII','OVI']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
#removed g01q0.05 to compare and fit all on one line
#model_gqs = ['g01q01', 'g01q0.1','g01q0.5','g01q1','g01q2','g01q10','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
model_gqs = ['g1q01','g1q1','g1q10']#,'g1q01','g1q1','g1q10','g1q01','g1q1','g1q10']
model_mid = '/frbz_1kpc_z02_'
#model_mid = '/frbFace_default_z02_'

werk_data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
l,u = 10.,20.
model_width = 0.0
nbins,ndraws = 500,10
max_r = 160.

#fileout = 'scatter_and_likelihood_HIGHhden.png'
fileout = 'scatter_matrix.png'
xlen,ylen = 3,3 #6#2,6
fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey='row')
fig.set_size_inches(8,8)#,16)
plt.subplots_adjust(.1,.1,.9,.9,0,0.1)
ax = ax.flat
i = 0

f = open('percent_above.dat', 'w')

model_gqs = ['g1q01','g1q1','g1q10']


for ion in ions:
	count = 0
	#while count < 2:
	while count < 3:
		print ion, count, model_gqs[count],i

	#	if count > 0:
	#		modelname = model_beg+model_gqs[count]+model_mid+ion+'_dens_HIGHhden.cpkl'
	#		profile_names = [model_beg+model_gqs[count]+'/frbx_1kpc_z02_'+ion+'_dens_HIGHhden.cpkl',model_beg+model_gqs[count]+'/frby_1kpc_z02_'+ion+'_dens_HIGHhden.cpkl',model_beg+model_gqs[count]+'/frbz_1kpc_z02_'+ion+'_dens_HIGHhden.cpkl']
	#	else:
		modelname = model_beg+model_gqs[count]+model_mid+ion+'dens.cpkl'
                profile_names = [model_beg+model_gqs[count]+'/frbx_1kpc_z02_'+ion+'dens.cpkl',model_beg+model_gqs[count]+'/frby_1kpc_z02_'+ion+'dens.cpkl',model_beg+model_gqs[count]+'/frbz_1kpc_z02_'+ion+'dens.cpkl']

		percent_above = full_scatter_plot(modelname,profile_names,ion,ax[i],werk_data,max_r)
		
	#	lines = Line(ion,modelname,l,u,model_width,nbins,ndraws)
	#	lines.total_probability()	
	#	idx = np.where(lines.label=='u')[0]
	#	if len(idx) > 0:
	#		ax[i+xlen].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.4,uplims=True,fmt=None,ecolor='m',capsize=5,elinewidth=2,mew=0)

	#	idx = np.where(lines.label=='l')[0]
        #	if len(idx) > 0:
        #	        ax[i+xlen].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.4,lolims=True,fmt=None,ecolor='Crimson',capsize=5,elinewidth=2,mew=0)
	
	#	idx = np.where(lines.label=='n')[0]
        #	if len(idx) > 0:
        #	        ax[i+xlen].errorbar(lines.rperp[idx],lines.likelihoods[idx],yerr=0.001,fmt='o',color='DarkOrange',ecolor='m')
			
	#	ax[i+xlen].axis([0.,160.,-12,2])

	#	printout = modelname+'\t'+str(percent_above)+'\n'
	#	f.write(printout)	

		count = count + 1
		i = i + 1	

	#i = i  + xlen

for j in range(len(model_gqs)):
	ax[j].set_title(model_gqs[j])

for k in range(len(ions)):
	print k*xlen+2, ions[k]
	ax[k*xlen+2].text(100,17.05,ions[k])


plt.tight_layout()
plt.savefig(fileout)
plt.close()
f.close()



