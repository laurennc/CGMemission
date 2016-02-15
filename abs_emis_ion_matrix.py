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
import matplotlib.gridspec as gridspec
#from Line_ORG import *
#from Galaxy import *
import re
from lauren_holding import make_SB_profile
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

def plot_scatter_points(ax,frbarr,r,dr,nrad):#,rpmean):
        mslope = (0.015-0.07)/nrad
	total_pts, num_above = 0.,0.
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
        	#For the percentage calculation
		total_pts = total_pts + len(frbarr[thisindex])
		alphahere = mslope*irad + 0.07
        	ax.plot(r[thisindex],np.log10(frbarr[thisindex]),'.',alpha=alphahere,color='Gray')
        return #percent_above

def full_scatter_plot(modelname_emis,profile_names_emis,ion1,modelname_abs,profile_names_abs,ion2,ax,werk_data,max_r):
	r,dr,nrad = make_radius_array()
	frbarr_emis = np.array(cPickle.load(open(modelname_emis,'rb')))
	frbarr_abs = np.array(cPickle.load(open(modelname_abs,'rb')))
	xL = np.linspace(-160,160,320)
        #rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(profile_names[0],profile_names[1],profile_names[2],xL)
	frbarr = frbarr_emis/(frbarr_abs**2.0)
	percent_above = plot_scatter_points(ax,frbarr,r,dr,nrad)#,rpmean)
	#idr = np.where(rpr <= max_r)[0]
	#ax.plot(rpr[idr],np.log10(rpmedian[idr]),'-',linewidth=1.7,color='Black')
	ax.set_xticks(range(0,160,30))
	#ax.axis([0.,160.,8,18])
	return percent_above
#########################################################################

ions1 = ['CIII_977','SiIV','CIV']
ions2 = ['CIII','SiIV','CIV']
model_beg_emis = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/z02/'
model_beg_abs = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/coldens/z02/'

model_gqs = ['g1q01','g1q1','g1q10']
model_mid = '/frbx_1kpc_500kpc_z02_'

werk_data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
l,u = 10.,20.
model_width = 0.0
nbins,ndraws = 500,10
max_r = 160.

fileout = 'emis_abs_scatter_matrix_z02_nozscale.png'
xlen,ylen = 3,3 #6#2,6
fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
#fig.set_size_inches(12,4)
fig.set_size_inches(8,8)#,16)
ax = ax.flat
i = 0
r,dr,nrad = make_radius_array()
xL = np.linspace(-160,160,320)


model_gqs = ['g1q01','g1q1','g1q10']
i2 = -1

for ion in ions1:
	count = 0
	i2 = i2 + 1
	while count < 3:
		print ion, count, model_gqs[count],i
		ion2 = ions2[i2]

		modelname_emis = model_beg_emis+model_gqs[count]+model_mid+ion+'.cpkl'
                profile_names_emis = [model_beg_emis+model_gqs[count]+model_mid+ion+'.cpkl',model_beg_emis+model_gqs[count]+model_mid+ion+'.cpkl',model_beg_emis+model_gqs[count]+model_mid+ion+'.cpkl']

		modelname_abs = model_beg_abs+model_gqs[count]+model_mid+ion2+'dens.cpkl'
                profile_names_abs = [model_beg_abs+model_gqs[count]+model_mid+ion2+'dens.cpkl',model_beg_abs+model_gqs[count]+model_mid+ion2+'dens.cpkl',model_beg_abs+model_gqs[count]+model_mid+ion2+'dense.cpkl']


		full_scatter_plot(modelname_emis,profile_names_emis,ion,modelname_abs,profile_names_abs,ion2,ax[i],werk_data,max_r)

		
		count = count + 1
		i = i + 1	


for j in range(len(model_gqs)):
	ax[j].set_title(model_gqs[j])

for k in range(len(ions1)):
	ax[k*xlen+2].text(100,17.05,ions1[k]+'/'+ions2[k])

ax[7].set_xlabel('Impact Parameter [kpc]')
ax[3].set_ylabel('Emission / (Column Density)^2')


plt.tight_layout()
plt.savefig(fileout)#,transparent=True)
plt.close()
f.close()



