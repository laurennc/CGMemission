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

def plot_scatter_points(ax,frbarr1,frbarr2,temperature,r,dr,nrad):#,rpmean):
        mslope = (0.015-0.07)/nrad
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
		frbarr1here,frbarr2here,temphere,rhere = frbarr1[thisindex],frbarr2[thisindex],temperature[thisindex],r[thisindex]
		idx = np.where((frbarr1here > 1.0) & (frbarr2here > 1.0))
		alphahere = mslope*irad + 0.07
        	frbarrplot = np.log10(frbarr1here[idx]/frbarr2here[idx])
		ax.scatter(rhere[idx],frbarrplot,c=np.log10(temphere[idx]),alpha=alphahere)
	#print float(num_above),float(total_pts)
	#percent_above = float(num_above)/float(total_pts)
        return #percent_above

def full_scatter_plot(modelname1,profile_names1,ion1,modelname2,profile_names2,ion2,ax,werk_data,max_r):
	r,dr,nrad = make_radius_array()
	frbarr1 = np.array(cPickle.load(open(modelname1,'rb')))
	frbarr2 = np.array(cPickle.load(open(modelname2,'rb')))
	temperature = np.array(cPickle.load(open('/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/z02/g1q1//frbx_1kpc_500kpc_z02_Temperature_Density.cpkl','rb')))
	density = np.array(cPickle.load(open('/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/z02/g1q1//frbx_1kpc_500kpc_z02_Density.cpkl','rb')))

	###FIX THIS###
	xL = np.linspace(-160,160,320)
        #rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(profile_names[0],profile_names[1],profile_names[2],xL)
	#frbarr = frbarr1/frbarr2
	percent_above = plot_scatter_points(ax,frbarr1,frbarr2,density,r,dr,nrad)#,rpmean)
	#idr = np.where(rpr <= max_r)[0]
	#ax.plot(rpr[idr],np.log10(rpmedian[idr]),'-',linewidth=1.7,color='Black')
	#plot_Werk_ColDens(ax,werk_data,ion,'Rperp',xmax=max_r)
	#ax.text(100,17.05,re.split('galquas/|/frb',modelname)[1])
	ax.set_xticks(range(0,160,30))
	#ax.axis([0.,160.,8,18])
	return percent_above
#########################################################################

#ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
#ions = ['SiIV','CIII','OVI']
ions1 = ['CIII_977','CIII_977','CIII_977']
ions2 = ['SiIV','CIV','OVI']
#ions = ['SiIV']
#ions = ['HI','MgII','SiII']
#ions = ['SiIII','CIV','OVI']
#model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/coldens/z02/'
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/z02/'

model_gqs = ['g1q01','g1q1','g1q10']#,'g1q01','g1q1','g1q10','g1q01','g1q1','g1q10']
#model_mid = '/frbz_1kpc_500kpc_z02_'
model_mid = '/frbx_1kpc_500kpc_z02_'
#model_mid = '/frbFace_default_z02_'

werk_data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
l,u = 10.,20.
model_width = 0.0
nbins,ndraws = 500,10
max_r = 160.

#fileout = 'paper1X_scatter_matrix_Zfixed_500kpc_sSFR_wemisprof.png'
#fileout = 'paper1_scatter_matrix_Zfixed_500kpc_sSFR.png'
fileout = 'ionfrac_scatter_emis_z02_nozscale_01limit_dens.png'
xlen,ylen = 3,3 #6#2,6
fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
#fig.set_size_inches(12,4)
fig.set_size_inches(8,8)#,16)
#gs1 = gridspec.GridSpec(4, 4)
#gs1.update(wspace=0.00, hspace=0.05)
#plt.subplots_adjust(.1,.1,.9,.9,0,0.1)
ax = ax.flat
i = 0
r,dr,nrad = make_radius_array()
xL = np.linspace(-160,160,320)

f = open('percent_above.dat', 'w')

model_gqs = ['g1q01','g1q1','g1q10']
i2 = -1

for ion in ions1:
	count = 0
	i2 = i2 + 1
	while count < 3:
		print ion, count, model_gqs[count],i
		ion2 = ions2[i2]

		modelname1 = model_beg+model_gqs[count]+model_mid+ion+'.cpkl'
                profile_names1 = [model_beg+model_gqs[count]+model_mid+ion+'.cpkl',model_beg+model_gqs[count]+model_mid+ion+'.cpkl',model_beg+model_gqs[count]+model_mid+ion+'.cpkl']

		modelname2 = model_beg+model_gqs[count]+model_mid+ion2+'.cpkl'
                profile_names2 = [model_beg+model_gqs[count]+model_mid+ion2+'.cpkl',model_beg+model_gqs[count]+model_mid+ion2+'.cpkl',model_beg+model_gqs[count]+model_mid+ion2+'.cpkl']


		full_scatter_plot(modelname1,profile_names1,ion,modelname2,profile_names2,ion2,ax[i],werk_data,max_r)

		
		count = count + 1
		i = i + 1	


for j in range(len(model_gqs)):
	ax[j].set_title(model_gqs[j])

for k in range(len(ions1)):
	#ax[k*xlen+2].text(100,17.05,ions1[k]+'/'+ions2[k])
	ax[k*xlen+2].text(100,4.05,ions1[k]+'/'+ions2[k])

ax[7].set_xlabel('Impact Parameter [kpc]')
#ax[3].set_ylabel('Column Density [cm^-2]')
ax[3].set_ylabel('Surface Brightness')

plt.tight_layout()
plt.savefig(fileout)#,transparent=True)
plt.close()
f.close()



