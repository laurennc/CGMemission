import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import itertools
import matplotlib.gridspec as gridspec
import re
from lauren_holding import make_SB_profile
from radial_data_lauren import *
from plotting_routines import *
from matplotlib import colors



def basic_profile_vals(modelnames,xL,z=0.):
        rpr,rpmean,rpmedian,rpmax,rpmin,rp_std = make_SB_profile(modelnames[0],modelnames[1],modelnames[2],xL,z=z)
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

def plot_profile(modelnames,ion,ax,res_key,max_r,znow,comoving=False,color='k'):
        xL,r,dr,nrad = make_radius_array(res_key)
        ##EDITED HERE FROM 0 TO 1
        frbarr = np.array(cPickle.load(open(modelnames[0],'rb')))
        rpr,rpmean,rpmedian,rpmax,rpmin = basic_profile_vals(modelnames,xL,z=znow)

        idr = np.where(rpr <= 160.0)
        if comoving:
                ax.plot(rpr[idr]*(1.+znow),np.log10(rpmedian[idr]),color='k',linewidth=2.2)
        else:
                ax.plot(rpr[idr],np.log10(rpmedian[idr]),color=color,linewidth=2.2,label=ion)#previously Goldenrod  
        print np.log10(rpmedian[99]),np.log10(rpmean[99])
        #print np.max(rpr[idr]),np.max(np.log10(rpmedian[idr]))
        ax.set_xticks(range(0,160,30))
        #ax.set_xlim(0,160)
        ax.axis([0.,160.,-6.,6.5])
        return


ions = ['SiIV','CIII_977','CIV','SiIII_1207','OVI']
redshift_key = 'z02'
znow = 0.

emis = cPickle.load(open('cloudywerk.cpkl','rb'))
gals = cPickle.load(open('werk_galaxy_properties.cpkl','rb'))

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'+redshift_key+'/' ##CHANGED FORM Z02
model_gqs = ['g1q01','g1q1','g1q10']
res_keys = ['1kpc','1kpc','1kpc']
max_r = 160.


fileout = 'medprofs_'+redshift_key+'_nozscale_SiIII_1207.pdf'
xlen,ylen = 3,1
fig,ax = plt.subplots(ylen,xlen,sharey=True)
fig.set_size_inches(9,3)
ax = ax.flat
i = 0
colors = ['b','g','r','k','m','c']

for ion in ions:
	#if ion == 'CIII_977':
	#	ionhere  = 'CIII'
	#else:
	#	ionhere = ion

	count = 0
	while count < 3:
		print ion, count, model_gqs[count],i
                modelnames = [model_beg+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_key+'_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frby_'+res_keys[count]+'_500kpc_'+redshift_key+'_'+ion+'.cpkl',model_beg+model_gqs[count]+'/frbz_'+res_keys[count]+'_500kpc_'+redshift_key+'_'+ion+'.cpkl']


		plot_profile(modelnames,ion,ax[count],res_keys[count],max_r,znow,comoving=False,color=colors[i])
		count = count + 1

	i = i + 1

ax[1].set_xlabel('Impact Parameter [kpc]')
ax[0].set_ylabel('log(SB)')
ax[2].legend()

plt.tight_layout()
plt.savefig(fileout)
plt.close()


