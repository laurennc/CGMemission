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
	elif key=='13kpc':
		xL = np.linspace(-160,160,24)
        else:
                print 'Failure in the key'
        maxnow = 160.
        xL2, yL = np.meshgrid(xL,xL)
        r = abs(xL2+1j*yL)
        dr = np.abs([xL2[0,0] - xL2[0,1]])
        radial = np.arange(maxnow/dr)*dr + dr/2
        nrad = len(radial)
        return xL,r, dr, nrad

def full_scatter_plot(modelnames,ax,res_key,max_r,znow,color,comoving=False):
        xL,r,dr,nrad = make_radius_array(res_key)
        ##EDITED HERE FROM 0 TO 1
        frbarr = np.array(cPickle.load(open(modelnames[0],'rb')))
        rpr,rpmean,rpmedian,rpmax,rpmin = basic_profile_vals(modelnames,xL,z=znow)

        idr = np.where(rpr <= 160.0)
        if comoving:
                ax.plot(rpr[idr]*(1.+znow),np.log10(rpmedian[idr]),color='k',linewidth=2.2)
        else:
                #ax.plot(rpr[idr],np.log10(rpmedian[idr]),color=color,linewidth=2.2,label=res_key) 
        	ax.plot(rpr[idr],np.log10(rpmean[idr]),color=color,linewidth=2.2,label=res_key) 
	return

model_gqs = ['g1q1']
res_keys = ['1kpc','5kpc','13kpc','25kpc']
colors = ['#e7298a','#1b9e77','#d95f02','#7570b3','#e7298a']
redshift_keys = ['z02']
znow = [0.]
ions = ['CIII_977']
ionsC = ['CIII']
xlen,ylen = 1,2
figxlen,figylen = 6,3
fileout = 'res_test_z02_CIII_avg.pdf'
model_beg_emis = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'+redshift_keys[0]+'/'
model_beg_col = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/coldens/'+redshift_keys[0]+'/' 
max_r = 160.

xlen,ylen = 2,1
fig,ax = plt.subplots(ylen,xlen)
fig.set_size_inches(6,3)
ax = ax.flat


i = 0
while i < len(res_keys):
	print i, res_keys[i]
	
	modelnames = [model_beg_emis+model_gqs[0]+'/frbx_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ions[0]+'.cpkl',model_beg_emis+model_gqs[0]+'/frby_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ions[0]+'.cpkl',model_beg_emis+model_gqs[0]+'/frbz_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ions[0]+'.cpkl']

	full_scatter_plot(modelnames,ax[1],res_keys[i],max_r,znow[0],colors[i])

	modelnames = [model_beg_col+model_gqs[0]+'/frbx_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ionsC[0]+'dens.cpkl',model_beg_col+model_gqs[0]+'/frby_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ionsC[0]+'dens.cpkl',model_beg_col+model_gqs[0]+'/frbz_'+res_keys[i]+'_500kpc_'+redshift_keys[0]+'_'+ionsC[0]+'dens.cpkl']


	full_scatter_plot(modelnames,ax[0],res_keys[i],max_r,znow[0],colors[i])


	i = i + 1


ax[0].set_ylabel('log(Column Density')
ax[1].set_ylabel('log(Surface Brightness')
ax[0].set_xlabel('Radius')
ax[1].set_xlabel('Radius')
ax[0].set_xticks([0,50,100,150],['0','50','100','150'])
ax[1].set_xticks([0,50,100,150],['0','50','100','150'])

plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(fileout)
plt.close()





