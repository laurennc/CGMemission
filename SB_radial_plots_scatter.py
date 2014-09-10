import matplotlib as mpl
mpl.use('agg')
import numpy as np
import cPickle
from lauren import make_SB_profile
from lauren import plot_scatter_percentile
from radial_data_lauren import *
import matplotlib.pyplot as plt
import matplotlib
import pylab
import re

font = {'weight' : 'bold', 'size' : 20}
matplotlib.rc('font',**font)
#F = pylab.gcf()
#F.set_size_inches(8,5.5)


#res = [200,200]
#res = [40,40]
#xL2 = np.arange(-100,100,0.25)
xL2 = np.linspace(-160,160,320)


patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/g1q01/frb"
#patt =  "/u/10/l/lnc2115/vega/repos/CGMemission/euvb_frbs/emis/frb"
pattends = ["_1kpc_z02_CIII.cpkl","_1kpc_z02_CIII_977.cpkl","_1kpc_z02_CIV.cpkl","_1kpc_z02_OVI.cpkl","_1kpc_z02_MgII.cpkl","_1kpc_z02_SiII.cpkl","_1kpc_z02_SiIII_1207.cpkl","_1kpc_z02_SiIV.cpkl","_1kpc_z02_HAlpha.cpkl"]#,"_1kpc_z02_SiIII_1883.cpkl","_1kpc_z02_SiIV.cpkl","_1kpc_z02_HAlpha.cpkl"]

#iax = 221
iax = 331
fig = plt.figure(figsize=(13.5,13.5))
i = 0

colors = ['r','b','g','m']

for i in range(len(pattends)):
	pattend = pattends[i]
	print pattend
	#if (i%3 == 0) * (i!=9):
	#	ax1 = fig.add_subplot(iax)
	#	ax1.set_xscale('log')
	#	ax1.set_xlabel('Radius (kpc)')
	#	ax1.set_ylabel('log(Surface Brightness) (LU)')
	#	ax1.axis([8.,200.,-4,7])
	#	iax = iax + 1
	
	#colorhere = colors[i%3]
	#if i == 9: colorhere = colors[3]

	ax1 = fig.add_subplot(iax)

	rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)

	idr = np.where(rpr <= 160.0)

	plot_scatter_percentile(ax1,cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.05,'k.',maxr=160.)
	ax1.plot(rpr[idr],np.log10(rpmedian[idr]),color='Goldenrod',linewidth=2.2,label=re.split('z02_|.cp',pattend)[1])

	#ax1.set_xscale('log')
        #ax1.set_xlabel('Radius (kpc)')
        #ax1.set_ylabel('log(Surface Brightness) (LU)')
        #ax1.axis([8.,200.,-4,7])
        ax1.legend(loc='upper right',prop={'size':12})
	ax1.axis([0.,160.,-8,8])
	iax = iax + 1
	i = i + 1

plt.gcf().subplots_adjust(bottom=0.15)
#plt.savefig('euvb_frbs/euvb_scatter_profile.png')
#plt.savefig('euvb_frbs/euvb_scatter_profile2.png')
plt.savefig('bertone_frbs/g1q01_scatter_limits.png')
plt.close()


