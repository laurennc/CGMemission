import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data_lauren import *
import matplotlib.pyplot as plt
import cPickle
import matplotlib
import pylab

patt = "/u/10/l/lnc2115/vega/data/Ryan/frbs/frb"
pattend = "_5kpc_CIV_Scaled.cpkl"

rpCIVr,rpCIVmean,rpCIVmedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)

plt.plot(rpCIVr,np.log10(rpCIVmean),'r-',linewidth=1.2,label='CIV Mean')
plt.plot(rpCIVr,np.log10(rpCIVmedian),'r--',linewidth=1.2,label='CIV Median')

pattend = "_5kpc_OVI_Scaled.cpkl"

rpOVIr,rpOVImean,rpOVImedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)

plt.plot(rpOVIr,np.log10(rpOVImean),'b-',linewidth=1.2,label='OVI Mean')
plt.plot(rpOVIr,np.log10(rpOVImedian),'b--',linewidth=1.2,label='OVI Median')

plt.savefig('metal_radial_profiles_Scaled.png')
