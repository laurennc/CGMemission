import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data_lauren import *
import matplotlib.pyplot as plt
import cPickle
import matplotlib
import pylab

#patt = "/u/10/l/lnc2115/vega/data/Ryan/frbs/frb"
patt = "/u/10/l/lnc2115/vega/data/Ryan/frbs/frb"
pattend = "_5kpc_CIV_Scaled_ncut.cpkl"

rpCIVr,rpCIVmean,rpCIVmedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,1.28e-11)

plt.plot(rpCIVr,np.log10(rpCIVmean),'r-',linewidth=1.2,label='CIV Mean')
plt.plot(rpCIVr,np.log10(rpCIVmedian),'r--',linewidth=1.2,label='CIV Median')

pattend = "_5kpc_OVI_Scaled_ncut.cpkl"

rpOVIr,rpOVImean,rpOVImedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,1.92e-11)

plt.plot(rpOVIr,np.log10(rpOVImean),'b-',linewidth=1.2,label='OVI Mean')
plt.plot(rpOVIr,np.log10(rpOVImedian),'b--',linewidth=1.2,label='OVI Median')

pattend = "_5kpc_CIII_977_Scaled_ncut.cpkl"

rpCIVr,rpCIVmean,rpCIVmedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,1.28e-11)

plt.plot(rpCIVr,np.log10(rpCIVmean),'g-',linewidth=1.2,label='CIII 977 Mean')
plt.plot(rpCIVr,np.log10(rpCIVmedian),'g--',linewidth=1.2,label='CIII 977 Median')

pattend = "_5kpc_HaR.cpkl"

rpHar,rpHamean,rpHamedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,3.028e-12)

plt.plot(rpHar,np.log10(rpHamean),'k-',linewidth=1.2,label="Ha Mean")
plt.plot(rpHar,np.log10(rpHamedian),'k--',linewidth=1.2,label="Ha Median")



plt.xlabel('r (kpc)')
plt.ylabel('SB (LU)')
plt.legend()
plt.xlim(5,300)
plt.ylim(-1,6)
plt.xscale('log')
plt.savefig('metal_radial_profiles_Scaled_ncut.png')



