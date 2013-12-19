import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data_lauren import *
import matplotlib.pyplot as plt
import cPickle
import matplotlib
import pylab

font = {'weight' : 'bold',
	'size' : 20}

matplotlib.rc('font',**font)
F = pylab.gcf()
F.set_size_inches(8,5.5)

#fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

#res = [200,200]
res = [40,40]

xL = np.arange(-20,20)*10.0
xL, yL = np.meshgrid(xL,xL)
r = abs(xL+1j*yL)
xL2 = np.arange(-20,20)*10.0

patt = "/u/10/l/lnc2115/vega/data/Ryan/frbs/frb"

pattend = "_5kpc_HaR.cpkl"
rpHar,rpHamean,rpHamedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,3.028e-12)

plt.plot(rpHar,np.log10(rpHamean),'k-',linewidth=1.2,label='Ha')
plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'ko',change_units=True,energy=3.028e-12)


pattend = "_5kpc_CIV_Scaled_ncut.cpkl"
rpCIVr,rpCIVmean,rpCIVmedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,1.28e-11)

plt.plot(rpCIVr,np.log10(rpCIVmean),'b-',linewidth=1.2,label='CIV')
plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'bo',change_units=True,energy=1.28e-11)


pattend = "_5kpc_OVI_Scaled_ncut.cpkl"
rpOVIr,rpOVImean,rpOVImedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,1.92e-11)

plt.plot(rpOVIr,np.log10(rpOVImean),'r-',linewidth=1.2,label='OVI')
plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'ro',change_units=True,energy=1.92e-11)



pattend = "_5kpc_CIII_977_Scaled_ncut.cpkl"
#don't want to scale because already correct. So!
energy_noscale = ( 5.7e-18)/(4.*np.pi*1.87e-12)
rpCIII977r,rpCIII977mean,rpCIII977median = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,energy_noscale)

plt.plot(rpCIII977r,np.log10(rpCIII977mean),'g-',linewidth=1.2,label='CIII 977')
plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'go',change_units=False)


pattend = "_5kpc_CIII_Scaled_ncut.cpkl"
#don't want to scale because already correct. So!
energy_noscale = ( 5.7e-18)/(4.*np.pi*1.87e-12)
rpCIIIr,rpCIIImean,rpCIIImedian = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend,energy_noscale)

plt.plot(rpCIIIr,np.log10(rpCIIImean),'m-',linewidth=1.2,label='CIII')
plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'mo',change_units=False)


#THIS was when I was plottini all the points and comparing top ten percent
#rflat = r.flatten()
#rHalpha = frbx.flatten()
#plt.plot(rflat,np.log10(rHalpha),'bo',markersize=1.0)



plt.xscale('log')
plt.xlabel('Radius (kpc)',fontsize=19.,fontweight='bold')
plt.ylabel('log(Surface Brightness) (LU)',fontsize=19.,fontweight='bold')
plt.axis([8.,200.,-1,6])
plt.legend(loc='lower left',prop={'size':15})
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('all_lines_scatter_profile.png')
plt.close()


