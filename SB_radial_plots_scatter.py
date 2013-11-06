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

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

#res = [200,200]
res = [40,40]

xL = np.arange(-20,20)*10.0
xL, yL = np.meshgrid(xL,xL)
r = abs(xL+1j*yL)
xL2 = np.arange(-20,20)*10.0

frbx = cPickle.load(open('frbx_5kpc_HaR.cpkl','rb'))
frby = cPickle.load(open('frby_5kpc_HaR.cpkl','rb'))
frbz = cPickle.load(open('frbz_5kpc_HaR.cpkl','rb'))

rp_Ralx = radial_data(frbx,x=xL,y=yL)
rp_Raly = radial_data(frby,x=xL,y=yL)
rp_Ralz = radial_data(frbz,x=xL,y=yL)

rp_Ralmean = (rp_Ralx.mean + rp_Raly.mean + rp_Ralz.mean)/3.0
rp_Ralmed  = (rp_Ralx.median + rp_Raly.median + rp_Ralz.median)/3.0


frbVx = cPickle.load(open('frbx_5kpc_HaVR.cpkl','rb'))
frbVy = cPickle.load(open('frby_5kpc_HaVR.cpkl','rb'))
frbVz = cPickle.load(open('frbz_5kpc_HaVR.cpkl','rb'))

rp_RalVx = radial_data(frbVx,x=xL,y=yL)
rp_RalVy = radial_data(frbVy,x=xL,y=yL)
rp_RalVz = radial_data(frbVz,x=xL,y=yL)

rp_RalVmean = (rp_RalVx.mean + rp_RalVy.mean + rp_RalVz.mean)/3.0
rp_RalVmed  = (rp_RalVx.median + rp_RalVy.median + rp_RalVz.median)/3.0


vs_x = (np.array([0.0025,0.022,0.045,0.146])/0.73)*10**3.0
vs_y = np.array([10**4.,10**3.,10**2.,10**1.])/((10**6.0)/(4.*np.pi))

rflat = r.flatten()
rHalpha = frbx.flatten()
plt.plot(rflat,np.log10(rHalpha),'bo',markersize=1.0)
plt.plot(rp_Ralx.r,np.log10(rp_Ralmean),'b-',linewidth=1.2,label='Mean')
plt.plot(rp_Ralx.r,np.log10(rp_Ralmed),'r-',linewidth=1.2,label='Median')
plt.plot(rp_RalVx.r,np.log10(rp_RalVmean),'m-',linewidth=1.2,label='No Dense Gas: Mean')
plt.plot(vs_x,np.log10(vs_y),'m--',linewidth=1.75,label='van de Voort and Schaye (2013)')
plot_scatter_percentile(frbx,xL2,xL2,0.1)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Radius (kpc)',fontsize=20.,fontweight='bold')
plt.ylabel('log(Surface Brightness) (R)',fontsize=20.,fontweight='bold')
#plt.axis([8,200,10**-5,10])
plt.axis([8.,200.,-5,1])
#plt.legend(loc='upper right',prop={'size':15})
plt.legend(loc='lower left',prop={'size':15})
plt.tight_layout()
#plt.axes().set_aspect('equal')
#plt.locator_params(axis='x',nbins=3)
plt.savefig('scatter_5kpc_x_voort2L55.png')
plt.close()
