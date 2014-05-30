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


res = [200,200]
#res = [40,40]

#xL = np.arange(-20,20)*10.0
#xL, yL = np.meshgrid(xL,xL)
#r = abs(xL+1j*yL)
xL2 = np.arange(-100,100,0.25)

patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/frb"
pattends = ["_250pc_CIII.cpkl","_250pc_CIII_977.cpkl","_250pc_CIV.cpkl","_250pc_OVI.cpkl","_250pc_MgII.cpkl","_250pc_SiII.cpkl","_250pc_SiIII_1207.cpkl","_250pc_SiIII_1883.cpkl","_250pc_SiIV.cpkl","_250pc_HAlpha.cpkl"]

for pattend in pattends:
	rpr,rpmean,rpmedian,rpmax,rpmin = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)

	plt.plot(rpr,np.log10(rpmedian),'k-',linewidth=1.2,label=re.split('c_|.cp',pattend)[1])

	plot_scatter_percentile(cPickle.load(open(patt+"x"+pattend,'rb')),xL2,xL2,0.1,'ko')

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
plt.savefig('bertone_frbs/bertone_scatter_profile.png')
plt.close()


