import matplotlib as mpl
mpl.use('agg')
from lauren import make_SB_profile
from radial_data_lauren import *
import matplotlib.pyplot as plt
import cPickle
import matplotlib
import pylab
import re
from plotting_routines import *

def basic_profile_vals(patt,pattend):
	rpr,rpmean,rpmedian,rpmax,rpmin = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)
	return rpr,rpmean,rpmedian,rpmax,rpmin

def make_plots():
	print "Don't forget to make sure radial_profile code in lauren.py is set to the same array size as the frbs I'm passing in!"
	
	font = {'weight' : 'bold','size' : 18}
	matplotlib.rc('font',**font)

#	patt = "/u/10/l/lnc2115/vega/repos/CGMemission/euvb_frbs/coldens/frb"
	patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/factor10/frb"
	pattends = ["_1kpc_z02_CIIIdens.cpkl","_1kpc_z02_CIVdens.cpkl","_1kpc_z02_OVIdens.cpkl","_1kpc_z02_MgIIdens.cpkl","_1kpc_z02_SiIIdens.cpkl","_1kpc_z02_SiIIIdens.cpkl","_1kpc_z02_SiIVdens.cpkl","_1kpc_z02_HIdens.cpkl"]	

	iax = 331
	fig = plt.figure(figsize=(13.5,13.5))	

	data = cPickle.load(open('werk_coldens_data.cpkl','rb'))

	#xL = np.arange(-100,100,0.25)
	xL = np.linspace(-160,160,320)
	maxnow = 160.
	xL, yL = np.meshgrid(xL,xL)
	r = abs(xL+1j*yL)
	dr = np.abs([xL[0,0] - xL[0,1]])
	#radial = np.arange(r.max()/dr)*dr + dr/2.
	radial = np.arange(maxnow/dr)*dr + dr/2
	nrad = len(radial)
	mslope = (0.0075-0.07)/nrad

	for pattend in pattends:
		ax1 = fig.add_subplot(iax)
               	ax1.axis([0.,160.,8,18])
		
	#loop through to plot all of the points, decreasing alpha
		frbarr = np.array(cPickle.load(open(patt+"x"+pattend,'rb')))
		for irad in range(nrad):
			minrad = irad*dr
			maxrad = minrad + dr
			thisindex = (r>=minrad) * (r<maxrad)
			alphahere = mslope*irad + 0.07
			ax1.plot(r[thisindex],np.log10(frbarr[thisindex]),'k.',alpha=alphahere)

		print pattend

		rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)
		idr = np.where(rpr <= maxnow)[0]
		ax1.plot(rpr[idr],np.log10(rpmedian[idr]),'c-',linewidth=2.2,label=re.split('z02_|dens',pattend)[1])
		
		key = re.split('z02_|dens',pattend)[1]
		plot_Werk_ColDens(ax1,data,key,'Rperp',xmax=maxnow)
		ax1.text(120,17.05,re.split('z02_|dens',pattend)[1])
		ax1.set_xticks(range(0,160,30))
		#ax1.legend()
		#ax1.set_xlim(0,160)
		#ax1.set_xscale('log')
		iax = iax + 1
		
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.savefig('bertColDens_factor10_limits.png')
	plt.close()


