import matplotlib as mpl
mpl.use('agg')
from lauren import *
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
	
	font = {'weight' : 'bold',
	'size' : 18}
	
	matplotlib.rc('font',**font)
	#F = pylab.gcf()
	#F.set_size_inches(8.75,5.5)

	patt = "/u/10/l/lnc2115/vega/repos/CGMemission/ColDensFRBs/frb"

	pattends = ["_250pc_CIII_coldens.cpkl","_250pc_CIV_coldens.cpkl","_250pc_OVI_coldens.cpkl","_250pc_MgII_coldens.cpkl","_250pc_SiII_coldens.cpkl","_250pc_SiIII_coldens.cpkl","_250pc_SiIV_coldens.cpkl","_250pc_HI_coldens.cpkl"]	

	ions = ['C III','C IV','O VI','Mg II',]

	iax = 331
	fig = plt.figure(figsize=(13.5,13.5))	

	data = cPickle.load(open('werk_coldens_data.cpkl','rb'))

	for pattend in pattends:
		print pattend
		rpr,rpmean,rpmedian,rpmax,rpmin,rpstd = make_SB_profile(patt+"x"+pattend,patt+"y"+pattend,patt+"z"+pattend)
		ax1 = fig.add_subplot(iax)
		ax1.plot(rpr,np.log10(rpmedian),'k-',linewidth=1.2,label=re.split('c_|_c',pattend)[1])
		key = re.split('c_|_c',pattend)[1]
		ax1.fill_between(rpr,np.log10(rpmedian-rpstd),np.log10(rpmedian+rpstd),color='gray',alpha='0.30')
		plot_Werk_ColDens(ax1,data,key,'Rperp')
		ax1.set_xscale('log')
		#ax1.set_ylim(12,20)
		ax1.set_xlim(8,200)
		ax1.text(12,15.05,re.split('c_|_c',pattend)[1])
		#ax1.set_xlabel('Radius (kpc)')
		#ax1.set_ylabel('Column Density [cm^-2]')	
		iax = iax + 1
		
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.savefig('coldens_profs_werk_std.png')
	plt.close()


