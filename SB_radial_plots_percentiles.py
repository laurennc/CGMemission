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

def combine_datasets_perline(patt,middles,pattend):
	first=True
	for mid in middles:
		filein = patt+mid+pattend
                datain = cPickle.load(open(filein,'rb'))
		if first==True:
			data=datain	
			first=False
		else:
			data =  np.concatenate((data,datain),axis=1)

	return data

def create_radii_grid(middles,maxr,resolution=1):
	xL = np.linspace(-maxr,maxr,maxr*2*resolution)
	xL,yL = np.meshgrid(xL,xL)
	for i in range(len(middles)):
		if i == 0:
			xL2,yL2 = xL,yL
		else:
			xL2 = np.concatenate((xL2,xL),axis=1)
			yL2 = np.concatenate((yL2,yL),axis=1)
		i = i + 1
	return xL2,yL2

def make_plots():
	print "Don't forget to make sure radial_profile code in lauren.py is set to the same array size as the frbs I'm passing in!"
	
	font = {'weight' : 'bold','size' : 18}
	matplotlib.rc('font',**font)

	#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/control/frb"
	patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q01/frb"
	pattends = ["_1kpc_z02_CIIIdens.cpkl","_1kpc_z02_CIVdens.cpkl","_1kpc_z02_OVIdens.cpkl","_1kpc_z02_MgIIdens.cpkl","_1kpc_z02_SiIIdens.cpkl","_1kpc_z02_SiIIIdens.cpkl","_1kpc_z02_SiIVdens.cpkl","_1kpc_z02_HIdens.cpkl"]	

	iax = 331
	fig = plt.figure(figsize=(13.5,13.5))	

	data = cPickle.load(open('werk_coldens_data.cpkl','rb'))

	maxnow = 160.
	middles = ["x","y","z"]
	
	xL,yL = create_radii_grid(middles,maxnow)

	for pattend in pattends:
		ax1 = fig.add_subplot(iax)
               	ax1.axis([0.,160.,8,18])
		print pattend
	
		line_data = combine_datasets_perline(patt,middles,pattend)
		print xL.shape,yL.shape,line_data.shape
		rp = radial_data(line_data,x=xL,y=yL)

		idr = np.where(rp.r <= maxnow)[0]
		#ax1.plot(rp.r[idr],np.log10(rp.median[idr]),'c-',linewidth=2.2,label=re.split('z02_|dens',pattend)[1])
		ax1.plot(rp.r[idr],np.log10(rp.median[idr]),'-',color='SlateGray',linewidth=2.2,label=re.split('z02_|dens',pattend)[1])
		ax1.fill_between(rp.r[idr],np.log10(rp.q25[idr]),np.log10(rp.q75[idr]),color='SlateGray',alpha=0.35)

		key = re.split('z02_|dens',pattend)[1]
		plot_Werk_ColDens(ax1,data,key,'Rperp',xmax=maxnow)
		ax1.text(120,17.05,re.split('z02_|dens',pattend)[1])
		ax1.set_xticks(range(0,160,30))
		iax = iax + 1
	print key		
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.savefig('g01q01_percentiles.png')
	plt.close()


