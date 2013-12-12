from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from radial_data_lauren import *

def ergs_sr_TO_raleighs(data_arr):
	return data_arr*3.304e11/79577.4715459

def distance_from_center(x,y,z,center):
        return ((x-center[0])**2.0+(y-center[1])**2.0+(z-center[2])**2.0)**0.5

def load_Ryan_data():
	fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
	pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
	pos = [0.40328598,0.47176743,0.46131516]
	rad = 108.0/pf['kpc']
	data = pf.h.sphere(pos,rad)
	return pf, data

def scale_by_metallicity(values,assumed_Z,wanted_Z):
	wanted_ratio = (10.**(wanted_Z))/(10.**(assumed_Z))
	return values*wanted_ratio


def plot_scatter_percentile(data,x,y,percentile,symbol,change_units=False,energy=1.):
	#data is the frb array data
	x,y = np.meshgrid(x,y)
	working_mask = np.ones(data.shape,bool)
	r = abs(x+1j*y)
	rmax = r[working_mask].max()
        dr = np.abs([x[0,0] - x[0,1]]) #* annulus_width
        radial = np.arange(rmax/dr)*dr + dr/2.
        nrad = len(radial)	

	if (change_units):
		data = data*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)

	for irad in range(nrad):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (r>=minrad) * (r<maxrad) * working_mask
		rhere = r[thisindex].flatten()
		datahere = data[thisindex].flatten()
		meanhere = data[thisindex].mean()
		lenperc = int(len(datahere)*percentile)
		idSort = np.argsort(datahere)[::-1]
		wanted = idSort[0:lenperc]
		if len(wanted) == 0:
			wanted = np.where(datahere == datahere.max())
		plt.plot(rhere[wanted],np.log10(datahere[wanted]),symbol,markersize=3.5)#,markersize=0.75)				

def make_Cloudy_table(table_index):
	hden_n_bins, hden_min, hden_max = 15, -6, 1
	T_n_bins, T_min, T_max = 51, 3, 8
#	patt = "/hpc/astrostats/astro/users/lnc2115/codes/cloudy_yt/all_lines/all_emissivity_run%i.dat"
	patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/all_lines/all_emissivity_run%i.dat"
	hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
	T=numpy.linspace(T_min,T_max, T_n_bins)
	table = np.zeros((hden_n_bins,T_n_bins))
	for i in range(hden_n_bins):
		print i
		table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]
	return hden,T,table

def make_SB_profile(filex,filey,filez,energy):
	xL = np.arange(-20,20)*10.0
	xL, yL = np.meshgrid(xL,xL)
	r = abs(xL+1j*yL)

	frbx = cPickle.load(open(filex,'rb'))
	frby = cPickle.load(open(filey,'rb'))
	frbz = cPickle.load(open(filez,'rb'))
	
	frbx = frbx*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)
	frby = frby*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)
	frbz = frbz*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)

	rp_Ralx = radial_data(frbx,x=xL,y=yL)
	rp_Raly = radial_data(frby,x=xL,y=yL)
	rp_Ralz = radial_data(frbz,x=xL,y=yL)
	
	rp_mean = (rp_Ralx.mean + rp_Raly.mean + rp_Ralz.mean)/3.0
	rp_median  = (rp_Ralx.median + rp_Raly.median + rp_Ralz.median)/3.0

	return rp_Ralx.r, rp_mean, rp_median

