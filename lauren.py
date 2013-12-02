from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

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


def plot_scatter_percentile(data,x,y,percentile):
	x,y = np.meshgrid(x,y)
	working_mask = np.ones(data.shape,bool)
	r = abs(x+1j*y)
	rmax = r[working_mask].max()
        dr = np.abs([x[0,0] - x[0,1]]) #* annulus_width
        radial = np.arange(rmax/dr)*dr + dr/2.
        nrad = len(radial)	

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
		plt.plot(rhere[wanted],np.log10(datahere[wanted]),'go',markersize=3.5)#,markersize=0.75)				

def make_Cloudy_table(table_index):
	hden_n_bins, hden_min, hden_max = 15, -6, 1
	T_n_bins, T_min, T_max = 51, 3, 8
	patt = "/hpc/astrostats/astro/users/lnc2115/codes/cloudy_yt/all_lines/all_emissivity_run%i.dat"
	hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
	T=numpy.linspace(T_min,T_max, T_n_bins)
	table = np.zeros((hden_n_bins,T_n_bins))
	for i in range(hden_n_bins):
		table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]
	return hden,T,table

def dummy():
	return	


