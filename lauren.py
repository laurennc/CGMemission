from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

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


def HAlpha_Emission_Arc(t,hden):
	#assuming that n_e/n_H = 1.15 and n_HII/n_H = 0.9936
	Tfour = (10**t)/1.0e4
	exponent = -0.942-0.031*np.log(Tfour)
	return (2.82e-26)*(1.15*0.9936)*(2.3528e-11)*((10**hden)**2.0)*Tfour**(exponent)

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


