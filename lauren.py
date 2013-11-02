from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def ergs_sr_TO_raleighs(data_arr):
	return data_arr*3.304e11/79577.4715459

def distance_from_center(x,y,z,center):
        return ((x-center[0])**2.0+(y-center[1])**2.0+(z-center[2])**2.0)**0.5

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

#	print 'nrad is',nrad	

	for irad in range(nrad):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (r>=minrad) * (r<maxrad) * working_mask
		rhere = r[thisindex].flatten()
		datahere = data[thisindex].flatten()
		meanhere = data[thisindex].mean()
		idx = np.where(datahere > meanhere)[0]
		lenperc = int(len(idx)*percentile)
		idSort = np.argsort(datahere)[::-1]
		wanted = idSort[0:lenperc]
		if len(wanted) == 0:
			wanted = np.where(datahere == datahere.max())
#		print wanted
		plt.plot(rhere[wanted],datahere[wanted],'go')#,markersize=0.75)				


