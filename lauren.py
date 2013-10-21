from yt.mods import *
import numpy as np

def ergs_sr_TO_raleighs(data_arr):
	return data_arr*3.304e11/79577.4715459

def distance_from_center(x,y,z,center):
        return ((x-center[0])**2.0+(y-center[1])**2.0+(z-center[2])**2.0)**0.5

#def make_plots(xs,labels=None, **kwargs):
#	K = len(xs)
#	factor = 2.0 # size of one side of one panel
#	lbdim = 0.5 * factor # size of left/bottom margin
#	trdim = 0.05 * factor # size of top/right margin
#	whspace = 0.05 # w/hspace size
#	plotdim = factor * K + factor * (K - 1.) * whspace
#	dim = lbdim + plotdim + trdim
#
#	fig, axes = plt.subplots(3, 6, figsize=(dim, dim))
#
#	lb = lbdim / dim
#	tr = (lbdim + plotdim) / dim
#	fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
#                       wspace=whspace, hspace=whspace)
#
