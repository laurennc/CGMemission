#from yt.mods import *
import numpy as np
import numpy
import cPickle
import matplotlib.pyplot as plt
from scipy import interpolate
from radial_data_lauren import *
import triangle


def scale_by_metallicity(values,assumed_Z,wanted_Z):
	wanted_ratio = (10.**(wanted_Z))/(10.**(assumed_Z))
	return values*wanted_ratio


def make_Cloudy_table(table_index):
	hden_n_bins, hden_min, hden_max = 15, -6, 1
	T_n_bins, T_min, T_max = 51, 3, 8
	patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/grid_galquas/emis/z0/g1q1/g1q1_run%i.dat"

	hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
	T=numpy.linspace(T_min,T_max, T_n_bins)
	table = np.zeros((hden_n_bins,T_n_bins))
	for i in range(hden_n_bins):
		table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]
	return hden,T,table



