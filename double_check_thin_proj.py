import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import itertools

###OK WHAT I WANT THIS TO DO IS PLOT THE FRBS SIDE BY SIDE WITH THE SAME COLOR BAR SO I CAN GET BETTER IDEA AS TO WHATS HAPPENING

fileout = 'frb_comparison_g1q1.png'
xlen,ylen = 2,3
fig,ax = plt.subplots(ylen,xlen)
ax = ax.flat

i = 0

model_gqs = ['g1q1']
ions = ['SiIV','CIII','OVI']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
model_mid1 = '/frbz_1kpc_320kpc_z02_'
model_mid2 = '/frbz_1kpz_z02_'





