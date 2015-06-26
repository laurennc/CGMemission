import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import itertools
#from Line_ORG import *
#from Galaxy import *
import re
from lauren import make_SB_profile
from radial_data_lauren import *
from plotting_routines import *
from matplotlib import colors

cmap = colors.ListedColormap(['Gray','DarkTurquoise'])
bounds = [-6,-2,2]
norm = colors.BoundaryNorm(bounds,cmap.N)

#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"
pf = load(fn, file_style="%s.grid.cpu%%04i")

val, pos = pf.h.find_max('Density')

sp = SlicePlot(pf,'z','H_NumberDensity',center=pos,width=(320.,'kpc'))

sp.set_cmap('all',cmap)
sp.save('self_shielding_z1')


