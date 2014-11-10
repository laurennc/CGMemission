#Class designed to hold the Lines of each model for plotting and other calculations
#Created by: Lauren
#Created on: 10/27/2014
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
from Line import *
from Galaxy import *

qs = [0.01,1,10.]
gs = [0.01,1,10.]
l,u,k = 0.,20.,8.
model_width = 0.5
models = []

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
#model_gqs = ['g01q0.05','g01q01','g01q0.1','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
model_gqs = ['g01q0.05','g01q01','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
model_mid = '/frbz_1kpc_z02_'
ions = ['HI','MgII','SiII','SiIII','SiIV','CIII']#,'OVI']

data  = cPickle.load(open('werk_coldens_data.cpkl','rb'))

galaxies = []

count = 0
for ID in np.unique(data['ID']):
	print count, ID
	gal = Galaxy(ID,l,u,model_width)
	for i in range(len(model_gqs)):
		model_name = model_beg+model_gqs[i]+model_mid
		gal.add_model(109.,109.,model_name,gal)

		for ion in ions:
			idx = gal.models[i].line_indices(ion)
			gal.models[i].add_line(ion,idx)
	
		gal.models[i].total_model_likelihood()
	
	gal.find_best_model()
	galaxies = np.append(galaxies,gal)

	count = count + 1

#print 'done'
#print galaxies

best_count = []
for i in range(len(galaxies)):
	best_count = np.append(best_count,galaxies[i].best_model)

c,bins,p = plt.hist(best_count,bins=len(model_gqs),alpha=0.5,facecolor='Plum')#'MediumBlue')
bin_centers = 0.5 * np.diff(bins) + bins[:-1]
plt.xticks(bin_centers,model_gqs,rotation='vertical')
#plt.margins(0.02)
plt.subplots_adjust(bottom=0.15)
plt.savefig('best_model_per_galaxy_sinOVI.png')
plt.close()





