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

create_file = True
fileout = 'galaxy_models_minusNpts.cpkl'

if create_file == True:
	
	qs = [0.01,0.05,0.5,1.,2.,10.,0.05,0.01,0.1,0.5,1.,2.,10.,0.01,1.,10.]
	gs = [0.01,0.01,0.01,0.01,0.01,0.01,1.,1.,1.,1.,1.,1.,1.,10.,10.,10.]
	
	l,u,k = 10.,20.,8.
	nbins,ndraws,model_width = 500,10,0.5
	models = []
	fileout = 'galaxy_models.cpkl'
	
	
	model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'
	model_gqs = ['g01q0.05','g01q01','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
	model_mid = '/frbz_1kpc_z02_'
	ions = ['HI','MgII','SiII','SiIII','SiIV','CIII']#,'OVI']
	data  = cPickle.load(open('werk_coldens_data.cpkl','rb'))
	
	galaxies = []
	
	count = 0
	for ID in np.unique(data['ID']):
		gal = Galaxy(ID,l,u,model_width,nbins,ndraws)
		print count, ID, len(gal.coldens)
		if len(gal.coldens) != 0.0:
			for i in range(len(model_gqs)):
				model_name = model_beg+model_gqs[i]+model_mid
				#print model_name
				gal.add_model(gs[i],qs[i],model_name,gal)
			
				for ion in ions:
					idx = gal.models[i].line_indices(ion)
					if len(idx)>0.0:
					#	print ion
						gal.models[i].add_line(ion,idx)
			
				gal.models[i].total_model_likelihood()
			
			#	for line in gal.models[i].lines:
			#		print line.ion, line.likelihood	
		
			gal.find_best_model()
			galaxies = np.append(galaxies,gal)

		count = count + 1
	print fileout	
	cPickle.dump(galaxies,open(fileout,'wb'),protocol=-1)
#Make plot of best model distribution

if create_file == False:
	galaxies = cPickle.load(open(fileout,'rb'))


image_of_models = np.zeros((7, 3))
q_map,g_map = {},{}
q_all = [0.01,0.05,0.1,0.5,1.,2.,10.]
g_all = [0.01,1.,10.]

for i in range(len(q_all)):
	q_map[q_all[i]] = i
for i in range(len(g_all)):
	g_map[g_all[i]] = i

best_count = []
for i in range(len(galaxies)):
	best_g = g_map[galaxies[i].best_g]
	best_q = q_map[galaxies[i].best_q]

	image_of_models[best_q,best_g] += 1.

print image_of_models

plt.imshow(image_of_models,interpolation='none',cmap='Greys',origin="lower")#,extent=(0,10,0.05,10),aspect='auto')
plt.xlabel('g')
plt.ylabel('q')
plt.xticks(range(len(g_all)),g_all)
plt.yticks(range(len(q_all)),q_all)
plt.colorbar()
pltout = 'galaxy_models_nbins500_ndraws10_sinOVI_detectonly.png'
print pltout
plt.savefig(pltout)
plt.close()

#c,bins,p = plt.hist(best_count,bins=len(model_gqs),alpha=0.5,facecolor='Plum')#'MediumBlue')
#bin_centers = 0.5 * np.diff(bins) + bins[:-1]
#plt.xticks(bin_centers,model_gqs,rotation='vertical')
#plt.margins(0.02)
#plt.subplots_adjust(bottom=0.15)
#plt.savefig('best_model_per_galaxy_sinOVI.png')
#plt.close()





