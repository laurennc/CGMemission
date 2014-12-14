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

create_file = False
fileout = 'galaxy_models_detectonly.cpkl'

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
	ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
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
	print 'in false'
	fileout = 'galaxy_models.cpkl'
	galaxies = cPickle.load(open(fileout,'rb'))

model_gqs = ['g01q0.05','g01q01','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']
ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']


count = 0
while count < len(galaxies):
	print count
	gal = galaxies[count]
	
	fig,axs = plt.subplots(3,3)
	fig.set_size_inches(13.5,13.5)
	ax2 = axs.flat

	i = 0
	while i < len(ions):
		ion = ions[i]
		ax2[i].set_title(ion)
	
		if len(np.where(gal.ions == ion)[0]) > 0:
			
			x = np.arange(len(model_gqs))+0.5
			y = []
			for model in gal.models:
				for line in model.lines:
					if line.ion == ion:
						y = np.append(y,line.likelihood)
			
			ax2[i].plot(x,y,'ob-')
			ax2[i].set_ylabel('Log Likelihood')
			ax2[i].set_xticks(x, model_gqs)

		i = i + 1

	plt.savefig('galaxy'+str(count)+'_ion_likelihoods.png')
	plt.close()

	count = count + 1

