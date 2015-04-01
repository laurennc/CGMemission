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
import pylab
from Line import *
from Galaxy import *

create_file = True
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
	zlen = len(model_gqs)

	IDs = np.unique(data['ID'])
	idlen = len(IDs)
#	ID_HI,ID_MgII,ID_SiII,ID_SiIII,ID_SiIV,ID_CIII,ID_OVI = [],[],[],[],[],[],[]

	imHI,imMgII,imSiII,imSiIII,imSiIV,imCIII,imOVI = np.zeros((zlen,idlen)),np.zeros((zlen,idlen)),np.zeros((zlen,idlen)),np.zeros((zlen,idlen)),np.zeros((zlen,idlen)),np.zeros((zlen,idlen)),np.zeros((zlen,idlen))

	imArrays = [imHI,imMgII,imSiII,imSiIII,imSiIV,imCIII,imOVI]
#	IDararys = [ID_HI,ID_MgII,ID_SiII,ID_SiIII,ID_SiIV,ID_CIII,ID_OVI]

	count = 0
	for ID in np.unique(data['ID']):
		gal = Galaxy(ID,l,u,model_width,nbins,ndraws)
		print count, ID, len(gal.coldens)
		if len(gal.coldens) != 0.0:
			for i in range(len(model_gqs)):
				model_name = model_beg+model_gqs[i]+model_mid
				#print model_name
				gal.add_model(gs[i],qs[i],model_name,gal)
				for j in range(len(ions)):
					idx = gal.models[i].line_indices(ions[j])
					if len(idx)>0.0:
						gal.models[i].add_line(ions[j],idx)
						#The line I want should be the last one that was added so
						#This is the line I'm least confident in I guess...
						#imArrays[j][i,count] = gal.models[i].lines[-1].likelihood
						for line in gal.models[i].lines:
							if line.ion == ions[j]:
								imArrays[j][i,count] = line.likelihood	

				gal.models[i].total_model_likelihood()
			
		
			gal.find_best_model()
			galaxies = np.append(galaxies,gal)

		count = count + 1
	print fileout	
	cPickle.dump(galaxies,open(fileout,'wb'),protocol=-1)

if create_file == False:
	galaxies = cPickle.load(open(fileout,'rb'))


gq_labels = ['g01q0.05','g01q01','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']

counting = 0
while counting < len(imArrays):
	F = pylab.gcf()
	F.set_size_inches(10,10)
	plt.imshow(imArrays[counting],interpolation='none',cmap='Greys_r',origin="lower")#,extent=(0,10,0.05,10),aspect='auto')
#	plt.xlabel('g/q values')
#	plt.ylabel('Galaxy')
#	plt.xticks(range(len(gq_labels)),gq_labels)
#	plt.yticks(range(len(IDs)),IDs)
	plt.colorbar()
	pltout = 'line_likelihoods_'+ions[counting]+'.png'
	plt.savefig(pltout)
	plt.close()
	counting = counting + 1

#END#





