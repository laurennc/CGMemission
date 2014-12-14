#Class designed to hold the Lines of each model for plotting and other calculations
#Created by: Lauren
#Created on: 10/27/2014
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from Line_ORG import *

class Model:
	def __init__(self,ions,model_beg,model_mid,g,q,l,u,model_width,nbins,ndraws): #,k) for tanh
		self.lines = []
		self.l = l
		self.u = u
		#self.k = k
		self.g = g
		self.q = q
		self.nbins,self.ndraws = nbins,ndraws
		self.model_width = model_width
		self.ions = ions
		self.model_probability = -999.
		for ion in ions:
			model_name = model_beg+'z'+model_mid+ion+'dens.cpkl'
			#print ion,model_name
			line = Line(ion,model_name,l,u,model_width,nbins,ndraws)#,k)
			line.total_probability()
			self.lines = np.append(self.lines,line)
			self.model_probability = self.model_probability + line.total_likelihood
		
		self.model_probability = self.model_probability - np.log(len(self.lines))

##END##



