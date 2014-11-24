#Class designed to estimate the probability of a certain ionizing model given the data of Werk et al. 2013
#This uses a convultion of a Gaussian and a step function to represent the upper and lower limits
#It can also handle adding uncertainty to the data to represent physical scatter if desired
#Created by: Lauren
#Created on: 11//4/2014
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special
import matplotlib.patches

class Line:
	def __init__(self,ion,idx,model_name,galaxy,radii):
		self.ion = ion
                self.idx = idx
		self.model_name = model_name+ion+'dens.cpkl'	
		self.model = np.log10(cPickle.load(open(self.model_name,'rb')))
		self.galaxy = galaxy
		self.radii = radii
		self.nbins = galaxy.nbins
	
		#So the functions run without having to edit them
		self.l,self.u,self.model_width = galaxy.l,galaxy.u,galaxy.model_width
		self.coldens = galaxy.coldens[idx]
		self.rperp = galaxy.rperp[idx]		
		self.label = galaxy.label[idx]
		self.err = galaxy.err[idx]
		#Finally, create a variable that will hold the total likelihood of the model for this line!
		self.likelihood = 0.0
	

	#l and u are the upper and lower limits of the model range that we are currently considering
	def upper_limit(self,bin_idx,hist,del_x):
		return sum(hist[:bin_idx+1]*del_x)

	def lower_limit(self,bin_idx,hist,del_x):
		return sum(hist[bin_idx:]*del_x)

	def model_probability(self,data_idx):
		idr = np.where((self.radii >= self.rperp[data_idx]-2.5) & (self.radii <= self.rperp[data_idx]+2.5))
		model_pts = self.model[idr]
		#Here I'm just binning the data and not adding in any measurement error or smoothing error
		#perhaps one way to not add gaussian is to add +x for real bin and +y for a given smoothing to whatever number of bins we want
		hist,edges = np.histogram(model_pts,bins=self.nbins,range=(self.l,self.u),density=True)
		
		#del_x = (self.u-self.l)/self.nbins
		del_x = edges[1]-edges[0]
		bin_want = int((self.coldens[data_idx]-self.l)/del_x)

		if self.label[data_idx] == 'n':
			prob = hist[bin_want]
		elif self.label[data_idx] == 'l':
			prob = self.lower_limit(bin_want,hist,del_x)
		elif self.label[data_idx] == 'u':
			prob = self.upper_limit(bin_want,hist,del_x)
		else:
			print 'WARNING: SOMETHING IS WRONG WITH THE DATA'
		if prob == 0.0:
			prob = 1e-45

		#print self.ion,prob
		return np.log(prob)

	def total_probability(self):
		for data_idx in range(len(self.idx)):
			self.likelihood = self.likelihood + self.model_probability(data_idx)
		self.likelihood = self.likelihood - np.log(len(self.idx))
		return
		
	
	
	def plot_pdf(self,data_idx,ax):
		data_idx = 0
		idr = np.where((self.radii >= self.rperp[data_idx]-2.5) & (self.radii <= self.rperp[data_idx]+2.5))
                model_pts = self.model[idr]
                hist,edges = np.histogram(model_pts,bins=self.nbins,range=(self.l,self.u),density=True)

		if self.label[data_idx] == 'u':
			rect1 = matplotlib.patches.Rectangle((edges[0],0), self.coldens[data_idx]-edges[0], np.max(hist), color='DeepPink',alpha=0.15)
			ax.add_patch(rect1)
		if self.label[data_idx] == 'l':
			rect1 = matplotlib.patches.Rectangle((self.coldens[data_idx],0), edges[-1]-self.coldens[data_idx], np.max(hist), color='DeepPink',alpha=0.15)
			ax.add_patch(rect1)

		ax.plot(edges[1:],hist,color='Gray')
                ax.axvline(x=self.coldens[data_idx],color='DeepPink',linewidth=2.2,linestyle='-.')

		prob = self.model_probability(data_idx)	
		prob = np.exp(prob)

		ax.set_xlabel('Values')
		ax.set_ylabel('Density')
		ax.set_title('p='+str(prob))
		
		return
