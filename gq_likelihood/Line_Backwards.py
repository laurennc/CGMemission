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

class Line:
	def __init__(self,ion,idx,model_name,galaxy,radii):
		self.ion = ion
                self.idx = idx
		self.model_name = model_name+ion+'dens.cpkl'	
		self.model = np.log10(cPickle.load(open(self.model_name,'rb')))
		self.galaxy = galaxy
		self.radii = radii
	
		#So the functions run without having to edit them
		self.l,self.u,self.model_width = galaxy.l,galaxy.u,galaxy.model_width
		self.coldens = galaxy.coldens[idx]
		self.rperp = galaxy.rperp[idx]		
		self.label = galaxy.label[idx]
		self.err = galaxy.err[idx]
		#Finally, create a variable that will hold the total likelihood of the model for this line!
		self.likelihood = 0.0
	

	#l and u are the upper and lower limits of the model range that we are currently considering
	def normalize_upper_limit(self,data_pt,sigma):
		a = (self.u-data_pt)/2.
		b = (self.l-data_pt)/2.
		c = ((self.u-data_pt)/(np.sqrt(2.)*sigma))*scipy.special.erf((self.u-data_pt)/(np.sqrt(2)*sigma)) + np.exp(-((self.u-data_pt)/(np.sqrt(2.)*sigma))**2.0)/np.sqrt(np.pi)
		d = ((self.l-data_pt)/(np.sqrt(2.)*sigma))*scipy.special.erf((self.l-data_pt)/(np.sqrt(2)*sigma)) + np.exp(-((self.l-data_pt)/(np.sqrt(2.)*sigma))**2.0)/np.sqrt(np.pi)
		c = c*np.sqrt(2.)*sigma/2.
		d = d*np.sqrt(2.)*sigma/2.
		return (1.-(1e-18*(self.u-self.l)))/(a-b+d-c)

	def upper_limit(self,model_x,data_pt,sigma):
		K = self.normalize_upper_limit(data_pt,sigma)
		return (K/2.0)*(1.-scipy.special.erf((model_x-data_pt)/(np.sqrt(2)*sigma)))

	def normalize_lower_limit(self,data_pt,sigma):
		a = -(data_pt-self.u)/2.
		b = -(data_pt-self.l)/2.
		c = ((data_pt-self.u)/(np.sqrt(2.)*sigma))*scipy.special.erf((data_pt-self.u)/(np.sqrt(2)*sigma)) + np.exp(-((data_pt-self.u)/(np.sqrt(2.)*sigma))**2.0)/np.sqrt(np.pi)
		d = ((data_pt-self.l)/(np.sqrt(2.)*sigma))*scipy.special.erf((data_pt-self.l)/(np.sqrt(2)*sigma)) + np.exp(-((data_pt-self.l)/(np.sqrt(2.)*sigma))**2.0)/np.sqrt(np.pi)
		c = -1.*c*np.sqrt(2.)*sigma/2.
		d = -1.*d*np.sqrt(2.)*sigma/2.
		return (1.-(1e-18*(self.u-self.l)))/(a-b+d-c)

	def lower_limit(self,model_x,data_pt,sigma):
		K = self.normalize_lower_limit(data_pt,sigma)
		return (K/2.0)*(1.-scipy.special.erf((-model_x+data_pt)/(np.sqrt(2)*sigma)))

	def detection(self,model_x,data_pt,data_err):
		return np.exp(-np.power(model_x-data_pt,2.)/(2.*np.power(data_err,2.)))

#
	def model_probability(self,data_idx):
		#print self.rperp[data_idx]
		idr = np.where((self.radii >= self.rperp[data_idx]-0.05) & (self.radii <= self.rperp[data_idx]+0.05))
		#idr = np.where((self.radii >= self.rperp[data_idx]-0.5) & (self.radii <= self.rperp[data_idx]+0.5))
		model_pts = self.model[idr]
		#print len(model_pts)
		if self.label[data_idx] == 'u':
			#print 'in upper'
			sigma = np.sqrt((0.3**2.0)+(self.model_width)**2.0)
               		probs = self.upper_limit(model_pts,self.coldens[data_idx],sigma)
			id0 = np.where(probs ==0.0)
			probs[id0] = 1e-18
			#print probs
		elif self.label[data_idx] == 'l':
			#print 'in lower'
			sigma = np.sqrt((0.3**2.0)+(self.model_width)**2.0)
			probs = self.lower_limit(model_pts,self.coldens[data_idx],sigma)
                        id0 = np.where(probs ==0.0)
			probs[id0] = 1e-18
			#print probs
		elif self.label[data_idx] == 'n':
			#print 'in normal'
			sigma = np.sqrt((self.err[data_idx])**2.0+(self.model_width)**2.0)
			probs = self.detection(model_pts,self.coldens[data_idx],sigma)
                        id0 = np.where(probs ==0.0)
			probs[id0] = 1e-18
			#print probs
		else:
                       print 'SOMETHING IS GOING WRONG....'
		return np.sum(np.log(probs))/len(probs)

	def total_probability(self):
		for data_idx in range(len(self.idx)):
			self.likelihood = self.likelihood + self.model_probability(data_idx)
		return
		
	
	

