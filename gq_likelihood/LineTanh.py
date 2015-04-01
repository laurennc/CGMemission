#Class designed to estimate the probability of a certain ionizing model given the data of Werk et al. 2013
#Created by: Lauren
#Created on: 10/21/2014
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

class Line:
	def __init__(self,ion,model_name,l,u,k):
		self.ion = ion
		self.model_name = model_name

		#Load the data from Werk et al. 2013
		data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
		#Choosing the data that I want for this ion, within the radius range David and I decided upon
		#logNA > 0.1 ---- double check what this is ---- right now copying from other code
		idx = np.where((data['ion'] == ion) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30))[0]
		
		#Keep the data that I want for this line
		self.coldens = data['logNA'][idx]
		self.rperp = data['Rperp'][idx]
	        self.label = data['l_logNA'][idx]
	        self.err = data['e_logNA'][idx]
		
		#Parameters for the tanh functions
		self.l,self.u,self.k = l,u,k
		#Indexes for upper and lower limits and detections
        	#self.upper = np.where(l_logNA == 'u')[0]
        	#self.lower = np.where(l_logNA == 'l')[0]
        	#self.detect  = np.where(l_logNA == 'n')[0]

		#Now loading in information for the model
		self.model = np.log10(cPickle.load(open(self.model_name,'rb')))
		xL = np.linspace(-160,160,320)
		xL,yL = np.meshgrid(xL,xL)
		self.radii = abs(xL+1j*yL)
		
		#Finally, create a variable that will hold the total likelihood of the model for this line!
		self.likelihood = 0.0
	
	#l and u are the upper and lower limits of the model range that we are currently considering
	def normalize(self,limit_type,data_pt):
		if limit_type == 'upper':
			a = self.u - self.l - (1./self.k)*np.log(np.cosh(self.k*(self.u-data_pt))/np.cosh(self.k*(self.l-data_pt)))
			a = 1./a
		elif limit_type == 'lower':
			a = self.u - self.l + (1./self.k)*np.log(np.cosh(self.k*(self.u-data_pt))/np.cosh(self.k*(self.l-data_pt)))
			a = 1./a
		else:
			a = 0.
		return a

	def upper_limit(self,model_x,a,data_pt):
		return a - a*np.tanh(self.k*(model_x-data_pt))

	def lower_limit(self,model_x,a,data_pt):
		return a + a*np.tanh(self.k*(model_x-data_pt))

	def detection(self,model_x,data_pt,data_err):
		return np.exp(-np.power(model_x-data_pt,2.)/(2.*np.power(data_err,2.)))

#	def upper_probability(model_x,l,a,k,data_pt):
#	return integrate.quad(upper_limit,l,model_x,args=(a,k,data_pt))
#
#	def lower_probability(model_x,u,a,k,data_pt):
#		return integrate.quad(lower_limit,data_pt,u,args=(a,k,data_pt))
#
	def model_probability(self,data_idx):
		#print self.rperp[data_idx]
		idr = np.where((self.radii >= self.rperp[data_idx]-0.05) & (self.radii <= self.rperp[data_idx]+0.05))
		#idr = np.where((self.radii >= self.rperp[data_idx]-0.5) & (self.radii <= self.rperp[data_idx]+0.5))
		model_pts = self.model[idr]
		#print len(model_pts)
		if self.label[data_idx] == 'u':
			#print 'in upper'
			a = self.normalize('upper',self.coldens[data_idx])
               		probs = self.upper_limit(model_pts,a,self.coldens[data_idx])
			id0 = np.where(probs == 0.0)
			probs[id0] = self.detection(model_pts[id0],self.coldens[data_idx],0.2)
			id0 = np.where(probs == 0.0)
                        probs[id0] = 1e-250
			#print probs
		elif self.label[data_idx] == 'l':
			#print 'in lower'
			a = self.normalize('lower',self.coldens[data_idx])
			probs = self.lower_limit(model_pts,a,self.coldens[data_idx])
			id0 = np.where(probs == 0.0)
                        probs[id0] = self.detection(model_pts[id0],self.coldens[data_idx],0.2)
			id0 = np.where(probs == 0.0)
			probs[id0] = 1e-250
			#print probs
		elif self.label[data_idx] == 'n':
			#print 'in normal'
			probs = self.detection(model_pts,self.coldens[data_idx],self.err[data_idx])
			id0 = np.where(probs == 0.0)
			probs[id0] = 1e-250
			#print probs
		else:
                       print 'SOMETHING IS GOING WRONG....'
		return np.sum(np.log(probs))/len(probs)

	def total_probability(self):
		bah = len(self.coldens)
		for data_idx in range(len(self.coldens)):
			#print data_idx
			#print data_idx,self.ion,len(self.coldens)
			self.likelihood = self.likelihood + self.model_probability(data_idx)
	#	print self.likelihood,len(self.coldens)
		self.likelihood = self.likelihood/float(bah)
		return
		
	
	

