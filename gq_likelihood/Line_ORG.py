#Class designed to estimate the probability of a certain ionizing model given the data of Werk et al. 2013
#Follows Lemonias et al. 2012 in calculating the probabilities of the data, upper, and lower limits
#It can also handle adding uncertainty to the data to represent physical scatter if desired
#Created by: Lauren
#Created on: 11//4/2014
#Edited on: 12/4/2014 -- can now store individual likelihoods
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special

class Line:
	def __init__(self,ion,model_name,l,u,model_width,nbins,ndraws):
		self.ion = ion
		self.model_name = model_name
		self.nbins,self.ndraws = nbins, ndraws

		#Load the data from Werk et al. 2013
		data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
		idx = np.where((data['ion'] == ion) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30))[0]# & (data['l_logNA']=='n'))[0]
		
		#Keep the data that I want for this line
		self.coldens = 10.**data['logNA'][idx]
		self.rperp = data['Rperp'][idx]
	        self.label = data['l_logNA'][idx]
	        self.err = 10.**data['e_logNA'][idx]
		
		#Lower and upper limits for the normalization
		self.l,self.u,self.model_width = 10.**l,10.**u,10.**model_width

		#Now loading in information for the model
		self.model = cPickle.load(open(self.model_name,'rb'))
		xL = np.linspace(-160,160,320)
		xL,yL = np.meshgrid(xL,xL)
		self.radii = abs(xL+1j*yL)
		
		#Finally, create a variable that will hold the total likelihood of the model for this line!
		self.total_likelihood = 0.0
		self.likelihoods = []
#
	def upper_limit(self,bin_idx,hist,del_x):
                return sum(hist[:bin_idx+1]*del_x)

        def lower_limit(self,bin_idx,hist,del_x):
                return sum(hist[bin_idx:]*del_x)

	def create_histogram(self,model_pts,err):
                hist_pts = []
                for model_pt in model_pts:
                        draws = np.random.normal(model_pt,err,self.ndraws)
                        hist_pts = np.append(hist_pts,draws)

                #Making sure all the points below the lower limit get counted
                id_l = np.where(hist_pts < self.l)[0]
                hist_pts[id_l] = self.l
                #and now for the upper limit
                id_u = np.where(hist_pts > self.u)[0]
                hist_pts[id_u] = self.u

                hist,edges = np.histogram(np.log10(hist_pts),bins=self.nbins,range=(np.log10(self.l),np.log10(self.u)),density=True)
                return hist, edges, len(hist_pts)


	def model_probability(self,data_idx):
		#print self.rperp[data_idx]
		idr = np.where((self.radii >= self.rperp[data_idx]-2.5) & (self.radii <= self.rperp[data_idx]+2.5))
		model_pts = self.model[idr]
		
		if self.label[data_idx] == 'n':
                        hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx]*(self.err[data_idx]-1.))
                        del_x = edges[1]-edges[0]
                        bin_want = int((np.log10(self.coldens[data_idx])-np.log10(self.l))/del_x)
                        prob = hist[bin_want]
                elif self.label[data_idx] == 'l':
                        hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx])
                        del_x = edges[1]-edges[0]
                        bin_want = int((np.log10(self.coldens[data_idx])-np.log10(self.l))/del_x)
                        prob = self.lower_limit(bin_want,hist,del_x)
                elif self.label[data_idx] == 'u':
                        hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx])
                        del_x = edges[1]-edges[0]
                        bin_want = int((np.log10(self.coldens[data_idx])-np.log10(self.l))/del_x)
                        prob = self.upper_limit(bin_want,hist,del_x)
                else:
                        print 'WARNING: SOMETHING IS WRONG WITH THE DATA'
                if prob == 0.0:
                        prob = 1./npts

		return np.log(prob)

	def total_probability(self):
		for data_idx in range(len(self.coldens)):
			temp_val = self.model_probability(data_idx)
			self.likelihoods = np.append(self.likelihoods,temp_val)
			self.total_likelihood = self.total_likelihood + temp_val
		self.total_likelihood = self.total_likelihood - np.log(len(self.coldens))
		return
		
	
	

