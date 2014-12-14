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
		self.model = cPickle.load(open(self.model_name,'rb'))
		self.galaxy = galaxy
		self.radii = radii
		self.nbins = galaxy.nbins
		self.ndraws = galaxy.ndraws	

		#So the functions run without having to edit them
		self.l,self.u,self.model_width = 10.**galaxy.l,10.**galaxy.u,galaxy.model_width
		self.coldens = 10.**galaxy.coldens[idx]
		self.rperp = galaxy.rperp[idx]		
		self.label = galaxy.label[idx]
		self.err = 10.**galaxy.err[idx]
		#Finally, create a variable that will hold the total likelihood of the model for this line!
		self.likelihood = 0.0
	

	#l and u are the upper and lower limits of the model range that we are currently considering
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
		idr = np.where((self.radii >= self.rperp[data_idx]-2.5) & (self.radii <= self.rperp[data_idx]+2.5))
		model_pts = self.model[idr]
		#Here I'm just binning the data and not adding in any measurement error or smoothing error
		#perhaps one way to not add gaussian is to add +x for real bin and +y for a given smoothing to whatever number of bins we want
		#hist_only,edges_only = np.histogram(model_pts,bins=self.nbins,range=(self.l,self.u),density=True)
	
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
		hist_org,edges_org = np.histogram(np.log10(model_pts),bins=np.linspace(np.log10(self.l),np.log10(self.u),num=self.nbins),density=True)

		if self.label[data_idx] == 'u':
			hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx])
			rect1 = matplotlib.patches.Rectangle((edges[0],0), np.log10(self.coldens[data_idx])-edges[0], np.max(hist), color='DeepPink',alpha=0.15)
			ax.add_patch(rect1)
		if self.label[data_idx] == 'l':
			hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx])
			rect1 = matplotlib.patches.Rectangle((np.log10(self.coldens[data_idx]),0), edges[-1]-np.log10(self.coldens[data_idx]), np.max(hist), color='DeepPink',alpha=0.15)
			ax.add_patch(rect1)
		if self.label[data_idx] == 'n':
			hist,edges,npts = self.create_histogram(model_pts,self.coldens[data_idx]*(self.err[data_idx]-1.))
			boxwidth = np.log10(2.*self.coldens[data_idx]*(self.err[data_idx]-1.))
			#rect1 = matplotlib.patches.Rectangle((np.log10(self.coldens[data_idx]-0.5*10.**boxwidth),0),boxwidth,np.max(hist),color='DeepPink',alpha=0.15)
			#ax.add_patch(rect1)
	
		rect2 = matplotlib.patches.Rectangle((edges[0],0), 0.3, hist[0], color='DodgerBlue',alpha=0.15)
                ax.add_patch(rect2)

		ax.plot(edges_org[1:],hist_org,color='Gray')
		ax.plot(edges[1:],hist,color='Indigo',linestyle='--')
                ax.axvline(x=np.log10(self.coldens[data_idx]),color='DeepPink',linewidth=2.2,linestyle='-.')

		prob = self.model_probability(data_idx)	
		prob = np.exp(prob)

		ax.set_xlabel('Values')
		#ax.set_ylabel('Density')
		ax.set_title('p='+"{0:.3f}".format(prob))
		
		return
