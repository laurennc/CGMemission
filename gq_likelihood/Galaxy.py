import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from Line import *

class Galaxy(object):
	def __init__(self,ID,l,u,model_width,nbins,ndraws):
		self.ID,self.u,self.l,self.model_width,self.nbins,self.ndraws = ID, u, l, model_width, nbins, ndraws
	
		data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
#Choosing the data that I want for this ion, within the radius range David and I decided upon
                #logNA > 0.1 ---- double check what this is ---- right now copying from other code
        
		####DETECTIONS ONLY####
	        idx = np.where((data['ID'] == ID) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30) & (data['l_logNA']=='n'))[0]
		#idx = np.where((data['ID'] == ID) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30))[0]
	

		self.ions = data['ion'][idx]
		self.coldens = data['logNA'][idx]
		self.rperp = data['Rperp'][idx]
		self.label = data['l_logNA'][idx]
		self.err = data['e_logNA'][idx]
		
		self.models = []
		self.best_model = -1
		self.best_q,self.best_g = -1, -1

	def add_model(self,g,q,model_name,galaxy):
		mnow = Model(g,q,model_name,galaxy)
		self.models = np.append(self.models,mnow)
		return

	def find_best_model(self):
		likehold = -99999999999
		likeidx = -1
		for i in range(len(self.models)):
			if self.models[i].total_likelihood > likehold:
				likehold = self.models[i].total_likelihood
				likeidx = i
				self.best_q = self.models[i].q
				self.best_g = self.models[i].g
		self.best_model = likeidx
		return

	def plot_pdfs_oneline_allmodels(self,ion,fileout):
		data_idx = np.where(self.ions ==ion)[0][0]
		fig,axs = plt.subplots(4,4)
		fig.set_size_inches(13.5,13.5)
		ax2 = axs.flat
		i = 0
		for model in self.models:
			for line in model.lines:
				if line.ion == ion:
					line.plot_pdf(data_idx,ax2[i])
					ax2[i].set_ylabel('g:'+str(model.g)+' q:'+str(model.q))
			i = i + 1

		plt.tight_layout()
		plt.savefig(fileout)
		plt.close()
		return
	

class Model():
        def __init__(self,g,q,model_name,galaxy): #,k) for tanh
                self.lines = []
                self.likelihood = 0.0
		self.galaxy = galaxy
		self.nlines = 0	
		self.model_name = model_name
		self.g, self.q = g, q		

		xL = np.linspace(-160,160,320)
                xL,yL = np.meshgrid(xL,xL)
                self.radii = abs(xL+1j*yL)
		
		self.total_likelihood = 0.0

	def line_indices(self,ion):
		return np.where(self.galaxy.ions == ion)[0]

	def add_line(self,ion,idx):
		lnow = Line(ion,idx,self.model_name,self.galaxy,self.radii)
		lnow.total_probability()
		self.lines = np.append(self.lines,lnow)
		self.nlines = self.nlines + 1
		return	

	def total_model_likelihood(self):
		for line in self.lines:
			self.nlines = self.nlines + 1
			self.total_likelihood = self.total_likelihood + line.likelihood
		self.total_likelihood - np.log10(float(self.nlines))
		return




