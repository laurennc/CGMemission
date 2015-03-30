import numpy as np
import matplotlib.pyplot as plt
import cPickle
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from astropy.io import fits
import os 

class ImageProducer(object):
	"""ImageProducer receives one of the frb image arrays that I've produced.
	   It then convolves the image with a Gaussian of determinable size in order
	   to build the wavelength dimension. It will also scale the image for the 
	   correct redshift of interest.

	   The output has units of photons cm^-2 s^-1 sr^-1 A^-1 """

	def __init__(self, frbarray,lambda_center,lambda_range,lambda_res,lambda_line,step_width,z_observed):
	#"""The inputs are as follows:
	#   frbarray = output projection of the simulation
        #   lambda_center = central wavelength of the desired wavelength window
        #   lambda_range = +/- this value covers the full wavelength range desired
        #   lambda_res = resolution of spectral dimension
        #   lambda_line = wavelength of the line being observed
        #   step_width_z = how wide should the step function be that approximates the simulation signal in resolution units
        #   z_observed = the redshift of the observation for scaling the flux
#
        #   The output variables are as follows:
        #   img = array that will contain the convolved imaged. Initialized as 0  """

		self.frb = frbarray
		self.wavelengths = np.arange(lambda_center-lambda_range,lambda_center+lambda_range+lambda_res,lambda_res)[:-1]
		#self.wavelengths = np.arange(0.204-0.0012,0.204+0.0012+3.4e-6,3.4e-6)
		
		#build the step function that will lead to the convolver shape
		self.step_function = np.zeros(len(self.wavelengths))
		idx = (np.abs(self.wavelengths-lambda_line)).argmin()
		self.step_function[idx-step_width/2.:idx+step_width/2.+1] = 1.0/(step_width*lambda_res)

		self.z_observed = z_observed

		self.img = 0

	def convolve_image(self,stddev,scale=True):
       # """The inputs are as follows:
       #    stddev = standard deviation of the Gaussian the step function will be convolved with
       #    scale = True -- will scale the flux by the redshift unless set to False
#
       #    The output variables are as follows:
       #    Nothing is returned but the object property, img, is updated to hold the convolved image
       #  """		

##CURRENTLY BUILDING IMAGE AS LAMBA,X,Y
		##Because writing the fits file transposes the x and z axes for some reason
		self.img = np.zeros((len(self.wavelengths),len(self.frb[0,:]),len(self.frb[1,:])))


		gauss = Gaussian1DKernel(stddev=stddev,x_size=len(self.wavelengths))
		zconvolve = convolve(self.step_function,gauss.array,boundary='extend')
		
		for i in range(len(self.frb[0,:])):
			for j in range(len(self.frb[1,:])):
				if scale:
					self.img[:,i,j] = zconvolve*(self.frb[i,j]/((1.+self.z_observed)**4.0))
				else:
					self.img[:,i,j] = zconvolve*self.frb[i,j]

		return

	def save_as_fits(self,outputname):
	#Save the convolution as a fits file
		try:
			os.remove(outputname)
		except OSError:
			pass
		hdu = fits.PrimaryHDU(self.img)
		hdu.writeto(outputname)		
		return

	
