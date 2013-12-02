from lauren import *

#I want to use this file to be able to generalize all my mock observation code!
def area_effective(diameter,efficiency):
        return efficiency*np.pi*0.25*(diameter*1.0e2)**2.0

def plate_scale(focal_length):
        print 'focal length must be in mm!'
        return 206265./focal_length

def pixel_size(frb_res,galaxy_dist):
        print 'give both in pc!'
        return frb_res/galaxy_dist

def pixel_counts(Ipix,Aeff,Opix,Texp):
        #return Ipix*(4.*np.pi/1e6)*Aeff*Opix*Texp
	return Ipix*(1e6/(4.*np.pi))*Aeff*Opix*Texp

def sky_background(Isky,Aeff,Opix,Texp,delLambda):
	return Isky*(4.*np.pi/1.e6)*Aeff*Opix*Texp*delLambda

def mock_observation(frb,resolution,gal_dist,texp,telescope_diameter,efficiency,delLambda):
	Aeff = area_effective(telescope_diameter,efficiency)
	print 'Aeff is ', Aeff
	Opix = pixel_size(resolution,gal_dist)
	Opix = Opix**2.0
	print 'Opix is ',Opix
	image = pixel_counts(frb,Aeff,Opix,texp)
#	image = image + np.sqrt(image)
	return image

#def plot_observation(image,filename):
	
