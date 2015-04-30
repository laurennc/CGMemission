from ImageProducer import *

model_gqs = ['g1q1']
model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/z1/'
model_mid = '/frbz_6kpc_1Mpc_z1_'
ions = ['OVI']

##PARAMETERS FOR THE IMO###
lambda_center,lambda_range,lambda_res = 0.204,0.0006,3.4e-6 
z_obs = 1
lambda_emit = np.array([0.1032])
#lambda_obs = lambda_emit*(1.+z_obs)
lambda_obs = np.array([0.2040])
stddev = (100./3.e5)*lambda_emit/lambda_res
step_width = 10.#5.

#i=0
for i in range(len(ions)):
	filename = model_beg+model_gqs[0]+model_mid+ions[i]+'.cpkl'
	frb = cPickle.load(open(filename,'rb'))
	
	ip = ImageProducer(frb,lambda_center,lambda_range,lambda_res,lambda_obs[i],step_width,z_obs)
		
	ip.convolve_image(stddev[i])

	outputfile = 'IMOcube_'+ions[i]+'_z'+str(z_obs)+'step'+str(int(step_width))+'.fits'	
	ip.save_as_fits(outputfile)
	
	

