import numpy as np

def load_emis_table():
	filename = 'werk.emis'
	
	lines = open(filename,'r').readlines()
	emis = {}
	count = 0
	for l in lines:
		if l.startswith('depth'):
			keys = l.split("\t")
			for i in np.arange(len(keys)):
			    keys[i] = keys[i].rstrip('\n')
			    emis[keys[i]] = []
			nlines = (len(lines)-1)/2
			data = np.empty((nlines,len(keys)))
		elif l.startswith('###########################'):
			pass
		else:
			data[count,:] = map(float,l.split())
			count = count + 1

	for i in np.arange(len(keys)):
		emis[keys[i]] = data[:,i]
	return emis

def return_heat_map_array(emis,key):
	metal_steps = [-3.,-2.5,-2.,-1.5,-1.,-0.5,-0.]
	dens_steps = [-5.,-4.,-3.,-2.,-1.]

	data = np.empty((len(dens_steps),len(metal_steps)))
	low = 0.
	high = len(metal_steps)
	i = 0
	while high <= len(emis[key])+1:
		data[i,:] = emis[key][low:high]
		low = low + len(metal_steps)
		high = high + len(metal_steps)
		i = i + 1
	return data

