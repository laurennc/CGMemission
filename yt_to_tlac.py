import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from yt.mods import *
import tlac_grid_create_hdf5 as tg
import tlac_analysis as ta

fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
#val, pos = pf.h.find_max('Density')
fields = ['H_NumberDensity','HI_NumberDensity','Temperature','Metallicity','x-velocity','y-velocity','z-velocity','x','y','z']
center = [ 0.39871597,  0.46913528,  0.46808243]
dims = [320,320,320]

level = 8 #pf.index.max_level
ncells = pf.domain_dimensions * pf.refine_by**level
full_boxsize = (pf.domain_width * ( 1. + 1. / ncells ))
cell_length = (full_boxsize[0]/ncells[0])  #0.90829222634295292 kpc

width = dims[0]*cell_length ##290.65351242974492 kpc
axis = fix_axis('x')
center = np.array(center)
LE,RE = center.copy(),center.copy()
LE[axis] -= width/2.0
RE[axis] += width/2.0

area_axes = [0,1,2]
i = area_axes.index(axis)
del area_axes[i]
LE[area_axes] -= width/2.0
RE[area_axes] += width/2.0

cg = pf.h.covering_grid(level, left_edge=LE, dims = dims)

region_size = (cg.right_edge-cg.left_edge)*pf['cm']
pos_coords = [np.array(cg['x'][:,0,0]),np.array(cg['y'][0,:,0]),np.array(cg['z'][0,0,:])]

vx = cg['x-velocity']
vy = cg['y-velocity']
vz = cg['z-velocity']
velocities = np.zeros(np.concatenate( (np.array(vx.shape),[3])))
velocities[:,:,:,0] = vx
velocities[:,:,:,1] = vy
velocities[:,:,:,2] = vz

#######INITIATE FILE############
output_file = 'grid_lauren.hdf5'
f = tg.open_file(output_file)
tg.write_header(f)
tg.write_grid_data(f,1,dims,region_size,0.0)

####ADD FIELDS OF INTEREST#########
##LOG OPTION
#tg.write_cell_data(f,pos_coords,None,np.log10(cg['Temperature']),np.log10(cg['HI_NumberDensity']),None,None,np.log10(cg['H_NumberDensity']),np.log10(cg['Metallicity']))
##NON LOG OPTION
tg.write_cell_data(f,pos_coords,None,cg['Temperature'],cg['HI_NumberDensity'],None,velocities,cg['H_NumberDensity'],cg['Metallicity'])

tg.write_cloud_data(f)
tg.close_file(f)

# Check
#print "Checking data..."
#ta.grid_summary( ta.grid_load(output_file) )
#print "...done!"
