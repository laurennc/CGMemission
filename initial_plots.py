from yt.mods import *
import numpy as np

#fn="/Users/laurennc/data/RyanSims/r0058_l10/redshift0058"
fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#what types of plots do I want to make?  Density projections, and HI projections that are all the same size, color maps, and limits as Ximena's data!

#Need to find the position of the maximum density!
#val, pos = pf.h.find_max('Density')
#Saving this output for the future:
val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

rad = 108.0/pf['kpc']
data = pf.h.sphere(pos,rad)

fields = ["Density","Temperature","Electron_Density","HAlphaEmissionArc"]

for dim in "xyz":
	pp = ProjectionPlot(pf,dim,fields,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'])
	pp.set_cmap('all','spectral')
	pp.save()

for dim in "xyz":
	pp = ProjectionPlot(pf,dim,'HI_Number_Density',center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'])
	pp.set_cmap('all','spectral')
	pp.set_zlim(10**14.0,10**23.5)
	pp.save()


#L = [-0.5934,1,1]
#L = [0.4066,1,1]
#L = [0.9774,1,1]
#prj = OffAxisProjectionPlot(pf,L,'Density',width=rad,center=pos)
#prj.save('testing_face')

#pp.set_zlim(10**14.0,10**23.5)
#pp.set_cmap('spectral')

