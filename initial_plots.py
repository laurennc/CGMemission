from yt.mods import *

fn="/Users/laurennc/data/RyanSims/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#what types of plots do I want to make?  Density projections, and HI projections that are all the same size, color maps, and limits as Ximena's data!

HI_Density

fields = {'Density','HHHH'}
zlim = {"Density": (),
	 "": ()}

#Need to find the position of the maximum density!
val, pos = pf.h.find_max('Density')

rad = 108.0/pf['kpc']

pp = ProjectionPlot(pf,'z',fields,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'])

pp.set_zlim(14.0,23.5)
pp.set_cmap('spectral')

#need to learn how to make off-axis projections.  Practice by making an edge-on and face-on disk for comparison with Ximena!
