from lauren import *

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 108.0/pf['kpc']

field = ['HAlpha_Emissivity2']

pp = ProjectionPlot(pf,'x',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'])

pp.set_cmap('all','spectral')
pp.save()
