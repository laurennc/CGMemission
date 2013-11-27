from lauren import *
import matplotlib

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 108.0/pf['kpc']

#field = ['HI_NumberDensity']
field = ['Emission_CIV','Emission_OVI']#,'HAlpha_Emissivity_R']

font = {'weight':'bold','size':22,'family':'sans-serif'}
matplotlib.rc('font',**font)

#pp = ProjectionPlot(pf,'y',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)
pp = ProjectionPlot(pf,'x',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
pp.set_font(font)
pp.set_cmap('all','spectral')
pp.set_zlim('all',10**-5.0,10**2.0)
#pp.set_zlim('all',10**14.0,10**23.5)
pp.save()
