from lauren import *
import matplotlib

#fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 108.0/pf['kpc']

data_source = pf.h.sphere(pos,rad)

#field = ['HI_NumberDensity']
#field = ['Emission_CIV','Emission_OVI']#,'HAlpha_Emissivity_R']
#field = ['Emission_CIV','Emission_CIV_Scaled','Emission_CIV_Scaled_ncut']
#field = ['Emission_OVI','Emission_OVI_Scaled','Emission_OVI_Scaled_ncut']
field = ['Temperature','Metallicity','SolarMetals']

font = {'weight':'bold','size':22,'family':'sans-serif'}
matplotlib.rc('font',**font)

#pp = ProjectionPlot(pf,'y',field,center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)
pp = ProjectionPlot(pf,'x',field,center=pos,width=(108.,'kpc'),weight_field="H_NumberDensity",axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
pp.set_font(font)
pp.set_cmap('all','spectral')
#pp.set_zlim('all',10**-5.0,10**2.0)
#pp.set_zlim('all',10**14.0,10**23.5)
#pp.plots['SolarMetals'].cb.set_ticks([0.1,0.3,0.5,0.7,1.0])
#pp.plots['SolarMetals'].cb.set_ticklabels([0.1,0.3,0.5,0.7,1.0])
pp.save()

#pp = ProjectionPlot(pf,'x','Metallicity',center=pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)#,fontweight='bold')
#pp = SlicePlot(pf,'x','Metallicity',pos,width=(108.,'kpc'),axes_unit=['kpc','kpc'],fontsize=22.)
#pp.set_font(font)
#pp.set_cmap('all','spectral')
#pp.set_zlim('all',10**-3.0,10**1.0)
#pp.save()


