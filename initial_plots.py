from lauren import *
#from yt.mods import *
#import numpy as np
import matplotlib.pyplot as plt
import triangle as triangle

#fn="/Users/laurennc/data/RyanSims/r0058_l10/redshift0058"
fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#Need to find the position of the maximum density!
#val, pos = pf.h.find_max('Density')
#Saving this output for the future:
val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

rad = 108.0/pf['kpc']
data = pf.h.sphere(pos,rad)

#First let's look at HAlpha vs HI Col Dens on a cell by cell basis
halphaArc = np.log10(data['HAlphaEmissionArc'])
em = np.log10(data['EmissionMeasurePC'])
hIcoldens = np.log10(data['HI_Number_Density']*data['dx'])
datain = np.zeros((len(hIcoldens),2))
datain[:,1] = em
datain[:,0] = hIcoldens
#datain[:,1] = halphaArc
#fileout = 'hIcoldens_vs_halphaArc.png'
#triangle.corner(datain,labels=['log(HI Column Density)','log(H Alpha Arc)']).savefig(fileout)
fileout = 'hIcoldens_vs_emPC.png'
triangle.corner(datain,labels=['log(HI Column Density) (cm^-2)','log(Emission Measure) (cm^-6 pc)']).savefig(fileout)


fields = ["Density","EmissionMeasure","Electron_Density_Squared","HAlphaEmissionSr","Temperature","Electron_Density","HAlphaEmissionArc","HI_Number_Density"]
fields = ["EmissionMeasureCM","EmissionMeasurePC","HI_Number_Density"]
fields = ['HAlphaEmissionRal','HAlphaEmissionSr','HAlphaEmissionArc']
fields = ['EmissionMeasurePC','EmissionMeasureCold','HAlphaEmissionSr','HAlphaEmissionRal','HI_Number_Density','Temperature']

fields = ['HI_Number_Density','HI_Number_Density_DX','EmissionMeasurePC','EmissionMeasureCold','EmissionMeasurePCDX','HAlphaEmissionArc','HAlphaEmissionRal']
proj = pf.h.proj(0,fields)
width = 108./pf['kpc']
res = [1000,1000]
frb = proj.to_frb(width,res,center=pos)


datain = np.zeros((len(frb['HI_Number_Density'].flatten()),4))
datain[:,0] = np.log10(frb['HI_Number_Density'].flatten())
#datain[:,1] = np.log10(frb['HAlphaEmissionArc'].flatten()*2.3582*10**-11.)
datain[:,1] = np.log10(frb['EmissionMeasurePC'].flatten())
datain[:,2] = np.log10(frb['HAlphaEmissionArc'].flatten())
datain[:,3] = np.log10(frb['Temperature'].flatten())
fileout = 'hIcoldens_em_halphaArc_temp_PROJ.png'
triangle.corner(datain,labels=['log(HI Column Density)','log(EM) (cm^-6 pc)','HAlpha Arc','log(Temperature)']).savefig(fileout)

datain = np.zeros((len(frb['Temperature'].flatten()),4))
datain[:,0] = np.log10(frb['Temperature'].flatten())
datain[:,1] = np.log10(frb['EmissionMeasureCold'].flatten())
fileout = 'temp_emcold_PROJ.png'
triangle.corner(datain,labels=['log(Temperature)','log(EM) (cm^-6 pc)']).savefig(fileout)

datain = np.zeros((len(idC),4))
datain[:,0] = np.log10(data['Temperature'][idC])
datain[:,1] = np.log10(data['EmissionMeasureCold'][idC])
dafgadf
fileout = 'hIcoldens_temp_emcold_CELL.png'
triangle.corner(datain,labels=['log(HI Column Density)','log(Temperature)','log(EM) (cm^-6 pc)']).savefig(fileout)


#can really only do it on a cell-by-cell basis since the projection seems so weird for the cold data
datain = np.zeros((len(idC),4))
datain[:,0] = np.log10(data['HI_Number_Density'][idC]*data['dx'][idC]*pf['cm'])
datain[:,1] = np.log10(data['EmissionMeasureCold'][idC])
datain[:,2] = np.log10(data['HAlphaEmissionArc'][idC])
datain[:,3] = np.log10(data['Temperature'][idC])
fileout = 'hIcoldens_emcold_halphaArc_temp_CELL.png'
triangle.corner(datain,labels=['HI Col Dens','EM Cold','HAlpha Arc','Temp']).savefig(fileout)

datain = np.zeros((len(data['Temperature']),4))
datain[:,0] = np.log10(data['HI_Number_Density']*data['dx']*pf['cm'])
datain[:,1] = np.log10(data['EmissionMeasurePC']*data['dx']*pf['pc'])
datain[:,2] = np.log10(data['HAlphaEmissionArc'])
datain[:,3] = np.log10(data['Temperature'])
fileout = 'hIcoldens_em_halphaArc_temp_CELL.png'
triangle.corner(datain,labels=['HI Col Dens','EM (cm^-6 pc)','HAlpha Arc','Temp']).savefig(fileout)




#max in the fbr:
i,j = np.unravel_index(frb['Density'].argmax(),frb['Density'].shape)
frbcenter = [i,j]

x1 = (np.arange(500)+1)*0.108
x2 = -1*x1
x2 = x2[::-1]
x2 = x2[1:]
x = np.concatenate((x2,[0]))
x = np.concatenate((x,x1))

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

