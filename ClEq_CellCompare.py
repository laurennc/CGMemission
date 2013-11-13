import matplotlib as mpl
mpl.use('agg')
from lauren import *
import matplotlib.pyplot as plt
import triangle as triangle

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]
fields = ["HAlpha_Emissivity", "HAlphaEmissionArc","Temperature"] 

rad = 108.0/pf['kpc']
data = pf.h.sphere(pos,rad)

datain = np.zeros((len(data['Temperature']),3))
for i in range(3):
	datain[:,i] = np.log10(data[fields[i]])

fileout='cloudy_eq_test.png'
triangle.corner(datain,labels=['Cloudy HAlpha','Equation HAlpha','Temperature']).savefig(fileout)
plt.close()

ratio = np.log10(data['HAlphaEmissionArc']/data["HAlpha_Emissivity"])
triangle.hist2d(np.log10(data['Temperature']),ratio)
plt.xlabel('Temperature')
plt.ylabel('Ratio Equation/Cloudy')
plt.savefig('ratio_equation_cloudy_temp.png')
plt.close()

idx = np.where(np.log10(data['Temperature']) > 4.5)
triangle.hist2d(np.log10(data['Temperature'][idx]),ratio[idx])
plt.xlabel('Temperature')
plt.ylabel('Ratio Equation/Cloudy')
plt.savefig('ratio_equation_cloudy_tempnocold.png')
plt.close()

x = np.linspace(-48,33,15)
triangle.hist2d(datain[:,0],datain[:,1])
plt.xlabel('Cloudy HAlpha')
plt.ylabel('Equation HAlpha')
plt.plot(x,x,'r-',linewidth='3')
plt.plot(x,(x+2),'b--',linewidth='3')
plt.plot(x,(x+3),'g--',linewidth='3')
plt.axis([-48,-33,-48,-33])
plt.axes().set_aspect('equal')
plt.savefig('equation_vs_cloudy.png')
plt.close()


