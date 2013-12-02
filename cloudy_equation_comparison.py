from lauren import *
import matplotlib.pyplot as plt

hden_n_bins,hden_min, hden_max = 15,-6,1
T_n_bins, T_min, T_max = 51,3,8
patt =  "/hpc/astrostats/astro/users/lnc2115/codes/cloudy_yt/yuan/yuan_hemissivity_run%i.dat"
from scipy import interpolate
Hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
T=numpy.linspace(T_min,T_max, T_n_bins)
table = numpy.zeros((hden_n_bins, T_n_bins))
for i in range(hden_n_bins):
       table[i,:]=[float(l.split()[-1]) for l in open(patt%(i+1)) if l[0] != "#"]

t,hden = np.meshgrid(T,Hden)

#for comparison with the equation predictions...
table2 = np.log10((10**table)*((10**hden)**2.0)*1.87e-12)
equation = np.log10(HAlpha_Emission_Arc(t,hden))

#for the ratio of the two should just be able to take the difference since they're in the log
ratio = table2-equation

plt.imshow(table2,extent=(3,8,1,-6),interpolation='none')
plt.title('Cloudy')
plt.xlabel("Temperature")
plt.ylabel('H Density')
plt.colorbar()
plt.savefig('halphaCloudy.png')
plt.close()

plt.imshow(equation,extent=(3,8,1,-6),interpolation='none')
plt.title('Equation')
plt.xlabel("Temperature")
plt.ylabel('H Density')
plt.colorbar()
plt.savefig('halphaEquation.png')
plt.close()

plt.imshow(ratio,extent=(3,8,1,-6),interpolation='none')
plt.title('Ratio')
plt.xlabel("Temperature")
plt.ylabel('H Density')
plt.colorbar()
plt.savefig('halphaRatio.png')
plt.close()


