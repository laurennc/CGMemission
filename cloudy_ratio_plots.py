import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from lauren import make_ion_table
from lauren import make_Cloudy_table
from skimage import measure

## pour tous
hden_yes, T_yes, table_HA = make_Cloudy_table(2)
hden_yes, T_yes = np.meshgrid(hden_yes,T_yes)
pts = np.array((hden_yes.ravel(),T_yes.ravel())).T


## CIII Emission
hden1,T1,table_CIII_977 = make_Cloudy_table(7)
sr_CIII_977 = table_CIII_977.T.ravel()
CIII_Emission = table_CIII_977

## CIV  Emission
hden1, T1, table_CIV_1 = make_Cloudy_table(3)
hden1, T1, table_CIV_2 = make_Cloudy_table(4)
sr_CIV_1 = table_CIV_1.T.ravel()
sr_CIV_2 = table_CIV_2.T.ravel()
CIV_Emission = 10.**(table_CIV_1) + 10.**(table_CIV_2)
CIV_Emission = np.log10(CIV_Emission)

## OVI  Emission
hden1, T1, table_OVI_1 = make_Cloudy_table(5)
hden1, T1, table_OVI_2 = make_Cloudy_table(6)
sr_OVI_1 = table_OVI_1.T.ravel()
sr_OVI_2 = table_OVI_2.T.ravel()
OVI_Emission = 10.**(table_OVI_1) + 10.**(table_OVI_2)
OVI_Emission = np.log10(OVI_Emission)

## NV   Emission
## Need to find the right line for NV and then run all of the Cloudy models....


## first just plotting them!
plt.imshow(CIII_Emission,origin='lower',extent=[3,8,-6,2],vmin=-50,interpolation='nearest')
plt.colorbar()
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.savefig('LineRatioPlots/CIII_Emission_g1q1.png')
plt.close()

plt.imshow(CIV_Emission,origin='lower',extent=[3,8,-6,2],vmin=-50,interpolation='nearest')
plt.colorbar()
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.savefig('LineRatioPlots/CIV_Emission_g1q1.png')
plt.close()

plt.imshow(OVI_Emission,origin='lower',extent=[3,8,-6,2],vmin=-50,interpolation='nearest')
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.colorbar()
plt.savefig('LineRatioPlots/OVI_Emission_g1q1.png')
plt.close()

## Now plotting the ratios that I want!
plt.imshow(CIII_Emission-CIV_Emission,origin='lower',extent=[3,8,-6,2],vmin=-25,vmax=25,interpolation='nearest')
plt.colorbar()
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.savefig('LineRatioPlots/CIII_CIV_Emission_g1q1.png')
plt.close()

plt.imshow(CIII_Emission-OVI_Emission,origin='lower',extent=[3,8,-6,2],vmin=-25,vmax=25,interpolation='nearest')
plt.colorbar()
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.savefig('LineRatioPlots/CIII_OVI_Emission_g1q1.png')
plt.close()

plt.imshow(CIV_Emission-OVI_Emission,origin='lower',extent=[3,8,-6,2],vmin=-25,vmax=25,interpolation='nearest')
plt.colorbar()
plt.xlabel('Temperature')
plt.ylabel('H Number Density')
plt.savefig('LineRatioPlots/CIV_OVI_Emission_g1q1.png')
plt.close()



