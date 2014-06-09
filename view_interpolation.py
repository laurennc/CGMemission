from lauren import *
from scipy import interpolate

hden1,T1,table1 = make_Cloudy_table(15)
hden1,T1 = np.meshgrid(hden1,T1)
pts1 = np.array((hden1.ravel(),T1.ravel())).T
sr1 = table1.T.ravel()
bl1 = interpolate.LinearNDInterpolator(pts1,sr1)

tnow = np.arange(3,8,0.005)
hnow = np.zeros(len(tnow))+1

iax = 331
fig = plt.figure(figsize=(13.5,13.5))

while hnow.min() > -7:
	z1 = bl1(hnow,tnow)
	ax1 = fig.add_subplot(iax)
	ax1.plot(tnow,z1)
	ax1.set_xlabel('Temperature')
	ax1.set_ylabel('Emissivity')
	ax1.text(3.5,-21,str(hnow.min()))
	iax = iax + 1
	hnow = hnow - 1.0

plt.savefig('test_interp.png')

