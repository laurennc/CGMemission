import numpy as np
import cPickle

filename = "werk13_galaxy_properties.dat"
fileout = "werk_galaxy_properties.cpkl"

f = open(filename)
lines = f.readlines()[5:]

ID = []
z = []
Rperp = []
mStar = []
lumin = []
urColor = []
urColorTag = []
SFR = []
SFRTag = []
abund = []

i = 0
for line in lines:
	vals = line.rstrip("\t\n").split("\t")
	ID = np.append(ID,vals[0])
	z = np.append(ID,float(vals[1]))
	Rperp = np.append(Rperp,float(vals[2]))
	mStar = np.append(mStar,float(vals[3]))
	lumin = np.append(lumin,float(vals[4]))


	if vals[5].startswith('>'):
		wanted = vals[5].lstrip('>')
		urColor = np.append(urColor,float(wanted))
		urColorTag = np.append(urColorTag,True)
	else:
		urColor = np.append(urColor,float(vals[5]))
		urColorTag = np.append(urColorTag,False)

        if vals[6].startswith('<'):
                wanted = vals[6].lstrip('<')
                SFR = np.append(SFR,float(wanted))
                SFRTag = np.append(SFRTag,True)
        else:
                SFR = np.append(SFR,float(vals[6]))
                SFRTag = np.append(SFRTag,False)

        if vals[7] == "****":
                urColor = np.append(abund,0.0)
        else:
                urColor = np.append(abund,float(vals[7]))


data = {}
data['ID'] = ID
data['z'] = z
data['Rperp'] = Rperp
data['mStar'] = mStar
data['lumin'] = lumin
data['urColor'] = urColor
data['urColorTag'] = urColorTag
data['SFR'] = SFR
data['SFRTag'] = SFRTag
data['abund'] = abund

cPickle.dump(data,open(fileout,'wb'),protocol=-1)

