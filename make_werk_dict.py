import numpy as np
import cPickle

filename = "werk13_coldens.dat"
fileout = "werk_coldens_data.cpkl"

f = open(filename)
lines = f.readlines()[41:]

ID = []
z = []
Rperp = []
ion = []
wavelength = []
vmin = []
vmax = []
l_Wrest = []
Wrest = []
e_Wrest = []
l_logN = []
logN = []
e_logN = []
Flg = []
l_logNA = []
logNA = []
e_logNA = []

i = 0
for line in lines:
	#print i
	ID = np.append(ID,line[0:17])
	z = np.append(z,float(line[18:24]))
	Rperp = np.append(Rperp,float(line[25:28]))
	#want to strip the end of these probs
	ion = np.append(ion,line[29:35].rstrip().replace(" ",""))
	wavelength = np.append(wavelength,int(line[36:40]))
	vmin = np.append(vmin,float(line[41:45]))
	vmax = np.append(vmax,float(line[46:50]))
	if line[51].isspace():
		l_Wrest = np.append(l_Wrest,'n')
		e_Wrest = np.append(e_Wrest,float(line[58:61]))
	else:
		if line[51] == '<':
			l_Wrest = np.append(l_Wrest,'u')
		if line[51] == '>':
			l_Wrest = np.append(l_Wrest,'l')
		e_Wrest = np.append(e_Wrest,0.0)
	Wrest = np.append(Wrest,float(line[53:57]))

	if line[62].isspace():
		l_logN = np.append(l_logN,'n')
		if line[70:74].isspace():
			e_logN = np.append(e_logN,0.0)
		else:
			e_logN = np.append(e_logN,float(line[70:74]))
	else:
		if line[62] == '<':
			l_logN = np.append(l_logN,'u')
		if line[62] == '>':
			l_logN = np.append(l_logN,'l')
		e_logN = np.append(e_logN,0.0)
	
	if line[64:69].isspace():
		logN = np.append(logN,0.0)
	else:
		logN = np.append(logN,float(line[64:69]))
	
	Flg = np.append(Flg,int(line[75:77]))
	
	if line[78].isspace():
		l_logNA = np.append(l_logNA,'n')
		if line[86:90].isspace():
                        e_logNA = np.append(e_logNA,0.0)
                else:
                        e_logNA = np.append(e_logNA,float(line[86:90]))
        else:
                if line[78] == '<':
                        l_logNA = np.append(l_logNA,'u')
                if line[78] == '>':
                        l_logNA = np.append(l_logNA,'l')
                e_logNA = np.append(e_logNA,0.0)

	if line[80:85].isspace():
                logNA = np.append(logNA,0.0)
        else:
                logNA = np.append(logNA,float(line[80:85]))
	i = i + 1

data = {}
data['ID'] = ID
data['z'] = z
data['Rperp'] = Rperp
data['ion'] = ion
data['wavelength'] = wavelength
data['vmin'] = vmin
data['vmax'] = vmax
data['l_Wrest'] = l_Wrest
data['Wrest'] = Wrest
data['e_Wrest'] = e_Wrest
data['l_logN'] = l_logN
data['logN'] = logN
data['e_logN'] = e_logN
data['Flg'] = Flg
data['l_logNA'] = l_logNA
data['logNA'] = logNA
data['e_logNA'] = e_logNA

cPickle.dump(data,open(fileout,'wb'),protocol=-1)

#return data


