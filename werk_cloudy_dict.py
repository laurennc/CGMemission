import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import cPickle

werk = np.genfromtxt('werk14_ionparams.dat',names=True,dtype=None)

ID, idf = [], []
T,hden = [], []
Lya, Ha = [], []
CIII,CIV,OVI,SiIV = [],[],[],[]
Lyamax,Lyamin = [],[]
Hamax,Hamin = [],[]
CIIImax,CIIImin = [],[]
CIVmax,CIVmin = [],[]
OVImax,OVImin = [],[]
SiIVmax,SiIVmin = [], []
depth_start,depth_end = [],[]
SiIII_1207,SiIII_1883 = [],[]
SiIII_1207max,SiIII_1207min = [],[]
SiIII_1883max,SiIII_1883min = [],[]

for i in range(len(werk['SDSSField'])):
	Lya_temp = []
	Ha_temp = []
	CIII_temp = []
	CIV_temp = []
	OVI_temp = []
	SiIV_temp = []
	SiIII_1207_temp = []
	SiIII_1883_temp = []	

	for j in "abcd":
		filein = 'WerkLoop/'+str(werk['SDSSField'][i])+'_'+str(werk['ID'][i])+'_'+j
		ovr = np.genfromtxt(filein+'.ovr',names=True,dtype=None)
		emis = np.genfromtxt(filein+'.emis',names=True,dtype=None,delimiter='\t')

		ID = np.append(ID,werk['SDSSField'][i]+'_'+werk['ID'][i])
		idf = np.append(idf,j)
		T = np.append(T,ovr['Te'][0])
		depth_start = np.append(depth_start,ovr['depth'][0])
		deoth_end = np.append(depth_end,ovr['depth'][-1])
		hden = np.append(hden,ovr['hden'][0])
		Lya = np.append(Lya, np.log10(simps(np.power(10.,emis['H__1__1216A']),emis['depth'])/(4.*np.pi*1.63e-11)))
		Ha  = np.append(Ha,  np.log10(simps(np.power(10.,emis['H__1__6563A']),emis['depth'])/(4.*np.pi*3.03e-12)))
		CIII = np.append(CIII,np.log10(simps(np.power(10.,emis['C__3_9770A']),emis['depth'])/(4.*np.pi*2.03e-11)))
		civ = simps(np.power(10.,emis['C__4__1548A']),emis['depth'])+simps(np.power(10.,emis['C__4__1551A']),emis['depth'])
		civ = np.log10(civ/(4.*np.pi*1.28e-11))
		CIV = np.append(CIV,civ)
		ovi = simps(np.power(10.,emis['O__6__1032A']),emis['depth'])+simps(np.power(10.,emis['O__6__1038A']),emis['depth'])
		ovi = np.log10(ovi/(4.*np.pi*1.92e-11))
		OVI = np.append(OVI,ovi)
		siiv = simps(np.power(10.,emis['Si_4__1394A']),emis['depth'])+simps(np.power(10.,emis['Si_4__1403A']),emis['depth'])
                siiv = np.log10(siiv/(4.*np.pi*1.42e-11))
                SiIV = np.append(SiIV,siiv)
		SiIII_1207 = np.append(SiIII_1207, np.log10(simps(np.power(10.,emis['Si_3__1207A']),emis['depth'])/(4*np.pi*1.635e-11)))
		SiIII_1883 = np.append(SiIII_1883, np.log10(simps(np.power(10.,emis['Si_3__1883A']),emis['depth'])/(4*np.pi*1.055e-11)))

		Lya_temp = np.append(Lya_temp,Lya[-1])
		Ha_temp = np.append(Ha_temp,Ha[-1])
		CIII_temp = np.append(CIII_temp,CIII[-1])
		CIV_temp = np.append(CIV_temp,CIV[-1])
		OVI_temp = np.append(OVI_temp,OVI[-1])
		SiIV_temp = np.append(SiIV_temp,SiIV[-1])
		SiIII_1207_temp = np.append(SiIII_1207_temp,SiIII_1207[-1])
		SiIII_1883_temp = np.append(SiIII_1883_temp,SiIII_1883[-1])		

		if j=='d':
			Lyamax,Lyamin = np.append(Lyamax,Lya_temp.max()),np.append(Lyamin,Lya_temp.min())
			Hamax,Hamin = np.append(Hamax,Ha_temp.max()),np.append(Hamin,Ha_temp.min())
			CIIImax,CIIImin = np.append(CIIImax,CIII_temp.max()),np.append(CIIImin,CIII_temp.min())
			CIVmax,CIVmin = np.append(CIVmax,CIV_temp.max()),np.append(CIVmin,CIV_temp.min())
			OVImax,OVImin = np.append(OVImax,OVI_temp.max()),np.append(OVImin,OVI_temp.min())				
			SiIVmax,SiIVmin = np.append(SiIVmax,SiIV_temp.max()),np.append(SiIVmin,SiIV_temp.min())
			SiIII_1207max,SiIII_1207min = np.append(SiIII_1207max,SiIII_1207_temp.max()),np.append(SiIII_1207min,SiIII_1207_temp.min())
			SiIII_1883max,SiIII_1883min = np.append(SiIII_1883max,SiIII_1883_temp.max()),np.append(SiIII_1883min,SiIII_1883_temp.min())

		
fini = {}
fini['ID'] = ID
fini['idf'] = idf
fini['T'] = T
fini['hden'] = hden
fini['Lya'] = Lya
fini['Lyamax'],fini['Lyamin'] = Lyamax,Lyamin
fini['Ha'] = Ha
fini['Hamax'],fini['Hamin'] = Hamax,Hamin
fini['CIII'] = CIII
fini['CIIImax'],fini['CIIImin'] = CIIImax,CIIImin
fini['CIV'] = CIV
fini['CIVmax'],fini['CIVmin'] = CIVmax,CIVmin
fini['OVI'] = OVI
fini['OVImax'],fini['OVImin'] = OVImax,OVImin
fini['SiIV'] = SiIV
fini['SiIVmax'],fini['SiIVmin'] = SiIVmax, SiIVmin
fini['SiIII_1207'] = SiIII_1207
fini['SiIII_1207max'],fini['SiIII_1207min'] = SiIII_1207max, SiIII_1207min
fini['SiIII_1883'] = SiIII_1883
fini['SiIII_1883max'],fini['SiIII_1883min'] = SiIII_1883max, SiIII_1883min
fini['depth_start'],fini['depth_end'] = depth_start,depth_end

fileout = 'cloudywerk.cpkl'
cPickle.dump(fini,open(fileout,'wb'),protocol=-1)


