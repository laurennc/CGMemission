import numpy as np
import sys
#sys.path.insert(0,'/u/10/l/lnc2115/astrostats/CGMemission')
#from lauren import *

#Some constants in cgs units
HI_mass = 1.6737236e-24 #g
electron_mass = 9.10938291e-28 #g
proton_mass = 1.67262178e-24 #g 
H2_mass = 3.3473e-24 #g
simunits_to_cm = 1.10202798453e+26 #cm
Z_solar = 0.0127

def _SolarMetals(field,data):
	return np.log10(data['Metallicity']) - np.log10(0.02)

add_field('SolarMetals',take_log=False,function=_SolarMetals)#,units='Solar')

def _HAlphaEmissionSr(field,data):
	Tfour = data['Temperature']/(1.0e4)
	ne = data['Electron_NumberDensity']
	nHII = data['HII_Density']/proton_mass
	exponent = -0.942-0.031*np.log(Tfour)
	return (2.82*10**(-26.0))*ne*nHII*Tfour**(exponent)

add_field('HAlphaEmissionSr',take_log=True,function=_HAlphaEmissionSr,units='ergs cm^{-3} s^{-1} sr^{-1}')

def _HAlphaEmissionArc(field,data):
	return data['HAlphaEmissionSr']*(2.3528*10**-11.0)

add_field('HAlphaEmissionArc',take_log=True,function=_HAlphaEmissionArc,units=r'ergs cm^{-3} s^{-1} arcsec^{-2}')

def _HAlphaEmissionRal(field,data):
	return (4151928.85099)*data['HAlphaEmissionSr']
#from David
#	return (5.7e-18)*data['HAlphaEmisisonArc']

add_field('HAlphaEmissionRal',take_log=True,function=_HAlphaEmissionRal,units='photons cm^{-3} s^{-1} sr^{-1}')

def _EmissionMeasureCM(field,data):
	return (data['Electron_Density']/electron_mass)*(data['HII_Density']/proton_mass)

add_field('EmissionMeasureCM',take_log=True,function=_EmissionMeasureCM,units='cm^{-6}')

def _EmissionMeasurePC(field,data):
	return (3.2408*10**-19.0)*data['EmissionMeasureCM']

add_field('EmissionMeasurePC',take_log=True,function=_EmissionMeasurePC,units='cm^{-6} pc')

def _EmissionMeasurePCDX(field,data):
	return data['EmissionMeasurePC']*(data['CellVolume']**(1.0/3.0))

add_field('EmissionMeasurePCDX',take_log=True,function=_EmissionMeasurePCDX,units='cm^{-6} pc')

def _EmissionMeasureCold(field,data):
	idx = np.where(np.log10(data['Temperature']) > 4.5)[0]
	em = data['EmissionMeasurePC']
	em[idx] = 0.0
	return em

add_field('EmissionMeasureCold',take_log=True,function=_EmissionMeasureCold,units='cm^{-5} pc')

def _HAlpha_Emissivity(field,data):
   hden_n_bins,hden_min, hden_max = 15,-6,1
   T_n_bins, T_min, T_max = 51,3,8
   patt =  "/hpc/astrostats/astro/users/lnc2115/codes/cloudy_yt/yuan/yuan_hemissivity_run%i.dat"  
   from scipy import interpolate
   hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
   T=numpy.linspace(T_min,T_max, T_n_bins)
   table = numpy.zeros((hden_n_bins, T_n_bins))
   for i in range(hden_n_bins):
       table[i,:]=[float(l.split()[-1]) for l in open(patt%(i+1)) if l[0] != "#"]
   sp=interpolate.RectBivariateSpline(hden,T,table)
   good=data["Temperature"].shape
   H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
   Temperature=numpy.log10(numpy.array(data["Temperature"]))
   H_N=H_N.reshape(H_N.size)
   Temperature=Temperature.reshape(Temperature.size)
   dia=sp.ev(H_N,Temperature)
   dia=dia.reshape(good)
   Halpha=(10.0**dia)*(data["H_NumberDensity"]**2.0)
   return Halpha*1.87e-12   # not yet

add_field("HAlpha_Emissivity",units=r"\rm{ergs s^{-1}cm^{-3}arcsec^{-2}",function=_HAlpha_Emissivity)

def _HAlpha_Emissivity_R(field,data):
	return data['HAlpha_Emissivity']/(5.7e-18)

add_field("HAlpha_Emissivity_R",units="R",function=_HAlpha_Emissivity_R)

def _Emission_HAlpha(field,data):
   from lauren import make_Cloudy_table
   from scipy import interpolate
   from lauren import scale_by_metallicity
   hden,T,table = make_Cloudy_table(2)
   sp=interpolate.RectBivariateSpline(hden,T,table)
   good=data["Temperature"].shape
   H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
   Temperature=numpy.log10(numpy.array(data["Temperature"]))
   H_N=H_N.reshape(H_N.size)
   Temperature=Temperature.reshape(Temperature.size)
   dia=sp.ev(H_N,Temperature)
   dia=dia.reshape(good)
   Halpha=(10.0**dia)*(data["H_NumberDensity"]**2.0)
   Halpha = Halpha/(4.*np.pi*3.03e-12)	
   return scale_by_metallicity(Halpha,0.0,np.log10(data['Metallicity']))

add_field("Emission_HAlpha",units=r"\rm{ergs s^{-1}cm^{-3}sr^{-1}",function=_Emission_HAlpha)

def _Emission_CIV(field,data):
	from lauren import make_Cloudy_table
	from scipy import interpolate
	from lauren import scale_by_metallicity
	hden1,T1,table1 = make_Cloudy_table(3)
	hden2,T2,table2 = make_Cloudy_table(4)
	sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
	sp2 = interpolate.RectBivariateSpline(hden2,T2,table2)
	good = data["Temperature"].shape
	H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
	Temperature=numpy.log10(numpy.array(data["Temperature"]))
	H_N=H_N.reshape(H_N.size)
	Temperature=Temperature.reshape(Temperature.size)
	dia1=sp1.ev(H_N,Temperature)
	dia1=dia1.reshape(good)
	dia2=sp2.ev(H_N,Temperature)
        dia2=dia2.reshape(good)
	CIV=((10.0**dia1)+(10**dia2))*(data["H_NumberDensity"]**2.0)
	CIV = CIV/(4.*np.pi*1.28e-11)
	return scale_by_metallicity(CIV,0.0,np.log10(data['Metallicity']))

add_field("Emission_CIV",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_CIV)

def _Emission_CIV_ncut(field,data):
	idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
	datawant = data['Emission_CIV']
	datawant[idx] = 0.0
	return datawant

add_field("Emission_CIV_ncut",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_CIV_ncut)

def _Emission_OVI(field,data):
        from lauren import make_Cloudy_table
	from scipy import interpolate
	from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(5)
        hden2,T2,table2 = make_Cloudy_table(6)
        #sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        #sp2 = interpolate.RectBivariateSpline(hden2,T2,table2)
        h1,T1 = np.meshgrid(hden1,T1)
	h2,T2 = np.meshgrid(hden2,T2)
	pts1 = np.array((h1.ravel(),T1.ravel())).T
	pts2 = np.array((h2.ravel(),T2.ravel())).T
        sr1 = table1.T.ravel()
	sr2 = table2.T.ravel()
	bl1 = interpolate.LinearNDInterpolator(pts1,sr1)
	bl2 = interpolate.LinearNDInterpolator(pts2,sr2)
	H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        dia1 = bl1(H_N,Temperature)
	dia2 = bl2(H_N,Temperature)
        idx = np.isnan(dia1)
	dia1[idx] = -300.0
	dia2[idx] = -300.0
	OVI=((10.0**dia1)+(10**dia2))*(data["H_NumberDensity"]**2.0)
	OVI =  OVI/(4.*np.pi*1.92e-11)
	return scale_by_metallicity(OVI,-0.3,np.log10(data['Metallicity']))

add_field("Emission_OVI",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_OVI)


def _Emission_OVI_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_OVI']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_OVI_ncut",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_OVI_ncut)


def _HAlpha_Voort_R(field,data):
	idx = np.where(np.log10(data['H_NumberDensity']) > -1)[0]
	halpha = data['HAlpha_Emissivity_R']
	halpha[idx] = 0
	return halpha

add_field("HAlpha_Voort_R",units="R",function=_HAlpha_Voort_R)

def _Emission_CIII_977(field,data):
   from lauren import make_Cloudy_table
   from scipy import interpolate
   from lauren import scale_by_metallicity
   hden1,T1,table1 = make_Cloudy_table(7)
   sp=interpolate.RectBivariateSpline(hden1,T1,table1)
   good=data["Temperature"].shape
   H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
   Temperature=numpy.log10(numpy.array(data["Temperature"]))
   H_N=H_N.reshape(H_N.size)
   Temperature=Temperature.reshape(Temperature.size)
   dia=sp.ev(H_N,Temperature)
   dia=dia.reshape(good)
   CIII_977=(10.0**dia)*(data["H_NumberDensity"]**2.0)
   CIII_977 = CIII_977/(4.*np.pi*2.03e-11)
   return scale_by_metallicity(CIII_977,0.0,np.log10(data['Metallicity']))

add_field("Emission_CIII_977",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_CIII_977)

def _Emission_CIII_977_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_CIII_977']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_CIII_977_ncut",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_CIII_977_ncut)

def _Emission_CIII(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
	from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(8)
        hden2,T2,table2 = make_Cloudy_table(9)        
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        sp2 = interpolate.RectBivariateSpline(hden2,T2,table2)
        good = data["Temperature"].shape        
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        dia2=sp2.ev(H_N,Temperature)
        dia2=dia2.reshape(good)
        CIII=((10.0**dia1)+(10**dia2))*(data["H_NumberDensity"]**2.0)        
        CIII = CIII/(4.*np.pi*1.04e-11)
	return scale_by_metallicity(CIII,0.0,np.log10(data['Metallicity']))


add_field("Emission_CIII",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_CIII)

def _Emission_CIII_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_CIII']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_CIII_ncut",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_CIII_ncut)

def _Emission_SiII(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
        from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(10)
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        SiII=(10.0**dia1)*(data["H_NumberDensity"]**2.0)
        SiII = SiII/(4.*np.pi*1.095e-11)
        return scale_by_metallicity(SiII,-0.3,np.log10(data['Metallicity']))

add_field("Emission_SiII",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiII)

def _Emission_SiII_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_SiII']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_SiII_ncut",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiII_ncut)

def _Emission_SiIII_1207(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
        from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(11)
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        SiIII_1207=(10.0**dia1)*(data["H_NumberDensity"]**2.0)
        SiIII_1207 = SiIII_1207/(4.*np.pi*1.65e-11)
        return scale_by_metallicity(SiIII_1207,-0.3,np.log10(data['Metallicity']))

add_field("Emission_SiIII_1207",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiIII_1207)

def _Emission_SiIII_1207_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_SiIII_1207']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_SiIII_1207_ncut",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiIII_1207_ncut)

def _Emission_SiIII_1883(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
        from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(12)
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        SiIII_1883=(10.0**dia1)*(data["H_NumberDensity"]**2.0)
        SiIII_1883 = SiIII_1883/(4.*np.pi*1.05e-11)
        return scale_by_metallicity(SiIII_1883,-0.3,np.log10(data['Metallicity']))

add_field("Emission_SiIII_1883",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiIII_1883)

def _Emission_SiIII_1883_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_SiIII_1883']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_SiIII_1883_ncut",units=r"\rm{photons s^{-1}cm^{-3}sr^{-1}",function=_Emission_SiIII_1883_ncut)

def _Emission_SiIV(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
        from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(13)
        hden2,T2,table2 = make_Cloudy_table(14)
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        sp2 = interpolate.RectBivariateSpline(hden2,T2,table2)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        dia2=sp2.ev(H_N,Temperature)
        dia2=dia2.reshape(good)
        SiIV=((10.0**dia1)+(10**dia2))*(data["H_NumberDensity"]**2.0)
        SiIV = SiIV/(4.*np.pi*1.42e-11)
        return scale_by_metallicity(SiIV,-0.3,np.log10(data['Metallicity']))

add_field("Emission_SiIV",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_SiIV)

def _Emission_SiIV_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_SiIV']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_SiIV_ncut",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_SiIV_ncut)

def _Emission_MgII(field,data):
        from lauren import make_Cloudy_table
        from scipy import interpolate
        from lauren import scale_by_metallicity
        hden1,T1,table1 = make_Cloudy_table(15)
        hden2,T2,table2 = make_Cloudy_table(16)
        sp1 = interpolate.RectBivariateSpline(hden1,T1,table1)
        sp2 = interpolate.RectBivariateSpline(hden2,T2,table2)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia1=sp1.ev(H_N,Temperature)
        dia1=dia1.reshape(good)
        dia2=sp2.ev(H_N,Temperature)
        dia2=dia2.reshape(good)
        MgII=((10.0**dia1)+(10**dia2))*(data["H_NumberDensity"]**2.0)
        MgII = MgII/(4.*np.pi*7.10e-12)
        return scale_by_metallicity(MgII,-0.3,np.log10(data['Metallicity']))

add_field("Emission_MgII",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_MgII)

def _Emission_MgII_ncut(field,data):
        idx = np.where(np.log10(data['H_NumberDensity']) < -4.)
        datawant = data['Emission_MgII']
        datawant[idx] = 0.0
        return datawant

add_field("Emission_MgII_ncut",units=r"\rm{ph \ s }^{-1} \rm{cm}^{-2} \rm{sr}^{-1}",function=_Emission_MgII_ncut)























def _CIII_Density(field,data):
	from lauren import make_ion_table
	from scipy import interpolate
	hden, T, table = make_ion_table('C',3)
	sp = interpolate.RectBivariateSpline(hden,T,table)
	good = data["Temperature"].shape
	H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
	return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(2.45e-4/Z_solar)

add_field("CIII_Density",units=r"\rm{cm^{-3}}",function=_CIII_Density)

def _CIV_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('C',4)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(2.45e-4/Z_solar)

add_field("CIV_Density",units=r"\rm{cm^{-3}}",function=_CIV_Density)

def _OVI_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('O',6)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(4.90e-4/Z_solar)

add_field("OVI_Density",units=r"\rm{cm^{-3}}",function=_OVI_Density)

def _MgII_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('Mg',2)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("MgII_Density",units=r"\rm{cm^{-3}}",function=_MgII_Density)


def _SiII_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('Si',2)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiII_Density",units=r"\rm{cm^{-3}}",function=_SiII_Density)

def _SiIII_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('Si',3)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiIII_Density",units=r"\rm{cm^{-3}}",function=_SiIII_Density)

def _SiIV_Density(field,data):
        from lauren import make_ion_table
        from scipy import interpolate
        hden, T, table = make_ion_table('Si',4)
        sp = interpolate.RectBivariateSpline(hden,T,table)
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=sp.ev(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiIV_Density",units=r"\rm{cm^{-3}}",function=_SiIV_Density)


