hden1, T1, table_CIIIdens = make_ion_table('C',3)
sr_CIIIdens = table_CIIIdens.T.ravel()
bl_CIIIdens = interpolate.LinearNDInterpolator(pts,sr_CIIIdens)

def _CIII_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_CIIIdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(2.45e-4/Z_solar)

add_field("CIII_Density",units=r"\rm{cm^{-3}}",function=_CIII_Density)

hden1, T1, table_CIVdens = make_ion_table('C',4)
sr_CIVdens = table_CIVdens.T.ravel()
bl_CIVdens = interpolate.LinearNDInterpolator(pts,sr_CIVdens)

def _CIV_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_CIVdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(2.45e-4/Z_solar)

add_field("CIV_Density",units=r"\rm{cm^{-3}}",function=_CIV_Density)

hden1, T1, table_OVIdens = make_ion_table('O',6)
sr_OVIdens = table_OVIdens.T.ravel()
bl_OVIdens = interpolate.LinearNDInterpolator(pts,sr_OVIdens)

def _OVI_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_OVIdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(4.90e-4/Z_solar)

add_field("OVI_Density",units=r"\rm{cm^{-3}}",function=_OVI_Density)

hden1, T1, table_MgIIdens = make_ion_table('Mg',2)
sr_MgIIdens = table_MgIIdens.T.ravel()
bl_MgIIdens = interpolate.LinearNDInterpolator(pts,sr_MgIIdens)

def _MgII_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_MgIIdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("MgII_Density",units=r"\rm{cm^{-3}}",function=_MgII_Density)

hden1, T1, table_SiIIdens = make_ion_table('Si',2)
sr_SiIIdens = table_SiIIdens.T.ravel()
bl_SiIIdens = interpolate.LinearNDInterpolator(pts,sr_SiIIdens)

def  _SiII_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_SiIIdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiII_Density",units=r"\rm{cm^{-3}}",function=_SiII_Density)

hden1, T1, table_SiIIIdens = make_ion_table('Si',3)
sr_SiIIIdens = table_SiIIIdens.T.ravel()
bl_SiIIIdens = interpolate.LinearNDInterpolator(pts,sr_SiIIIdens)

def _SiIII_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_SiIIIdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiIII_Density",units=r"\rm{cm^{-3}}",function=_SiIII_Density)

hden1, T1, table_SiIVdens = make_ion_table('Si',4)
sr_SiIVdens = table_SiIVdens.T.ravel()
bl_SiIVdens = interpolate.LinearNDInterpolator(pts,sr_SiIVdens)

def _SiIV_Density(field,data):
        good = data["Temperature"].shape
        H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
        Temperature=numpy.log10(numpy.array(data["Temperature"]))
        H_N=H_N.reshape(H_N.size)
        Temperature=Temperature.reshape(Temperature.size)
        dia=bl_SiIVdens(H_N,Temperature)
        dia=dia.reshape(good)
        return (10.0**dia)*data['Metallicity']*data['H_NumberDensity']*(3.47e-5/Z_solar)

add_field("SiIV_Density",units=r"\rm{cm^{-3}}",function=_SiIV_Density)


