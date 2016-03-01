import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure

def make_Cloudy_table(model_gq,table_index):
	hden_n_bins,hden_min,hden_max = 17, -6, 2
	T_n_bins,T_min,T_max = 51, 3, 8
	patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/grid_galquas/emis/z02/"+model_gq+"/"+model_gq+"_run%i.dat"
	hden=np.linspace(hden_min,hden_max,hden_n_bins)
        T=np.linspace(T_min,T_max, T_n_bins)
        table = np.zeros((hden_n_bins,T_n_bins))
        for i in range(hden_n_bins):
                table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]

        return hden,T,table

def make_ion_table(model_gq,ion,number):
	hden_n_bins, hden_min, hden_max = 17, -6, 2
        T_n_bins, T_min, T_max = 51, 3, 8
        patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/grid_galquas/z02/"+model_gq+"/bert_ion_run%i_"+ion+".dat"
        hden=np.linspace(hden_min,hden_max,hden_n_bins)
        T=np.linspace(T_min,T_max, T_n_bins)
        table = np.zeros((hden_n_bins,T_n_bins))
        for i in range(hden_n_bins):
                table[i,:]=[float(l.split()[number]) for l in open(patt%(i+1)) if l[0] != "#"]
        return hden,T,table

def load_emission_tables(model_gq):
	## CIII Emission
	hden1,T1,table_CIII_977 = make_Cloudy_table(model_gq,7)
	sr_CIII_977 = table_CIII_977.T.ravel()
	CIII_Emission = table_CIII_977
	
	## CIV  Emission
	hden1, T1, table_CIV_1 = make_Cloudy_table(model_gq,3)
	hden1, T1, table_CIV_2 = make_Cloudy_table(model_gq,4)
	sr_CIV_1 = table_CIV_1.T.ravel()
	sr_CIV_2 = table_CIV_2.T.ravel()
	CIV_Emission = 10.**(table_CIV_1) + 10.**(table_CIV_2)
	CIV_Emission = np.log10(CIV_Emission)
	
	## OVI  Emission
	hden1, T1, table_OVI_1 = make_Cloudy_table(model_gq,5)
	hden1, T1, table_OVI_2 = make_Cloudy_table(model_gq,6)
	sr_OVI_1 = table_OVI_1.T.ravel()
	sr_OVI_2 = table_OVI_2.T.ravel()
	OVI_Emission = 10.**(table_OVI_1) + 10.**(table_OVI_2)
	OVI_Emission = np.log10(OVI_Emission)
	
	return CIII_Emission,CIV_Emission,OVI_Emission

def load_ionfrac_tables(model_gq):
        ## CIII Emission
        hden1,T1,table_CIII = make_ion_table(model_gq,'C',3)
        ## CIV  Emission
        hden1, T1, table_CIV = make_ion_table(model_gq,'C',4)
        ## OVI  Emission
        hden1, T1, table_OVI = make_ion_table(model_gq,'O',6)
	## HI Fraction
	hden1, T1, table_HI = make_ion_table(model_gq,'H',1)
	return table_CIII,table_CIV,table_OVI,table_HI


def plot_emission_imshows():
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
	return

def plot_emission_ratio_imshow(): 
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
	
	return

def hden_line(x):
    m = (2.+6.)/(16.)
    b = -6.
    return m*x+b

def T_line(x):
    m = (8.-3.)/(50.)
    b = 3.
    return m*x+b

def plot_contours(contours,color='k'):
	for n,contour in enumerate(contours):
		cx_hden,cx_T = contour[:,0],contour[:,1]
                cy_hden,cy_T = hden_line(cx_hden),T_line(cx_T)
		plt.plot(cy_T,cy_hden,color=color,lw=2.0)
	
	return
		
def plot_emissivity_curves(cval):
	#model_gqs = ['g1q01','g1q1','g1q10']
	model_gqs = ['g1q1']
	i = 0
	while i < len(model_gqs):
		CIII_Emission,CIV_Emission,OVI_Emission = load_emission_tables(model_gqs[i])
		contour_CIII = measure.find_contours(CIII_Emission, cval)	
		contour_CIV  = measure.find_contours(CIV_Emission,cval)
		contour_OVI  = measure.find_contours(OVI_Emission,cval)

		plot_contours(contour_CIII,'b')
		plot_contours(contour_CIV,'g')
		plot_contours(contour_OVI,'r')	
	
		i = i + 1
	
	plt.xlabel('Temperature [k]')
	plt.ylabel('n_H [cm^-3')
	plt.title('Contour Level = '+str(cval))
	plt.savefig('hden_T_contours_'+str(int(cval))+'.png')
	plt.close()

	return


