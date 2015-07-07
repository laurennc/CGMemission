import numpy as np

def logU_to_loghden(x):
	##the Cloudy derived form for z=0.2:
	return -0.1689069978*x - 6.
	
def hden_for_cloudy(werk):
	hlow = logU_to_loghden(werk['Ulow'])
	hhigh = logU_to_loghden(werk['Uhigh'])
	return hlow,hhigh

def write_loop(letter,i,f,werk,hlow,hhigh):
	if letter == "a":
		f.write('hden '+str(hlow[i])+'\n')
                f.write('metals '+str(werk['Zlow'][i])+' log \n')
	elif letter == "b":
		f.write('hden '+str(hlow[i])+'\n')
                f.write('metals '+str(werk['Zhigh'][i])+' log \n')
	elif letter == "c":
		f.write('hden '+str(hhigh[i])+'\n')
                f.write('metals '+str(werk['Zlow'][i])+' log \n')
	elif letter == "d":
		f.write('hden '+str(hhigh[i])+'\n')
                f.write('metals '+str(werk['Zhigh'][i])+' log \n')
	else:
		print "Incorrect Letter Option"
	return

werk = np.genfromtxt('werk14_ionparams.dat',names=True,dtype=None)
hlow,hhigh = hden_for_cloudy(werk)

exeout = open('run_WerkLoop.exe','w')

for i in range(len(werk['SDSSField'])):
	#print str(i)
	for j in "abcd":
		exeout.write('echo '+str(i)+' '+j+' \n')
		exeout.write('/scratch/lauren/codes/cloudy/source/cloudy.exe -p '+str(werk['SDSSField'][i])+'_'+str(werk['ID'][i])+'_'+j+' \n')
		f = open('WerkLoop/'+str(werk['SDSSField'][i])+'_'+str(werk['ID'][i])+'_'+j+'.in','w')
		f.write('table HM05 redshift 0.2 \n')
		f.write('CMB redshift 0.2 \n')
		f.write('stop neutral column density '+str(werk['N_HI'][i]-np.log10(2.))+'\n')
		
		write_loop(j,i,f,werk,hlow,hhigh)		

		f.write('iterate to convergence \n')
		f.write('save last physical conditions ".phy" \n')
		f.write('save last overview ".ovr" \n')
		f.write('save last lines emissivity ".emis" \n')
		f.write('H  1  1216 \n')
		f.write ('H  1  6563 \n')
		f.write('C  3  977 \n')
		f.write('C  4  1548 \n')
		f.write('C  4  1551 \n')
		f.write('O  6  1032 \n')
		f.write('O  6  1038 \n')
		f.write('end of lines \n')
		f.write('save some column densities ".coldens" \n')
		f.write('hydrogen 1 \n')
		f.write('hydrogen 2 \n')
		f.write('carbon 3 \n')
		f.write('carbon 4 \n')
		f.write('oxygen 6 \n')
		f.write('end of column densities \n')
	
f.close()
exeout.close()

