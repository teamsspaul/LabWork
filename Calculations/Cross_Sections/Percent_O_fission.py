import numpy as np
import sys
import matplotlib.pyplot as plt
import math

#This program needs four inputs: t or f (true or false)
#The first corresponds to printing off the group cross sections
#The second corresponds to plotting cross sections and spectrum
#The third corresponds to printing the percent of fissions
#the fourth corresponds to printing the percent fissions in each energy group
#The 5th corresponds to printing flux times cross section summed over all energy groups

#There is an issue in the program, where the I would expect the fission spectrum
#To have the highest percent of U238 fissions, but it does not. Why is that?


#open Fission cross-sections
sigma_f_235 = np.genfromtxt('u235_fission.csv', delimiter=",")
sigma_f_238 = np.genfromtxt('u238_f3.csv', delimiter=",")
sigma_f_239 = np.genfromtxt('pu239_f.csv', delimiter=",")
sigma_f_240 = np.genfromtxt('pu240_fission.csv', delimiter=",")
sigma_f_241 = np.genfromtxt('pu241_fission.csv', delimiter=",")

sigma_a_235 = np.genfromtxt('u235_a.csv', delimiter=",")
sigma_a_238 = np.genfromtxt('u238_a.csv', delimiter=",")
sigma_a_239 = np.genfromtxt('pu239_ab.csv',delimiter=",")
sigma_a_240 = np.genfromtxt('pu240_a.csv',delimiter=",")
sigma_a_241 = np.genfromtxt('pu241_a.csv',delimiter=",")

#The fission and gamma cross sections for U239 and Np239 Am241
#sigma_f_239 = np.genfromtxt('u239_f.csv', delimiter=",")
#sigma_f_239 = np.genfromtxt('np239_f.csv', delimiter=",")
#sigma_a_239 = np.genfromtxt('u239_a.csv',delimiter=",")
#sigma_a_239 = np.genfromtxt('np239_a.csv',delimiter=",")
#sigma_a_241 = np.genfromtxt('am241_g.csv',delimiter=",")
#sigma_f_241 = np.genfromtxt('am241_f.csv', delimiter=",")

sigma_a_106 = np.genfromtxt('Ru_106_a.csv',delimiter=",")

#Create The Flux Spectrum (three of them)
DU_Samples = np.genfromtxt('DU_Samples.csv', delimiter=",")
FBR_Blanket = np.genfromtxt('FBR_Blanket.csv', delimiter=",")

#Determine atoms per gram for each isotope
Na=0.60221409
N235=3.368235853*Na*10**-9
N238=1170.465962*Na*10**-9
N239=18.14187504*Na*10**-9
N240=1.507603392*Na*10**-9
N241=0.624170052*Na*10**-9

#This is PPB like Dr. Chirayath wanted
#N235=791.6833976
#N238=278630.3382
#N239=4336.854352
#N240=361.9059344
#N241=150.4604638

#make interpolation functions
from scipy import interpolate
#Interpolation of cross sections
sig_f_235_interp = interpolate.interp1d(sigma_f_235[:,0], sigma_f_235[:,1],bounds_error=False, fill_value=sigma_f_235[-1,1])
sig_f_238_interp = interpolate.interp1d(sigma_f_238[:,0], sigma_f_238[:,1],bounds_error=False, fill_value=sigma_f_238[-1,1])
sig_f_239_interp = interpolate.interp1d(sigma_f_239[:,0], sigma_f_239[:,1],bounds_error=False, fill_value=sigma_f_239[-1,1])
sig_f_240_interp = interpolate.interp1d(sigma_f_240[:,0], sigma_f_240[:,1],bounds_error=False, fill_value=sigma_f_240[-1,1])
sig_f_241_interp = interpolate.interp1d(sigma_f_241[:,0], sigma_f_241[:,1],bounds_error=False, fill_value=sigma_f_241[-1,1])

sig_a_235_interp = interpolate.interp1d(sigma_a_235[:,0], sigma_a_235[:,1],bounds_error=False, fill_value=sigma_a_235[-1,1])
sig_a_238_interp = interpolate.interp1d(sigma_a_238[:,0], sigma_a_238[:,1],bounds_error=False, fill_value=sigma_a_238[-1,1])
sig_a_239_interp = interpolate.interp1d(sigma_a_239[:,0], sigma_a_239[:,1],bounds_error=False, fill_value=sigma_a_239[-1,1])
sig_a_240_interp = interpolate.interp1d(sigma_a_240[:,0], sigma_a_240[:,1],bounds_error=False, fill_value=sigma_a_240[-1,1])
sig_a_241_interp = interpolate.interp1d(sigma_a_241[:,0], sigma_a_241[:,1],bounds_error=False, fill_value=sigma_a_241[-1,1])

sig_a_106_interp = interpolate.interp1d(sigma_a_106[:,0], sigma_a_106[:,1],bounds_error=False, fill_value=sigma_a_106[-1,1])
#Interpolation of flux spectrum
# 0 = Energy, 1 = Flux MCNP, 2 = Flux at Power?, 3 = Flux/Lethargy MCNP, 4 = Flux/Lethargy at Power? 5 Energy 6 flux (E) - I think
# 7 is flux(E)*E #8 is the same as 2 #9 is the same as 6 (new spectrum) #10 is the same as 7 new spectrum
# 1 and 2 give the same answer, and 3 and 4 give the same answer. There is a difference between the two, but not a huge one
DU_Samples_interp = interpolate.interp1d(DU_Samples[:,5], DU_Samples[:,9],bounds_error=False, fill_value=DU_Samples[-1,1])
FBR_Blanket_interp = interpolate.interp1d(FBR_Blanket[:,5], FBR_Blanket[:,2],bounds_error=False, fill_value=DU_Samples[-1,1])
chi = lambda E:  0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)  #Fission Spectrum (1/MeV) I think this is flux(E)
#chi = lambda E:  0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)*E  #Fission Spectrum (I think this is flux(E)*E (like what all the other spectrums are)

#get the union of the energy grids
energies = np.union1d(sigma_f_235[:,0], sigma_f_238[:,0])
energies_new = np.union1d(energies,sigma_f_239[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_f_240[:,0])
energies = energies_new
#print(len(energies))
energies_new = np.union1d(energies,sigma_f_241[:,0])
energies = energies_new
#print(len(energies))
#Perform Integration
from scipy import integrate
from scipy.integrate import trapz

#Make Functions of things I want to integrate.
#Function with Fission Spectrum
X_f_235_phi=interpolate.interp1d(energies,chi(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi=interpolate.interp1d(energies,chi(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi=interpolate.interp1d(energies,chi(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi=interpolate.interp1d(energies,chi(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi=interpolate.interp1d(energies,chi(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

#Function with DU sample spectrum
X_f_235_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

#These are for gamma+fission = absorption
X_a_235_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_235_interp(energies)+sig_f_235_interp(energies)),fill_value=0,bounds_error=False)
X_a_238_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_238_interp(energies)+sig_f_238_interp(energies)),fill_value=0,bounds_error=False)
X_a_239_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_239_interp(energies)+sig_f_239_interp(energies)),fill_value=0,bounds_error=False)
X_a_240_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_240_interp(energies)+sig_f_240_interp(energies)),fill_value=0,bounds_error=False)
X_a_241_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_241_interp(energies)+sig_f_241_interp(energies)),fill_value=0,bounds_error=False)

#These are for gamma 
#X_a_235_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_235_interp(energies)),fill_value=0,bounds_error=False)
#X_a_238_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_238_interp(energies)),fill_value=0,bounds_error=False)
#X_a_239_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_239_interp(energies)),fill_value=0,bounds_error=False)
#X_a_240_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_240_interp(energies)),fill_value=0,bounds_error=False)
#X_a_241_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_a_241_interp(energies)),fill_value=0,bounds_error=False)



X_a_106_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_a_106_interp(energies),fill_value=0,bounds_error=False)
#Function with FBR spectrum
X_f_235_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

#Perform the integration for each group
#Groups=[4e-7,0.1,0.25,1,20]
#Groups=np.logspace(-10,1.30103,100)
#Groups=np.concatenate((np.logspace(-11,-5,60),np.logspace(-5.001,-1.69897,330),np.logspace(-1.69898,1.30103,30)),axis=1)
Groups=[20]
#Groups=[0.001,20] #This is for the percent of "thermal" and fast fissions

index=np.zeros(len(Groups)+1)

for jjj in range(1,len(Groups)+1):
	for xx in range(0,len(energies)):
		if(energies[xx]<=Groups[jjj-1]):
			index[jjj]=xx

#for i in range(0,len(Groups)):

#For the overall fissions
f239=np.zeros(3)
f238=np.zeros(3)
f235=np.zeros(3)
f240=np.zeros(3)
f241=np.zeros(3)
Tot=np.zeros(3)
#These are for determining individual fission percentages (thermal vs fast)
fi235=np.zeros((len(Groups),3))
fi238=np.zeros((len(Groups),3))
fi239=np.zeros((len(Groups),3))
fi240=np.zeros((len(Groups),3))
fi241=np.zeros((len(Groups),3))

#These are for flux times cross section
ff235=np.zeros(3)
ff238=np.zeros(3)
ff239=np.zeros(3)
ff240=np.zeros(3)
ff241=np.zeros(3)



fa235=np.zeros(3)
fa238=np.zeros(3)
fa239=np.zeros(3)
fa240=np.zeros(3)
fa241=np.zeros(3)





X_f_avg_235=np.zeros((3,len(Groups)))
X_f_avg_238=np.zeros((3,len(Groups)))
X_f_avg_239=np.zeros((3,len(Groups)))
X_f_avg_240=np.zeros((3,len(Groups)))
X_f_avg_241=np.zeros((3,len(Groups)))
phi_avg=np.zeros((3,len(Groups)))
for i in range(0,len(Groups)):
	#Perform the integral
	phi_int_g1=integrate.trapz(chi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_235_g1=integrate.trapz(X_f_235_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_238_g1=integrate.trapz(X_f_238_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_239_g1=integrate.trapz(X_f_239_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_240_g1=integrate.trapz(X_f_240_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_241_g1=integrate.trapz(X_f_241_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	
	phi_int_g12=integrate.trapz(DU_Samples_interp(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_235_g12=integrate.trapz(X_f_235_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_238_g12=integrate.trapz(X_f_238_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_239_g12=integrate.trapz(X_f_239_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_240_g12=integrate.trapz(X_f_240_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_241_g12=integrate.trapz(X_f_241_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	
	X_a_235_g12=integrate.trapz(X_a_235_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_a_238_g12=integrate.trapz(X_a_238_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_a_239_g12=integrate.trapz(X_a_239_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_a_240_g12=integrate.trapz(X_a_240_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_a_241_g12=integrate.trapz(X_a_241_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	
	X_a_106_g12=integrate.trapz(X_a_106_phi2(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	
	phi_int_g13=integrate.trapz(FBR_Blanket_interp(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_235_g13=integrate.trapz(X_f_235_phi3(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_238_g13=integrate.trapz(X_f_238_phi3(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_239_g13=integrate.trapz(X_f_239_phi3(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_240_g13=integrate.trapz(X_f_240_phi3(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_241_g13=integrate.trapz(X_f_241_phi3(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	
	#Store all flux values in one place
	if(energies[index[i+1]]-energies[index[i]]>0):
		phi_avg[0,i]=phi_int_g1/(energies[index[i+1]]-energies[index[i]])
		phi_avg[1,i]=phi_int_g12/(energies[index[i+1]]-energies[index[i]])
		phi_avg[2,i]=phi_int_g13/(energies[index[i+1]]-energies[index[i]])
	#print ("Group ", i+1, "Upper Energy ", energies[index[i+1]], "(MeV) Lower Energy ", energies[index[i]], "(MeV)")
	
	
	
	
	
	
	#Calculate fission fraction and the group cross sections
	for j in range(0,3):
		if(j==1):
			#Fission Fraction DUO2 Spectrum
			f239[j]=X_f_239_g12*N239+f239[j]  ### To calculate the fission fraction
			f238[j]=X_f_238_g12*N238+f238[j]
			f235[j]=X_f_235_g12*N235+f235[j]
			f240[j]=X_f_240_g12*N240+f240[j]
			f241[j]=X_f_241_g12*N241+f241[j]
			Tot[j]=Tot[j]+X_f_239_g12*N239+X_f_235_g12*N235+X_f_238_g12*N238+X_f_240_g12*N240+X_f_241_g12*N241
			
			fi239[i,j]=X_f_239_g12*N239
			fi238[i,j]=X_f_238_g12*N238
			fi235[i,j]=X_f_235_g12*N235
			fi240[i,j]=X_f_240_g12*N240
			fi241[i,j]=X_f_241_g12*N241
			
			ff239[j]=X_f_239_g12*10**-24+ff239[j]  ### To calculate flux times cross section
			ff238[j]=X_f_238_g12*10**-24+ff238[j]
			ff235[j]=X_f_235_g12*10**-24+ff235[j]
			ff240[j]=X_f_240_g12*10**-24+ff240[j]
			ff241[j]=X_f_241_g12*10**-24+ff241[j]
			
			fa239[j]=X_a_239_g12*10**-24+fa239[j]  ### To calculate flux times cross section
			fa238[j]=X_a_238_g12*10**-24+fa238[j]
			fa235[j]=X_a_235_g12*10**-24+fa235[j]
			fa240[j]=X_a_240_g12*10**-24+fa240[j]
			fa241[j]=X_a_241_g12*10**-24+fa241[j]
			
			
			
			
			
			
			
			
			
			#Cross Section
			if (phi_int_g12!=0):
				X_f_avg_235[j,i]=X_f_235_g12/phi_int_g12
				X_f_avg_238[j,i]=X_f_238_g12/phi_int_g12
				X_f_avg_239[j,i]=X_f_239_g12/phi_int_g12
				X_f_avg_240[j,i]=X_f_240_g12/phi_int_g12
				X_f_avg_241[j,i]=X_f_241_g12/phi_int_g12
			
		if(j==0):
			#Fission fraction (Fission Spectrum)
			f239[j]=X_f_239_g1*N239+f239[j]
			f238[j]=X_f_238_g1*N238+f238[j]
			f235[j]=X_f_235_g1*N235+f235[j]
			f240[j]=X_f_240_g1*N240+f240[j]
			f241[j]=X_f_241_g1*N241+f241[j]
			Tot[j]=Tot[j]+X_f_239_g1*N239+X_f_235_g1*N235+X_f_238_g1*N238+X_f_240_g1*N240+X_f_241_g1*N241
			
			fi239[i,j]=X_f_239_g1*N239
			fi238[i,j]=X_f_238_g1*N238
			fi235[i,j]=X_f_235_g1*N235
			fi240[i,j]=X_f_240_g1*N240
			fi241[i,j]=X_f_241_g1*N241
			
			
			#Cross Section
			if (phi_int_g1!=0):
				X_f_avg_235[j,i]=X_f_235_g1/phi_int_g1
				X_f_avg_238[j,i]=X_f_238_g1/phi_int_g1
				X_f_avg_239[j,i]=X_f_239_g1/phi_int_g1
				X_f_avg_240[j,i]=X_f_240_g1/phi_int_g1
				X_f_avg_241[j,i]=X_f_241_g1/phi_int_g1
		if(j==2):
			#Fission fraction FBR Spectrum
			f239[j]=X_f_239_g13*N239+f239[j]
			f238[j]=X_f_238_g13*N238+f238[j]
			f235[j]=X_f_235_g13*N235+f235[j]
			f240[j]=X_f_240_g13*N240+f240[j]
			f241[j]=X_f_241_g13*N241+f241[j]
			Tot[j]=Tot[j]+X_f_239_g13*N239+X_f_235_g13*N235+X_f_238_g13*N238+X_f_240_g13*N240+X_f_241_g13*N241
			
			
			fi239[i,j]=X_f_239_g13*N239
			fi238[i,j]=X_f_238_g13*N238
			fi235[i,j]=X_f_235_g13*N235
			fi240[i,j]=X_f_240_g13*N240
			fi241[i,j]=X_f_241_g13*N241
			
			#Cross Section FBR Specrum
			if (phi_int_g13!=0):
				X_f_avg_235[j,i]=X_f_235_g13/phi_int_g13
				X_f_avg_238[j,i]=X_f_238_g13/phi_int_g13
				X_f_avg_239[j,i]=X_f_239_g13/phi_int_g13
				X_f_avg_240[j,i]=X_f_240_g13/phi_int_g13
				X_f_avg_241[j,i]=X_f_241_g13/phi_int_g13
	

	
	if (sys.argv[1]=="t"):
		print (" ")
		print ("Group ", i+1, "Upper Energy ", energies[index[i+1]])
		if (phi_int_g1!=0):
			print ("Fission Cross Section")
			print ("Flux is ", phi_int_g1)
			print ("235 ", X_f_235_g1/phi_int_g1)
			print ("238 ", X_f_238_g1/phi_int_g1)
			print ("239 ", X_f_239_g1/phi_int_g1)
			print ("240 ", X_f_240_g1/phi_int_g1)
			print ("241 ", X_f_241_g1/phi_int_g1)
		
		if (phi_int_g12!=0):
			print (" ")
			print ("DUO2 Samples")
			print ("Flux is ", phi_int_g12)
			print ("235  ", X_f_235_g12/phi_int_g12)
			print ("238  ", X_f_238_g12/phi_int_g12)
			print ("239  ", X_f_239_g12/phi_int_g12)
			print ("240  ", X_f_240_g12/phi_int_g12)
			print ("241  ", X_f_241_g12/phi_int_g12)
			print (" ")
			print ("235a ", X_a_235_g12/phi_int_g12)
			print ("238a ", X_a_238_g12/phi_int_g12)
			print ("239a ", X_a_239_g12/phi_int_g12)
			print ("240a ", X_a_240_g12/phi_int_g12)
			print ("241a ", X_a_241_g12/phi_int_g12)
			print ("106  ", X_a_106_g12/phi_int_g12)
		
		if (phi_int_g13!=0):
			print (" ")
			print ("FBR Blanket")
			print ("Flux is ", phi_int_g13)
			print ("235 ", X_f_235_g13/phi_int_g13)
			print ("238 ", X_f_238_g13/phi_int_g13)
			print ("239 ", X_f_239_g13/phi_int_g13)
			print ("240 ", X_f_240_g13/phi_int_g13)
			print ("241 ", X_f_241_g13/phi_int_g13)
		print(" ")
	


#Plot Everything
if (sys.argv[2]=='t'):
	fig,ax1=plt.subplots()
	#Plot different spectrum
	#ax1.semilogx(energies,chi(energies)/np.sum(chi(energies)),'b--',linewidth=1,label="Fission Spectrum",rasterized=True)
	ax1.semilogx(energies,(DU_Samples_interp(energies))/np.sum(DU_Samples_interp(energies)),'b',linewidth=2,label="DU_Samples",rasterized=True)
	#ax1.semilogx(energies,FBR_Blanket_interp(energies)/np.sum(FBR_Blanket_interp(energies)),'b',linewidth=1,label="FBR_Blanket",rasterized=True)
	for i in range(0,len(Groups)): #Plots the energy groups and the average flux values for each spectrum
		ax1.semilogx(energies[index[i+1]]*(energies/energies),np.linspace(0,0.0016,len(energies)),'m',linewidth=2,rasterized=True)  #0 is fission, 1 is DU, 2 is FBR
		#ax1.semilogx(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),(phi_avg[0,i]*np.ones(len(energies[index[i]:index[i+1]])))/np.sum(chi(energies)),'b--',linewidth=1,rasterized=True)
		ax1.semilogx(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),phi_avg[1,i]*np.ones(len(energies[index[i]:index[i+1]]))/np.sum(DU_Samples_interp(energies)),'b',linewidth=2,rasterized=True)
		#ax1.semilogx(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),phi_avg[2,i]*np.ones(len(energies[index[i]:index[i+1]]))/np.sum(FBR_Blanket_interp(energies)),'b',linewidth=1,rasterized=True)
	#label some stuff
	ax1.set_xlabel("E (MeV)")
	ax1.set_ylabel("$\phi$(E)$\cdot$E Normalized (n/cm$^{2}$s)",color='b')
	ax1.legend(loc=1)
	for t1 in ax1.get_yticklabels():
		t1.set_color('b')
	ax2=ax1.twinx()
	#Plot the Cross Sections
	ax2.loglog(energies, sig_f_235_interp(energies),'r',label="$\sigma_\mathrm{fU235}$",rasterized=True)
	#ax2.loglog(energies, sig_f_238_interp(energies),'c',label="$\sigma_\mathrm{fU238}$",rasterized=True)
	ax2.loglog(energies, sig_a_238_interp(energies),'k',label="$\sigma_\mathrm{aU238}$",rasterized=True) #U238 absorption
	#for i in range(0,len(Groups)): #plots the average x section values for each spectrum in each energy group (U238) #0 is fission, 1 is DU, 2 is FBR
		#ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_238[0,i]*np.ones(len(energies[index[i]:index[i+1]])),'c',rasterized=True)
		#ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_238[1,i]*np.ones(len(energies[index[i]:index[i+1]])),'c',rasterized=True)
		#ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_238[2,i]*np.ones(len(energies[index[i]:index[i+1]])),'c',rasterized=True)
	ax2.loglog(energies, sig_f_239_interp(energies),'g',label="$\sigma_\mathrm{fPu239}$",rasterized=True)
	for i in range(0,len(Groups)): #Plots the average x section values for each spectrum in each energy group (Pu239) #0 is fission, 1 is DU, 2 is FBR
		#ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_239[0,i]*np.ones(len(energies[index[i]:index[i+1]])),'g',rasterized=True)
		ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_239[1,i]*np.ones(len(energies[index[i]:index[i+1]])),'g',rasterized=True)
		#ax2.loglog(np.logspace(math.log10(energies[index[i]]),math.log10(energies[index[i+1]]),len(energies[index[i]:index[i+1]])),X_f_avg_239[2,i]*np.ones(len(energies[index[i]:index[i+1]])),'g',rasterized=True)
	ax2.loglog(energies, sig_f_240_interp(energies),'k',label="$\sigma_\mathrm{fPu240}$",rasterized=True)
	ax2.loglog(energies, sig_f_241_interp(energies),'y',label="$\sigma_\mathrm{fPu241}$",rasterized=True)
	ax2.set_ylabel("$\sigma$ (barns)")
	
	ax2.legend(loc=3)
	plt.savefig("Spectrum.pdf")
	plt.clf()

	
#show the percentages
print (" ")
if (sys.argv[3]=='t'):
	for j in range(0,3): 
	#for j in range(1,2): #Only the DUO2 Sample
		if (j==0): print ("Fission Spectrum")
		if (j==1): print ("DU Spectrum")
		if (j==2): print ("FBR Spectrum")
		print ("f235 ",f235[j]/Tot[j])
		print ("f238 ",f238[j]/Tot[j])
		print ("f239 ",f239[j]/Tot[j])
		print ("f240 ",f240[j]/Tot[j])
		print ("f241 ",f241[j]/Tot[j])
		print ("Total", (f235[j]+f238[j]+f239[j]+f240[j]+f241[j])/Tot[j])
		print ("Real Total", Tot[j])
		print (" ")
		
#Show the other percentages (only for DUO2 Spectrum but could be for others)
print (" ")
if (sys.argv[4]=='t'):
	#for j in range(0,3): #Loop over spectra
	for j in range(1,2): #Only the DUO2 Sample
		if (j==0): print ("Fission Spectrum")
		if (j==1): print ("DU Spectrum")
		if (j==2): print ("FBR Spectrum")
		for i in range(0,len(Groups)):#Loop over groups
			print ("Group ", i+1, "Upper Energy ", energies[index[i+1]], "(MeV) Lower Energy ", energies[index[i]], "(MeV)")
			print ("f235 ",fi235[i,j]/f235[j])
			print ("f238 ",fi238[i,j]/f238[j])
			print ("f239 ",fi239[i,j]/f239[j])
			print ("f240 ",fi240[i,j]/f240[j])
			print ("f241 ",fi241[i,j]/f241[j])
		print (" ")
		
#Display the flux times cross section summed over all energy groups
if (sys.argv[5]=='t'):
	#for j in range(0,3): #Loop over spectra
	for j in range(1,2): #Only the DUO2 Sample
		if (j==0): print ("Fission Spectrum")
		if (j==1): print ("DU Spectrum")
		if (j==2): print ("FBR Spectrum")
		print ("Group ", i+1, "Upper Energy ", energies[index[i+1]], "(MeV) Lower Energy ", energies[index[i]], "(MeV)")
		print ("ff235 ",ff235[j])
		print ("ff238 ",ff238[j])
		print ("ff239 ",ff239[j])
		print ("ff240 ",ff240[j])
		print ("ff241 ",ff241[j])
		print (" ")
		print ("fa235 ",fa235[j])
		print ("fa238 ",fa238[j])
		print ("fa239 ",fa239[j])
		print ("fa240 ",fa240[j])
		print ("fa241 ",fa241[j])
		print (" ")