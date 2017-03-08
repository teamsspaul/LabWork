import numpy as np
import sys
import matplotlib.pyplot as plt

#open Fission cross-sections
sigma_f_235 = np.genfromtxt('u235_fission.csv', delimiter=",")
sigma_f_238 = np.genfromtxt('u238_fission.csv', delimiter=",")
sigma_f_239 = np.genfromtxt('pu239_fission.csv', delimiter=",")
sigma_f_240 = np.genfromtxt('pu240_fission.csv', delimiter=",")
sigma_f_241 = np.genfromtxt('pu241_fission.csv', delimiter=",")
#Create The Flux Spectrum (three of them)
DU_Samples = np.genfromtxt('DU_Samples.csv', delimiter=",")
FBR_Blanket = np.genfromtxt('FBR_Blanket.csv', delimiter=",")
chi = lambda E:  0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)  #Fission Spectrum

#Determine atoms per gram for each isotope
Na=0.60221409
N235=3.368235853*Na*10**-9
N238=1170.465962*Na*10**-9
N239=18.14187504*Na*10**-9
N240=1.507603392*Na*10**-9
N241=0.624170052*Na*10**-9



#make interpolation functions
from scipy import interpolate
#Interpolation of cross sections
sig_f_235_interp = interpolate.interp1d(sigma_f_235[:,0], sigma_f_235[:,1],bounds_error=False, fill_value=sigma_f_235[-1,1])
sig_f_238_interp = interpolate.interp1d(sigma_f_238[:,0], sigma_f_238[:,1],bounds_error=False, fill_value=sigma_f_238[-1,1])
sig_f_239_interp = interpolate.interp1d(sigma_f_239[:,0], sigma_f_239[:,1],bounds_error=False, fill_value=sigma_f_239[-1,1])
sig_f_240_interp = interpolate.interp1d(sigma_f_240[:,0], sigma_f_240[:,1],bounds_error=False, fill_value=sigma_f_240[-1,1])
sig_f_241_interp = interpolate.interp1d(sigma_f_241[:,0], sigma_f_241[:,1],bounds_error=False, fill_value=sigma_f_241[-1,1])
#Interpolation of flux spectrum
DU_Samples_interp = interpolate.interp1d(DU_Samples[:,0], DU_Samples[:,4],bounds_error=False, fill_value=DU_Samples[-1,1])
FBR_Blanket_interp = interpolate.interp1d(FBR_Blanket[:,0], FBR_Blanket[:,4],bounds_error=False, fill_value=DU_Samples[-1,1])

#get the union of the energy grids
energies = np.union1d(sigma_f_235[:,0], sigma_f_238[:,0])
energies_new = np.union1d(energies,sigma_f_239[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_f_240[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_f_241[:,0])

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
#Function with FBR spectrum
X_f_235_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

#Perform the integration for each group
Groups=[4e-7,0.1,0.25,1,20]
index=np.zeros(len(Groups)+1)

for jjj in range(1,len(Groups)+1):
	for xx in range(0,len(energies)):
		if(energies[xx]<=Groups[jjj-1]):
			index[jjj]=xx

for i in range(0,len(Groups)):
	phi_int_g1=integrate.trapz(chi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	X_f_235_g1=integrate.trapz(X_f_235_phi(energies[0:index]),energies[0:index])
	X_f_238_g1=integrate.trapz(X_f_238_phi(energies[0:index]),energies[0:index])
	X_f_239_g1=integrate.trapz(X_f_239_phi(energies[0:index]),energies[0:index])
	X_f_240_g1=integrate.trapz(X_f_240_phi(energies[0:index]),energies[0:index])
	X_f_241_g1=integrate.trapz(X_f_241_phi(energies[0:index]),energies[0:index])
	
	phi_int_g12=integrate.trapz(DU_Samples_interp(energies[0:index]),energies[0:index])
	X_f_235_g12=integrate.trapz(X_f_235_phi2(energies[0:index]),energies[0:index])
	X_f_238_g12=integrate.trapz(X_f_238_phi2(energies[0:index]),energies[0:index])
	X_f_239_g12=integrate.trapz(X_f_239_phi2(energies[0:index]),energies[0:index])
	X_f_240_g12=integrate.trapz(X_f_240_phi2(energies[0:index]),energies[0:index])
	X_f_241_g12=integrate.trapz(X_f_241_phi2(energies[0:index]),energies[0:index])
	
	phi_int_g13=integrate.trapz(FBR_Blanket_interp(energies[0:index]),energies[0:index])
	X_f_235_g13=integrate.trapz(X_f_235_phi3(energies[0:index]),energies[0:index])
	X_f_238_g13=integrate.trapz(X_f_238_phi3(energies[0:index]),energies[0:index])
	X_f_239_g13=integrate.trapz(X_f_239_phi3(energies[0:index]),energies[0:index])
	X_f_240_g13=integrate.trapz(X_f_240_phi3(energies[0:index]),energies[0:index])
	X_f_241_g13=integrate.trapz(X_f_241_phi3(energies[0:index]),energies[0:index])







#Perform the integration For Group 1 [0 to 0.4 ev]
EE=4e-7
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index=xx
phi_int_g1=integrate.trapz(chi(energies[0:index]),energies[0:index])
X_f_235_g1=integrate.trapz(X_f_235_phi(energies[0:index]),energies[0:index])
X_f_238_g1=integrate.trapz(X_f_238_phi(energies[0:index]),energies[0:index])
X_f_239_g1=integrate.trapz(X_f_239_phi(energies[0:index]),energies[0:index])
X_f_240_g1=integrate.trapz(X_f_240_phi(energies[0:index]),energies[0:index])
X_f_241_g1=integrate.trapz(X_f_241_phi(energies[0:index]),energies[0:index])

phi_int_g12=integrate.trapz(DU_Samples_interp(energies[0:index]),energies[0:index])
X_f_235_g12=integrate.trapz(X_f_235_phi2(energies[0:index]),energies[0:index])
X_f_238_g12=integrate.trapz(X_f_238_phi2(energies[0:index]),energies[0:index])
X_f_239_g12=integrate.trapz(X_f_239_phi2(energies[0:index]),energies[0:index])
X_f_240_g12=integrate.trapz(X_f_240_phi2(energies[0:index]),energies[0:index])
X_f_241_g12=integrate.trapz(X_f_241_phi2(energies[0:index]),energies[0:index])

phi_int_g13=integrate.trapz(FBR_Blanket_interp(energies[0:index]),energies[0:index])
X_f_235_g13=integrate.trapz(X_f_235_phi3(energies[0:index]),energies[0:index])
X_f_238_g13=integrate.trapz(X_f_238_phi3(energies[0:index]),energies[0:index])
X_f_239_g13=integrate.trapz(X_f_239_phi3(energies[0:index]),energies[0:index])
X_f_240_g13=integrate.trapz(X_f_240_phi3(energies[0:index]),energies[0:index])
X_f_241_g13=integrate.trapz(X_f_241_phi3(energies[0:index]),energies[0:index])

if (sys.argv[1]=="t"):
	print ("Group 1")
#	print ("Fission Cross Section")
#	print (X_f_235_g1/phi_int_g1)
#	print (X_f_238_g1/phi_int_g1)
#	print (X_f_239_g1/phi_int_g1)
#	print (X_f_240_g1/phi_int_g1)
#	print (X_f_241_g1/phi_int_g1)
#	
#	print (" ")
	print ("DUO2 Samples")
	print ("Flux is ", phi_int_g12)
	print ("235 ", X_f_235_g12/phi_int_g12)
	print ("238 ", X_f_238_g12/phi_int_g12)
	print ("239 ", X_f_239_g12/phi_int_g12)
	print ("240 ", X_f_240_g12/phi_int_g12)
	print ("241 ", X_f_241_g12/phi_int_g12)
	
#	print (" ")
#	print ("FBR Blanket")
#	print (X_f_235_g13/phi_int_g13)
#	print (X_f_238_g13/phi_int_g13)
#	print (X_f_239_g13/phi_int_g13)
#	print (X_f_240_g13/phi_int_g13)
#	print (X_f_241_g13/phi_int_g13)


#Perform the integration for group 2
EE=0.1
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index2=xx
phi_int_g2=integrate.trapz(chi(energies[index:index2]),energies[index:index2])
X_f_235_g2=integrate.trapz(X_f_235_phi(energies[index:index2]),energies[index:index2])
X_f_238_g2=integrate.trapz(X_f_238_phi(energies[index:index2]),energies[index:index2])
X_f_239_g2=integrate.trapz(X_f_239_phi(energies[index:index2]),energies[index:index2])
X_f_240_g2=integrate.trapz(X_f_240_phi(energies[index:index2]),energies[index:index2])
X_f_241_g2=integrate.trapz(X_f_241_phi(energies[index:index2]),energies[index:index2])

phi_int_g22=integrate.trapz(DU_Samples_interp(energies[index:index2]),energies[index:index2])
X_f_235_g22=integrate.trapz(X_f_235_phi2(energies[index:index2]),energies[index:index2])
X_f_238_g22=integrate.trapz(X_f_238_phi2(energies[index:index2]),energies[index:index2])
X_f_239_g22=integrate.trapz(X_f_239_phi2(energies[index:index2]),energies[index:index2])
X_f_240_g22=integrate.trapz(X_f_240_phi2(energies[index:index2]),energies[index:index2])
X_f_241_g22=integrate.trapz(X_f_241_phi2(energies[index:index2]),energies[index:index2])

phi_int_g23=integrate.trapz(FBR_Blanket_interp(energies[index:index2]),energies[index:index2])
X_f_235_g23=integrate.trapz(X_f_235_phi3(energies[index:index2]),energies[index:index2])
X_f_238_g23=integrate.trapz(X_f_238_phi3(energies[index:index2]),energies[index:index2])
X_f_239_g23=integrate.trapz(X_f_239_phi3(energies[index:index2]),energies[index:index2])
X_f_240_g23=integrate.trapz(X_f_240_phi3(energies[index:index2]),energies[index:index2])
X_f_241_g23=integrate.trapz(X_f_241_phi3(energies[index:index2]),energies[index:index2])

if (sys.argv[1]=="t"):
	print ("Group 2")
#	print ("Fission Cross Section")
#	print (X_f_235_g2/phi_int_g2)
#	print (X_f_238_g2/phi_int_g2)
#	print (X_f_239_g2/phi_int_g2)
#	print (X_f_240_g2/phi_int_g2)
#	print (X_f_241_g2/phi_int_g2)
#	
#	print (" ")
	print ("DUO2 Samples")
	print ("Flux is ", phi_int_g22)
	print ("235 ", X_f_235_g22/phi_int_g22)
	print ("238 ", X_f_238_g22/phi_int_g22)
	print ("239 ", X_f_239_g22/phi_int_g22)
	print ("240 ", X_f_240_g22/phi_int_g22)
	print ("241 ", X_f_241_g22/phi_int_g22)
	
#	print (" ")
#	print ("FBR Blanket")
#	print (X_f_235_g23/phi_int_g23)
#	print (X_f_238_g23/phi_int_g23)
#	print (X_f_239_g23/phi_int_g23)
#	print (X_f_240_g23/phi_int_g23)
#	print (X_f_241_g23/phi_int_g23)

#Perform the integration for group 3
EE=0.25
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index3=xx
phi_int_g3=integrate.trapz(chi(energies[index2:index3]),energies[index2:index3])
X_f_235_g3=integrate.trapz(X_f_235_phi(energies[index2:index3]),energies[index2:index3])
X_f_238_g3=integrate.trapz(X_f_238_phi(energies[index2:index3]),energies[index2:index3])
X_f_239_g3=integrate.trapz(X_f_239_phi(energies[index2:index3]),energies[index2:index3])
X_f_240_g3=integrate.trapz(X_f_240_phi(energies[index2:index3]),energies[index2:index3])
X_f_241_g3=integrate.trapz(X_f_241_phi(energies[index2:index3]),energies[index2:index3])
         
phi_int_g32=integrate.trapz(DU_Samples_interp(energies[index2:index3]),energies[index2:index3])
X_f_235_g32=integrate.trapz(X_f_235_phi2(energies[index2:index3]),energies[index2:index3])
X_f_238_g32=integrate.trapz(X_f_238_phi2(energies[index2:index3]),energies[index2:index3])
X_f_239_g32=integrate.trapz(X_f_239_phi2(energies[index2:index3]),energies[index2:index3])
X_f_240_g32=integrate.trapz(X_f_240_phi2(energies[index2:index3]),energies[index2:index3])
X_f_241_g32=integrate.trapz(X_f_241_phi2(energies[index2:index3]),energies[index2:index3])
         
phi_int_g33=integrate.trapz(FBR_Blanket_interp(energies[index2:index3]),energies[index2:index3])
X_f_235_g33=integrate.trapz(X_f_235_phi3(energies[index2:index3]),energies[index2:index3])
X_f_238_g33=integrate.trapz(X_f_238_phi3(energies[index2:index3]),energies[index2:index3])
X_f_239_g33=integrate.trapz(X_f_239_phi3(energies[index2:index3]),energies[index2:index3])
X_f_240_g33=integrate.trapz(X_f_240_phi3(energies[index2:index3]),energies[index2:index3])
X_f_241_g33=integrate.trapz(X_f_241_phi3(energies[index2:index3]),energies[index2:index3])

if (sys.argv[1]=="t"):
	print ("Group 3")
#	print ("Fission Cross Section")
#	print (X_f_235_g3/phi_int_g3)
#	print (X_f_238_g3/phi_int_g3)
#	print (X_f_239_g3/phi_int_g3)
#	print (X_f_240_g3/phi_int_g3)
#	print (X_f_241_g3/phi_int_g3)
#	
#	print (" ")
	print ("DUO2 Samples")
	print ("Flux is ", phi_int_g32)
	print ("235 ", X_f_235_g32/phi_int_g32)
	print ("238 ", X_f_238_g32/phi_int_g32)
	print ("239 ", X_f_239_g32/phi_int_g32)
	print ("240 ", X_f_240_g32/phi_int_g32)
	print ("241 ", X_f_241_g32/phi_int_g32)
#	
#	print (" ")
#	print ("FBR Blanket")
#	print (X_f_235_g33/phi_int_g33)
#	print (X_f_238_g33/phi_int_g33)
#	print (X_f_239_g33/phi_int_g33)
#	print (X_f_240_g33/phi_int_g33)
#	print (X_f_241_g33/phi_int_g33)

#Perform the integration for group 4
EE=1
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index4=xx
phi_int_g4=integrate.trapz(chi(energies[index3:index4]),energies[index3:index4])
X_f_235_g4=integrate.trapz(X_f_235_phi(energies[index3:index4]),energies[index3:index4])
X_f_238_g4=integrate.trapz(X_f_238_phi(energies[index3:index4]),energies[index3:index4])
X_f_239_g4=integrate.trapz(X_f_239_phi(energies[index3:index4]),energies[index3:index4])
X_f_240_g4=integrate.trapz(X_f_240_phi(energies[index3:index4]),energies[index3:index4])
X_f_241_g4=integrate.trapz(X_f_241_phi(energies[index3:index4]),energies[index3:index4])
         
phi_int_g42=integrate.trapz(DU_Samples_interp(energies[index3:index4]),energies[index3:index4])
X_f_235_g42=integrate.trapz(X_f_235_phi2(energies[index3:index4]),energies[index3:index4])
X_f_238_g42=integrate.trapz(X_f_238_phi2(energies[index3:index4]),energies[index3:index4])
X_f_239_g42=integrate.trapz(X_f_239_phi2(energies[index3:index4]),energies[index3:index4])
X_f_240_g42=integrate.trapz(X_f_240_phi2(energies[index3:index4]),energies[index3:index4])
X_f_241_g42=integrate.trapz(X_f_241_phi2(energies[index3:index4]),energies[index3:index4])
         
phi_int_g43=integrate.trapz(FBR_Blanket_interp(energies[index3:index4]),energies[index3:index4])
X_f_235_g43=integrate.trapz(X_f_235_phi3(energies[index3:index4]),energies[index3:index4])
X_f_238_g43=integrate.trapz(X_f_238_phi3(energies[index3:index4]),energies[index3:index4])
X_f_239_g43=integrate.trapz(X_f_239_phi3(energies[index3:index4]),energies[index3:index4])
X_f_240_g43=integrate.trapz(X_f_240_phi3(energies[index3:index4]),energies[index3:index4])
X_f_241_g43=integrate.trapz(X_f_241_phi3(energies[index3:index4]),energies[index3:index4])

if (sys.argv[1]=="t"):
	print ("Group 4")
#	print ("Fission Cross Section")
#	print (X_f_235_g4/phi_int_g4)
#	print (X_f_238_g4/phi_int_g4)
#	print (X_f_239_g4/phi_int_g4)
#	print (X_f_240_g4/phi_int_g4)
#	print (X_f_241_g4/phi_int_g4)
#	
#	print (" ")
	print ("DUO2 Samples")
	print ("Flux is ", phi_int_g42)
	print ("235 ", X_f_235_g42/phi_int_g42)
	print ("238 ", X_f_238_g42/phi_int_g42)
	print ("239 ", X_f_239_g42/phi_int_g42)
	print ("240 ", X_f_240_g42/phi_int_g42)
	print ("241 ", X_f_241_g42/phi_int_g42)
#	                            
#	print (" ")                 
#	print ("FBR Blanket")       
#	print (X_f_235_g43/phi_int_g43)
#	print (X_f_238_g43/phi_int_g43)
#	print (X_f_239_g43/phi_int_g43)
#	print (X_f_240_g43/phi_int_g43)
#	print (X_f_241_g43/phi_int_g43)

	
#Perform the integration for group 5
EE=20
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index5=xx
phi_int_g5=integrate.trapz(chi(energies[index4:index5]),energies[index4:index5])
X_f_235_g5=integrate.trapz(X_f_235_phi(energies[index4:index5]),energies[index4:index5])
X_f_238_g5=integrate.trapz(X_f_238_phi(energies[index4:index5]),energies[index4:index5])
X_f_239_g5=integrate.trapz(X_f_239_phi(energies[index4:index5]),energies[index4:index5])
X_f_240_g5=integrate.trapz(X_f_240_phi(energies[index4:index5]),energies[index4:index5])
X_f_241_g5=integrate.trapz(X_f_241_phi(energies[index4:index5]),energies[index4:index5])
         
phi_int_g52=integrate.trapz(DU_Samples_interp(energies[index4:index5]),energies[index4:index5])
X_f_235_g52=integrate.trapz(X_f_235_phi2(energies[index4:index5]),energies[index4:index5])
X_f_238_g52=integrate.trapz(X_f_238_phi2(energies[index4:index5]),energies[index4:index5])
X_f_239_g52=integrate.trapz(X_f_239_phi2(energies[index4:index5]),energies[index4:index5])
X_f_240_g52=integrate.trapz(X_f_240_phi2(energies[index4:index5]),energies[index4:index5])
X_f_241_g52=integrate.trapz(X_f_241_phi2(energies[index4:index5]),energies[index4:index5])
         
phi_int_g53=integrate.trapz(FBR_Blanket_interp(energies[index4:index5]),energies[index4:index5])
X_f_235_g53=integrate.trapz(X_f_235_phi3(energies[index4:index5]),energies[index4:index5])
X_f_238_g53=integrate.trapz(X_f_238_phi3(energies[index4:index5]),energies[index4:index5])
X_f_239_g53=integrate.trapz(X_f_239_phi3(energies[index4:index5]),energies[index4:index5])
X_f_240_g53=integrate.trapz(X_f_240_phi3(energies[index4:index5]),energies[index4:index5])
X_f_241_g53=integrate.trapz(X_f_241_phi3(energies[index4:index5]),energies[index4:index5])

if (sys.argv[1]=="t"):
	print ("Group 5")
#	print ("Fission Cross Section")
#	print (X_f_235_g5/phi_int_g5)
#	print (X_f_238_g5/phi_int_g5)
#	print (X_f_239_g5/phi_int_g5)
#	print (X_f_240_g5/phi_int_g5)
#	print (X_f_241_g5/phi_int_g5)
#	
#	print (" ")
	print ("DUO2 Samples")
	print ("Flux is ", phi_int_g52)
	print ("235 ", X_f_235_g52/phi_int_g52)
	print ("238 ", X_f_238_g52/phi_int_g52)
	print ("239 ", X_f_239_g52/phi_int_g52)
	print ("240 ", X_f_240_g52/phi_int_g52)
	print ("241 ", X_f_241_g52/phi_int_g52)
#	                            
#	print (" ")                 
#	print ("FBR Blanket")       
#	print (X_f_235_g53/phi_int_g53)
#	print (X_f_238_g53/phi_int_g53)
#	print (X_f_239_g53/phi_int_g53)
#	print (X_f_240_g53/phi_int_g53)
#	print (X_f_241_g53/phi_int_g53)
	
if (sys.argv[2]=='t'):
	#Plot Fission Spectrum with cross sections 
	#fig = plt.figure(figsize=(8,6), dpi=1600)
	fig,ax1=plt.subplots()
	#Plot the spectrum
	#ax1.semilogx(energies,chi(energies)/np.sum(chi(energies)),'b--',linewidth=1,label="Fission Spectrum",rasterized=True)
	ax1.semilogx(energies,DU_Samples_interp(energies)/np.sum(DU_Samples_interp(energies)),'b',linewidth=2,label="DU_Samples",rasterized=True)
	#ax1.semilogx(energies,FBR_Blanket_interp(energies)/np.sum(FBR_Blanket_interp(energies)),linewidth=1,label="FBR_Blanket",rasterized=True)
	#ax1.set_legend(loc=3) #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#Plot the energy bins
	ax1.semilogx(4e-7*(energies/energies),np.linspace(0,0.0017,len(energies)),'m',linewidth=2,rasterized=True)
	ax1.semilogx(0.1*(energies/energies),np.linspace(0,0.0017,len(energies)),'m',linewidth=2,rasterized=True)
	ax1.semilogx(0.25*(energies/energies),np.linspace(0,0.0017,len(energies)),'m',linewidth=2,rasterized=True)
	ax1.semilogx(1*(energies/energies),np.linspace(0,0.0017,len(energies)),'m',linewidth=2,rasterized=True)
	#label some stuff
	ax1.set_xlabel("E (MeV)")
	ax1.set_ylabel("Flux Normalized (n/cm$^{2}$s)",color='b')
	ax1.legend(loc=1)
	for t1 in ax1.get_yticklabels():
		t1.set_color('b')
	ax2=ax1.twinx()
	#Plot the Cross Sections
	ax2.loglog(energies, sig_f_235_interp(energies),'r',label="$\sigma_\mathrm{fU235}$",rasterized=True)
	ax2.loglog(energies, sig_f_238_interp(energies),'c',label="$\sigma_\mathrm{fU238}$",rasterized=True)
	ax2.loglog(energies, sig_f_239_interp(energies),'g',label="$\sigma_\mathrm{fPu239}$",rasterized=True)
	ax2.loglog(energies, sig_f_240_interp(energies),'k',label="$\sigma_\mathrm{fPu240}$",rasterized=True)
	ax2.loglog(energies, sig_f_241_interp(energies),'y',label="$\sigma_\mathrm{fPu241}$",rasterized=True)
	ax2.set_ylabel("$\sigma$ (barns)")
	
	ax2.legend(loc=3)

	plt.savefig("Cross_Section_Over_Spectrum_Energy_Bins.pdf")
	plt.clf()

