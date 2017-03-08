#This Program will Calculate the multi-group integrated flux*cross section values
#Given cross section data, and a flux spectrum.

#Type "python Flux_X_section.py Cross_Section_file.csv
#Make sure the cross section file has MeV first column, and barns second column

import numpy as np
import sys
import matplotlib.pyplot as plt
import math
from scipy import interpolate
from scipy import integrate
from scipy.integrate import trapz

#Open the flux spectrum (Make interp function): Use 9 for DU_Samples its flux(E) - I think
DU_Samples = np.genfromtxt('DU_Samples.csv', delimiter=",")
DU_Samples_interp = interpolate.interp1d(DU_Samples[:,5], DU_Samples[:,9],bounds_error=False, fill_value=DU_Samples[-1,1])

#open cross-section data 
sigma = np.genfromtxt(sys.argv[1], delimiter=",")
#Create an interpolation function (just for fun?)
sig_interp = interpolate.interp1d(sigma[:,0], sigma[:,1],bounds_error=False, fill_value=sigma[-1,1])
#Use all the energies
energies = sigma[:,0]
#Function with DU sample spectrum
X_phi=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_interp(energies),fill_value=0,bounds_error=False)

#These are for gamma+fission = absorption
check=1
try :sys.argv[2]
except IndexError:
	check=0
if (check==1):
	sigma2= np.genfromtxt(sys.argv[2], delimiter=",")
	sig_interp2 = interpolate.interp1d(sigma2[:,0], sigma2[:,1],bounds_error=False, fill_value=sigma[-1,1])
	energies = np.union1d(sigma[:,0], sigma2[:,0])
	X_phi=interpolate.interp1d(energies,DU_Samples_interp(energies)*(sig_interp(energies)+sig_interp2(energies)),fill_value=0,bounds_error=False)

#Groups=np.concatenate((np.logspace(-11,-5,120),np.logspace(-5.001,-1.69897,1000),np.logspace(-1.69898,1.30103,60)),axis=1)
#Groups=[0.001,20] #This is for the percent of "thermal" and fast
Groups=[20] #This is for a single group (but trapz over all the energies)
#Groups=np.logspace(-11,1.30103,10000)

index=np.zeros(len(Groups)+1)

for jjj in range(1,len(Groups)+1):
	for xx in range(0,len(energies)):
		if(energies[xx]<=Groups[jjj-1]):
			index[jjj]=xx


#For the overall fraction
f=np.zeros(1)
#These are for determining individual fission percentages (thermal vs fast)
fi=np.zeros((len(Groups),1))
#These are for flux times cross section
ff=np.zeros(1)

#For the average cross section
X_avg=np.zeros((1,len(Groups)))
#For the average flux 
phi_avg=np.zeros((1,len(Groups)))
phi_int_g=np.zeros((1,len(Groups)))


for i in range(0,len(Groups)):

	phi_int=integrate.trapz(DU_Samples_interp(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])
	phi_int_g[0,i]=phi_int
	#Integrated flux times cross section 
	X_g=integrate.trapz(X_phi(energies[index[i]:index[i+1]]),energies[index[i]:index[i+1]])

	#Store all flux values in one place
	if(energies[index[i+1]]-energies[index[i]]>0):
		phi_avg[0,i]=phi_int/(energies[index[i+1]]-energies[index[i]])

	#Calculate fission fraction and the group cross sections
	f[0]=X_g+f[0]  ### Sum up for total integral
	fi[i,0]=X_g   ### Sum up for part of integral (thermal and fast)
	
	ff[0]=X_g*10**-24+ff[0]  ### To calculate flux times cross section
	
	#Group averaged Cross Section
	if (phi_int!=0):
		X_avg[0,i]=X_g/phi_int

	

	
#Print out the total integral
#print ("Length of Energies",len(energies))
#print ("Total Integral",f[0])
for i in range(0,len(Groups)):#Loop over groups
	print ("Group ", i+1, "Upper Energy ", energies[index[i+1]], "(MeV) Lower Energy ", energies[index[i]], "(MeV)")
	if (phi_int_g[0,i]!=0):
		print ("Flux is ", phi_int_g[0,i], "Flux Average is",phi_avg[0,i])
		print ("Cross Section (b)", X_avg[0,i])
	#print ("Percent of Integral", fi[i,0]/f[0])
	
#print ("flux times cross section Summed over all energy groups",ff[0])