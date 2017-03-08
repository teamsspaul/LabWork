import numpy as np
import sys
import matplotlib.pyplot as plt
#open Fission cross-sections
sigma_f_235 = np.genfromtxt('u235_fission.csv', delimiter=",")
sigma_f_238 = np.genfromtxt('u238_fission.csv', delimiter=",")
sigma_f_239 = np.genfromtxt('pu239_fission.csv', delimiter=",")
sigma_f_240 = np.genfromtxt('pu240_fission.csv', delimiter=",")
sigma_f_241 = np.genfromtxt('pu241_fission.csv', delimiter=",")

sigma_tot_147=np.genfromtxt('Nd_147_tot.csv',delimiter=",")
sigma_el_147=np.genfromtxt('Nd_147_el.csv',delimiter=",")
sigma_inel_147=np.genfromtxt('Nd_147_inel.csv',delimiter=",")
sigma_2n_147=np.genfromtxt('Nd_147_2n.csv',delimiter=",")
sigma_a_147=np.genfromtxt('Nd_147_a.csv',delimiter=",")
sigma_a_148=np.genfromtxt('Nd_148_a.csv',delimiter=",")
sigma_a_137=np.genfromtxt('cs_137_a.csv',delimiter=",")

DU_Samples = np.genfromtxt('DU_Samples.csv', delimiter=",")
FBR_Blanket = np.genfromtxt('FBR_Blanket.csv', delimiter=",")



#create the fission spectrum 
chi = lambda E:  0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

#make interpolation functions
from scipy import interpolate
sig_f_235_interp = interpolate.interp1d(sigma_f_235[:,0], sigma_f_235[:,1],bounds_error=False, fill_value=sigma_f_235[-1,1])
sig_f_238_interp = interpolate.interp1d(sigma_f_238[:,0], sigma_f_238[:,1],bounds_error=False, fill_value=sigma_f_238[-1,1])
sig_f_239_interp = interpolate.interp1d(sigma_f_239[:,0], sigma_f_239[:,1],bounds_error=False, fill_value=sigma_f_239[-1,1])
sig_f_240_interp = interpolate.interp1d(sigma_f_240[:,0], sigma_f_240[:,1],bounds_error=False, fill_value=sigma_f_240[-1,1])
sig_f_241_interp = interpolate.interp1d(sigma_f_241[:,0], sigma_f_241[:,1],bounds_error=False, fill_value=sigma_f_241[-1,1])

DU_Samples_interp = interpolate.interp1d(DU_Samples[:,0], DU_Samples[:,4],bounds_error=False, fill_value=DU_Samples[-1,1])
FBR_Blanket_interp = interpolate.interp1d(FBR_Blanket[:,0], FBR_Blanket[:,4],bounds_error=False, fill_value=DU_Samples[-1,1])


sig_tot_147_interp=interpolate.interp1d(sigma_tot_147[:,0],sigma_tot_147[:,1],bounds_error=False,fill_value=0)
sig_el_147_interp=interpolate.interp1d(sigma_el_147[:,0],sigma_el_147[:,1],bounds_error=False,fill_value=0)
sig_inel_147_interp=interpolate.interp1d(sigma_inel_147[:,0],sigma_inel_147[:,1],bounds_error=False,fill_value=0)
sig_2n_147_interp=interpolate.interp1d(sigma_2n_147[:,0],sigma_2n_147[:,1],bounds_error=False,fill_value=0)
sig_a_147_interp=interpolate.interp1d(sigma_a_147[:,0],sigma_a_147[:,1],bounds_error=False,fill_value=0)
sig_a_148_interp=interpolate.interp1d(sigma_a_148[:,0],sigma_a_148[:,1],bounds_error=False,fill_value=0)
sig_a_137_interp=interpolate.interp1d(sigma_a_137[:,0],sigma_a_137[:,1],bounds_error=False,fill_value=0)


#get the union of the energy grids
energies = np.union1d(sigma_f_235[:,0], sigma_f_238[:,0])

#I am not sure what these Commands do; print(len(energies))
energies_new = np.union1d(energies,sigma_f_239[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_f_240[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_f_241[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_tot_147[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_el_147[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_inel_147[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_2n_147[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_a_147[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_a_148[:,0])
energies = energies_new
energies_new = np.union1d(energies,sigma_a_137[:,0])
energies = energies_new
#But Dr. MC did them...

#energies = np.union1d(energies,sigma_t_238[:,0])
#If I want to append energies I Should do this print(len(energies))

sig_ab_147_interp=sig_tot_147_interp(energies)-sig_el_147_interp(energies)-sig_inel_147_interp(energies)-sig_2n_147_interp(energies)


#chi(energies)

#Perform Integration
from scipy import integrate
from scipy.integrate import trapz
#Make Functions of things I want to integrate.
X_f_235_phi=interpolate.interp1d(energies,chi(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi=interpolate.interp1d(energies,chi(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi=interpolate.interp1d(energies,chi(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi=interpolate.interp1d(energies,chi(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi=interpolate.interp1d(energies,chi(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

X_tot_147_phi=interpolate.interp1d(energies,chi(energies)*sig_tot_147_interp(energies),fill_value=0,bounds_error=False)
X_el_147_phi=interpolate.interp1d(energies,chi(energies)*sig_el_147_interp(energies),fill_value=0,bounds_error=False)
X_inel_147_phi=interpolate.interp1d(energies,chi(energies)*sig_inel_147_interp(energies),fill_value=0,bounds_error=False)
X_ab_147_phi=interpolate.interp1d(energies,chi(energies)*sig_ab_147_interp,fill_value=0,bounds_error=False)
X_2n_147_phi=interpolate.interp1d(energies,chi(energies)*sig_2n_147_interp(energies),fill_value=0,bounds_error=False)
X_a_147_phi=interpolate.interp1d(energies,chi(energies)*sig_a_147_interp(energies),fill_value=0,bounds_error=False)
X_a_148_phi=interpolate.interp1d(energies,chi(energies)*sig_a_148_interp(energies),fill_value=0,bounds_error=False)
X_a_137_phi=interpolate.interp1d(energies,chi(energies)*sig_a_137_interp(energies),fill_value=0,bounds_error=False)

X_f_235_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

X_tot_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_tot_147_interp(energies),fill_value=0,bounds_error=False)
X_el_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_el_147_interp(energies),fill_value=0,bounds_error=False)
X_inel_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_inel_147_interp(energies),fill_value=0,bounds_error=False)
X_ab_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_ab_147_interp,fill_value=0,bounds_error=False)
X_2n_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_2n_147_interp(energies),fill_value=0,bounds_error=False)
X_a_147_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_a_147_interp(energies),fill_value=0,bounds_error=False)
X_a_148_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_a_148_interp(energies),fill_value=0,bounds_error=False)
X_a_137_phi2=interpolate.interp1d(energies,DU_Samples_interp(energies)*sig_a_137_interp(energies),fill_value=0,bounds_error=False)

X_f_235_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_235_interp(energies),fill_value=0,bounds_error=False)
X_f_238_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_238_interp(energies),fill_value=0,bounds_error=False)
X_f_239_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_239_interp(energies),fill_value=0,bounds_error=False)
X_f_240_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_240_interp(energies),fill_value=0,bounds_error=False)
X_f_241_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_f_241_interp(energies),fill_value=0,bounds_error=False)

X_tot_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_tot_147_interp(energies),fill_value=0,bounds_error=False)
X_el_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_el_147_interp(energies),fill_value=0,bounds_error=False)
X_inel_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_inel_147_interp(energies),fill_value=0,bounds_error=False)
X_ab_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_ab_147_interp,fill_value=0,bounds_error=False)
X_2n_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_2n_147_interp(energies),fill_value=0,bounds_error=False)
X_a_147_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_a_147_interp(energies),fill_value=0,bounds_error=False)
X_a_148_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_a_148_interp(energies),fill_value=0,bounds_error=False)
X_a_137_phi3=interpolate.interp1d(energies,FBR_Blanket_interp(energies)*sig_a_137_interp(energies),fill_value=0,bounds_error=False)





#Perform the integration For Group 1
EE=30
for xx in range(0, len(energies)):
    if (energies[xx]<=EE):
        index=xx
phi_int_g1=integrate.trapz(chi(energies[0:index]),energies[0:index])
X_f_235_g1=integrate.trapz(X_f_235_phi(energies[0:index]),energies[0:index])
X_f_238_g1=integrate.trapz(X_f_238_phi(energies[0:index]),energies[0:index])
X_f_239_g1=integrate.trapz(X_f_239_phi(energies[0:index]),energies[0:index])
X_f_240_g1=integrate.trapz(X_f_240_phi(energies[0:index]),energies[0:index])
X_f_241_g1=integrate.trapz(X_f_241_phi(energies[0:index]),energies[0:index])

X_tot_147_g1=integrate.trapz(X_tot_147_phi(energies[0:index]),energies[0:index])
X_el_147_g1=integrate.trapz(X_el_147_phi(energies[0:index]),energies[0:index])
X_inel_147_g1=integrate.trapz(X_inel_147_phi(energies[0:index]),energies[0:index])
X_ab_147_g1=integrate.trapz(X_ab_147_phi(energies[0:index]),energies[0:index])
X_2n_147_g1=integrate.trapz(X_2n_147_phi(energies[0:index]),energies[0:index])
X_a_147_g1=integrate.trapz(X_a_147_phi(energies[0:index]),energies[0:index])
X_a_148_g1=integrate.trapz(X_a_148_phi(energies[0:index]),energies[0:index])
X_a_137_g1=integrate.trapz(X_a_137_phi(energies[0:index]),energies[0:index])

phi_int_g12=integrate.trapz(DU_Samples_interp(energies[0:index]),energies[0:index])
X_f_235_g12=integrate.trapz(X_f_235_phi2(energies[0:index]),energies[0:index])
X_f_238_g12=integrate.trapz(X_f_238_phi2(energies[0:index]),energies[0:index])
X_f_239_g12=integrate.trapz(X_f_239_phi2(energies[0:index]),energies[0:index])
X_f_240_g12=integrate.trapz(X_f_240_phi2(energies[0:index]),energies[0:index])
X_f_241_g12=integrate.trapz(X_f_241_phi2(energies[0:index]),energies[0:index])

X_tot_147_g12=integrate.trapz(X_tot_147_phi2(energies[0:index]),energies[0:index])
X_el_147_g12=integrate.trapz(X_el_147_phi2(energies[0:index]),energies[0:index])
X_inel_147_g12=integrate.trapz(X_inel_147_phi2(energies[0:index]),energies[0:index])
X_ab_147_g12=integrate.trapz(X_ab_147_phi2(energies[0:index]),energies[0:index])
X_2n_147_g12=integrate.trapz(X_2n_147_phi2(energies[0:index]),energies[0:index])
X_a_147_g12=integrate.trapz(X_a_147_phi2(energies[0:index]),energies[0:index])
X_a_148_g12=integrate.trapz(X_a_148_phi2(energies[0:index]),energies[0:index])
X_a_137_g12=integrate.trapz(X_a_137_phi2(energies[0:index]),energies[0:index])

phi_int_g13=integrate.trapz(FBR_Blanket_interp(energies[0:index]),energies[0:index])
X_f_235_g13=integrate.trapz(X_f_235_phi3(energies[0:index]),energies[0:index])
X_f_238_g13=integrate.trapz(X_f_238_phi3(energies[0:index]),energies[0:index])
X_f_239_g13=integrate.trapz(X_f_239_phi3(energies[0:index]),energies[0:index])
X_f_240_g13=integrate.trapz(X_f_240_phi3(energies[0:index]),energies[0:index])
X_f_241_g13=integrate.trapz(X_f_241_phi3(energies[0:index]),energies[0:index])

X_tot_147_g13=integrate.trapz(X_tot_147_phi3(energies[0:index]),energies[0:index])
X_el_147_g13=integrate.trapz(X_el_147_phi3(energies[0:index]),energies[0:index])
X_inel_147_g13=integrate.trapz(X_inel_147_phi3(energies[0:index]),energies[0:index])
X_ab_147_g13=integrate.trapz(X_ab_147_phi3(energies[0:index]),energies[0:index])
X_2n_147_g13=integrate.trapz(X_2n_147_phi3(energies[0:index]),energies[0:index])
X_a_147_g13=integrate.trapz(X_a_147_phi3(energies[0:index]),energies[0:index])
X_a_148_g13=integrate.trapz(X_a_148_phi3(energies[0:index]),energies[0:index])
X_a_137_g13=integrate.trapz(X_a_137_phi3(energies[0:index]),energies[0:index])



if (sys.argv[1]=="t"):
	print ("Fission Cross Section")
	print (X_f_235_g1/phi_int_g1)
	print (X_f_238_g1/phi_int_g1)
	print (X_f_239_g1/phi_int_g1)
	print (X_f_240_g1/phi_int_g1)
	print (X_f_241_g1/phi_int_g1)
	
	print (" ")
	print ("DUO2 Samples")
	print (X_f_235_g12/phi_int_g12)
	print (X_f_238_g12/phi_int_g12)
	print (X_f_239_g12/phi_int_g12)
	print (X_f_240_g12/phi_int_g12)
	print (X_f_241_g12/phi_int_g12)
	
	print (" ")
	print ("FBR Blanket")
	print (X_f_235_g13/phi_int_g13)
	print (X_f_238_g13/phi_int_g13)
	print (X_f_239_g13/phi_int_g13)
	print (X_f_240_g13/phi_int_g13)
	print (X_f_241_g13/phi_int_g13)


print (" ")
print ("Absorption Nd147")
#print (X_tot_147_g1/phi_int_g1)
#print (X_el_147_g1/phi_int_g1)
#print (X_inel_147_g1/phi_int_g1)
#print (X_2n_147_g1/phi_int_g1)
#print (X_ab_147_g1/phi_int_g1)
print (X_a_147_g1/phi_int_g1)
print (X_a_147_g12/phi_int_g12)
print (X_a_147_g13/phi_int_g13)

print (" ")
print ("Absorption Nd148")
print (X_a_148_g1/phi_int_g1)
print (X_a_148_g12/phi_int_g12)
print (X_a_148_g13/phi_int_g13)
#print (X_tot_147_g12/phi_int_g12)
#print (X_el_147_g12/phi_int_g12)
#print (X_inel_147_g12/phi_int_g12)
#print (X_2n_147_g12/phi_int_g12)
#print (X_ab_147_g12/phi_int_g12)
#print (X_a_147_g12/phi_int_g12)
#print (X_a_148_g12/phi_int_g12)
print (" ")
print ("Absorption Cs137")
print (X_a_137_g1/phi_int_g1)
print (X_a_137_g12/phi_int_g12)
print (X_a_137_g13/phi_int_g13)



#Plot if I want to
if (sys.argv[1]=="t"):
	fig = plt.figure(figsize=(8,6), dpi=1600)
	plt.loglog(energies, sig_f_235_interp(energies), label="$\sigma_\mathrm{fU235}$",rasterized=True)
	plt.loglog(energies, sig_f_238_interp(energies), label="$\sigma_\mathrm{fU238}$",rasterized=True)
	plt.loglog(energies, sig_f_239_interp(energies), label="$\sigma_\mathrm{fPu239}$",rasterized=True)
	plt.loglog(energies, sig_f_240_interp(energies), label="$\sigma_\mathrm{fPu240}$",rasterized=True)
	plt.loglog(energies, sig_f_241_interp(energies), label="$\sigma_\mathrm{fPu241}$",rasterized=True)
	plt.legend(loc=3) #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.ylabel("$\sigma$ (barns)")
	plt.xlabel("E (MeV)")
	plt.savefig("Fission_Cross_Sections.pdf")
	plt.clf()
	#Plot the fission spectrum
	fig = plt.figure(figsize=(8,6), dpi=1600)
	plt.semilogx(energies,chi(energies),rasterized=True)
	plt.xlabel("E (MeV)")
	plt.ylabel("Probability (MeV$^{-1}$)")
	plt.savefig("chi.pdf")
	plt.clf()
	#Plot the other spectrums
	fig = plt.figure(figsize=(8,6), dpi=1600)
	plt.semilogx(energies,DU_Samples_interp(energies),label="DU_Samples",rasterized=True)
	plt.semilogx(energies,FBR_Blanket_interp(energies),label="FBR_Blanket",rasterized=True)
	plt.xlabel("E (MeV)")
	plt.ylabel("Flux (n/cm$^{2}$s)")
	plt.savefig("Matts_Flux_Magnitude.pdf")
	plt.clf()
	#Plot All Spectrum
	fig = plt.figure(figsize=(8,6), dpi=1600)
	plt.semilogx(energies,DU_Samples_interp(energies)/np.sum(DU_Samples_interp(energies)),label="DU_Samples",rasterized=True)
	plt.semilogx(energies,FBR_Blanket_interp(energies)/np.sum(FBR_Blanket_interp(energies)),label="FBR_Blanket",rasterized=True)
	plt.semilogx(energies,chi(energies)/np.sum(chi(energies)),label="Fission Spectrum",rasterized=True)
	plt.legend(loc=3) #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.xlabel("E (MeV)")
	plt.ylabel("Flux Normalized (n/cm$^{2}$s)")
	plt.savefig("Spectrums.pdf")
	plt.clf()
	
if (sys.argv[2]=='t'):
	#Plot Fission Spectrum with cross sections 
	#fig = plt.figure(figsize=(8,6), dpi=1600)
	fig,ax1=plt.subplots()
	ax1.semilogx(energies,chi(energies)/np.sum(chi(energies)),'b-',linewidth=2,label="Fission Spectrum",rasterized=True)
	ax1.semilogx(energies,DU_Samples_interp(energies)/np.sum(DU_Samples_interp(energies)),'b--',linewidth=1,label="DU_Samples",rasterized=True)
	ax1.semilogx(energies,FBR_Blanket_interp(energies)/np.sum(FBR_Blanket_interp(energies)),linewidth=1,label="FBR_Blanket",rasterized=True)
	#ax1.set_legend(loc=3) #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	ax1.set_xlabel("E (MeV)")
	ax1.set_ylabel("Flux Normalized (n/cm$^{2}$s)",color='b')
	ax1.legend(loc=1)
	for t1 in ax1.get_yticklabels():
		t1.set_color('b')
	ax2=ax1.twinx()
	ax2.loglog(energies, sig_f_235_interp(energies),'r',label="$\sigma_\mathrm{fU235}$",rasterized=True)
	ax2.loglog(energies, sig_f_238_interp(energies),'m',label="$\sigma_\mathrm{fU238}$",rasterized=True)
	ax2.loglog(energies, sig_f_238_interp(energies),'c',label="$\sigma_\mathrm{fU238}$",rasterized=True)
	ax2.loglog(energies, sig_f_239_interp(energies),'g',label="$\sigma_\mathrm{fPu239}$",rasterized=True)
	ax2.loglog(energies, sig_f_240_interp(energies),'k',label="$\sigma_\mathrm{fPu240}$",rasterized=True)
	ax2.loglog(energies, sig_f_241_interp(energies),'y',label="$\sigma_\mathrm{fPu241}$",rasterized=True)
	ax2.set_ylabel("$\sigma$ (barns)")
	ax2.legend(loc=3)
	
	plt.savefig("Cross_Section_Over_Spectrum.pdf")
	plt.clf()
	
if (sys.argv[3]=='t'):
	fig = plt.figure(figsize=(8,6), dpi=1600)
	plt.loglog(energies, sig_ab_147_interp, label="$\sigma_\mathrm{aNd147}$",rasterized=True)
	plt.loglog(energies, sig_tot_147_interp(energies), label="$\sigma_\mathrm{tNd147}$",rasterized=True)
	plt.loglog(energies, sig_el_147_interp(energies), label="$\sigma_\mathrm{elNd147}$",rasterized=True)
	plt.loglog(energies, sig_inel_147_interp(energies), label="$\sigma_\mathrm{inelNd147}$",rasterized=True)
	plt.loglog(energies, sig_2n_147_interp(energies), label="$\sigma_\mathrm{2nNd147}$",rasterized=True)
	plt.loglog(energies, sig_a_147_interp(energies), label="$\sigma_\mathrm{a2Nd147}$",rasterized=True)
	plt.legend(loc=3) #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.ylabel("$\sigma$ (barns)")
	plt.xlabel("E (MeV)")
	plt.savefig("Ab_Nd147_Cross_Sections.pdf")
	plt.clf()