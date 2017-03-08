#!/usr/bin/env python3

"""
My functions...to clean up the main code
"""

__author__     =  "Paul Mendoza"
__copyright__  =  "Copyright 2016, Planet Earth"
__credits__    = ["Ryan_McClarren"]
__license__    =  "GPL"
__version__    =  "1.0.1"
__maintainer__ =  "Paul Mendoza"
__email__      =  "paul.m.mendoza@gmail.com"
__status__     =  "Production"

################################################################
##################### Import packages ##########################
################################################################

import scipy.special as sps
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "monospace"
import matplotlib
matplotlib.rc('text',usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import random as rn
import matplotlib.mlab as mlab
import copy
import os
from scipy import interpolate
from scipy import integrate
from scipy.integrate import trapz

#############################################################
######################### Variables #########################
#############################################################

# Basic information
FigureSize = (11, 6)              # Dimensions of the figure
TypeOfFamily='monospace'          # This sets the type of font for text
font = {'family' : TypeOfFamily}  # This sets the type of font for text
LegendFontSize = 12
Lfont = {'family' : TypeOfFamily}  # This sets up legend font
Lfont['size']=LegendFontSize

Title = ''
TitleFontSize = 22
TitleFontWeight = "bold"  # "bold" or "normal"

Xlabel='E (eV)'   # X label
XFontSize=18          # X label font size
XFontWeight="normal"  # "bold" or "normal"
XScale="log"       # 'linear' or 'log'
XScale='log'

Xlimits=False                                                # Set xlimits?
XLim=[10**-5,10**8]                                            # Limits that will be set 

Ylimits=True                                         # Set Ylimits?
YLim=[0,0.4]                                     # Limits that will be set  

YFontSize=18                    # Y label font size
YFontWeight="normal"            # "bold" or "normal"
YScale="log"                 # 'linear' or 'log'
YScale='linear'

###################################################################################################
######################################## Markers ##################################################
###################################################################################################

# These are some colors that I found that are distinct (the last two are repeats)
# For coloring the points on the plot, these colors will be used
# http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
Colors=["aqua","gray","red","blue","black","green","magenta","indigo","lime","peru","steelblue",
                "darkorange","salmon","yellow","lime","black"]

# If you want to highlight a specific item, set its alpha value =1 and all others to 0.4
# You can also change the MarkSize (or just use the highlight option below)
Alpha_Value=[1  ,1  ,1  ,1  ,1  ,1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]
MarkSize=   [8  ,8  ,8  ,8  ,8  ,8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8]


Linewidth=[1  ,1  ,1  ,1  ,1  ,1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]

MarkerType=["8","s","p","D","*","H","h","d","^",">"] # Can change all these to "." or "" for nothing "x" isn't that good
# "*" "H" "h" "d" "^" ">" good ones
# http://matplotlib.org/1.4.1/api/markers_api.html for more options

# Set the line styles
# More options: http://matplotlib.org/api/lines_api.html#matplotlib.lines.Line2D.set_linestyle
# LineStyles=["solid","dashed","dash_dot","dotted","."]
LineStyles=["."]

LineWidth=[2,2,2,2]

################################################################
######################### Functions ############################
################################################################

def flux(E,Emt,Eme,E0,Ef):
    """
    Feed in energy as a function of ev not MeV
    """
    #Set constants
    C1=(E0**2)/(Emt**2)*np.exp(Emt/E0)
    C2=1
    C3=(Ef/Eme)*np.exp(Eme/Ef)*1/(np.sqrt(Eme/Ef))
    #Initialize flux
    F=np.zeros(len(E))
    FdE=np.zeros(len(E)-1)
    #Loop through all values of E and calculate Flux
    for i in range(0,len(E)):       
        if E[i]<=Emt:
            F[i]=C1*(E[i]/(E0**2))*np.exp(-E[i]/E0)
        if Emt<E[i] and E[i]<=Eme:
            F[i]=C2/E[i]
        if E[i]>Eme:
            F[i]=C3*(np.sqrt(E[i]/Ef)/Ef)*np.exp(-E[i]/Ef)
    
    return(F)

#def chi(E,C=0.4865,a=1,b=2):
def chi(E,C=0.4865,a=1,b=2):
    """
    Function for producing a chi spectrum (Energy is a function of MeV not ev)
    """
    #a ranges from 0.966 - 1.17 (according to MCNP)
    #b ranges from 1.4610 - 3.4005
    #Initialize Flux
    F=np.zeros(len(E))
    for i in range(0,len(E)):
        F[i]=C*np.exp(-E[i]/a)*np.sinh(np.sqrt(b*E[i]))
    return(F)

def PercentChi(C,a,b,E,Ethermal,Eepi,Efast):
    F=chi(E*10**-6,C,a,b)
    Phi_int=integrate.trapz(F,E*10**-6)


def VaryChi(E,C,a,b):
    """
    This function will produce a plot with varied Chi Values
    """
    fig=plt.figure(figsize=FigureSize)
    ax=fig.add_subplot(111)

    Ylabel='$\phi$(E)$\cdot$E Normalized (n/cm$^{2}$s)'    # Y label
    Ylabel='$\phi$ Normalized (n/cm$^{2}$s)'    # Y label

    Check=0
    for ci in C:
        for ai in a:
            for bi in b:
                F=chi(E*10**-6,ci,ai,bi) #Plug values into function
                #Perform the integral for Flux(E)
                Phi_int=integrate.trapz(F,E*10**-6) #Chi
                #Normalized phi
                Fnorm=F/Phi_int
                #This should be one
                Phi_int2=integrate.trapz(Fnorm,E*10**-6)
                #print values
                print("C=",ci,"a=",ai,"b=",bi,"Int=",Phi_int,";Int2=",Phi_int2)
                #Plot values
                label="C="+str(ci)+";a="+str(ai)+";b="+str(bi)+";I="+str(Phi_int)
                #(fig,ax)=f.plot(E,(Fnorm*E),ax,Check,label,fig,Ylabel,'log','linear',MarkerType,[1])
                (fig,ax)=plot(E,Fnorm,ax,Check,label,fig,Ylabel,'log','linear',MarkerType,[1])
                Check=Check+1
    ax=Legend(ax)
    plt.savefig("Figures/Flux_Spectra_Chi_Vary.pdf")

def VaryPhi(E,Emt,Eme,E0,Ef):
    """
    This function will produce a plot with varied phi values
    """
    fig=plt.figure(figsize=FigureSize)
    ax=fig.add_subplot(111)
    
    Ylabel='$\phi$(E)$\cdot$E Normalized (n/cm$^{2}$s)'    # Y label
    Ylabel='$\phi$ Normalized (n/cm$^{2}$s)'    # Y label

    Check=0
    for Emti in Emt:
        for Emei in Eme:
            for E0i in E0:
                for Efi in Ef:
                    #Calculate flux (yes we need E)
                    F=flux(E,Emti,Emei,E0i,Efi)
                    #Perform the integral for Flux(E)
                    Phi_int=integrate.trapz(F,E) #Chi
                    #Normalized phi
                    Fnorm=F/Phi_int
                    #This should be one
                    Phi_int2=integrate.trapz(Fnorm,E)
                    #print values
                    print("Emt=",Emti,"Eme=",Emei,"E0=",E0i,"Ef",Efi,"Int=",Phi_int,";Int2=",Phi_int2)
                    #Plot values
                    label="Emt="+str(Emti)+"Eme="+str(Emei)+"E0="+str(E0i)+"Ef"+str(Efi)+\
                    "Int="+str(Phi_int)
                    #(fig,ax)=f.plot(E,(Fnorm*E),ax,Check,label,fig,Ylabel,'log','log',MarkerType,[1])
                    (fig,ax)=plot(E,Fnorm,ax,Check,label,fig,Ylabel,'log','log',MarkerType,[1])
                    Check=Check+1
    ax=Legend(ax)
    plt.savefig("Figures/Flux_Spectra_Phi_Vary.pdf")

def PlotChiNPhi(E,Emt,Eme,E0,Ef,C,a,b):
    """
    This function plots both ChiNPhi on the same plot
    """
    fig=plt.figure(figsize=FigureSize)
    ax=fig.add_subplot(111)
    
    Ylabel='$\phi$(E)$\cdot$E Normalized (n/cm$^{2}$s)'    # Y label
    Ylabel='$\phi$ Normalized (n/cm$^{2}$s)'    # Y label
    
    #Calculate flux (yes we need E)
    FPhi=flux(E,Emt,Eme,E0,Ef)
    #Perform the integral for Flux(E)
    Phi_int=integrate.trapz(FPhi,E) 
    #Normalized phi
    FnormPhi=FPhi/Phi_int
    #This should be one
    Phi_int2=integrate.trapz(FnormPhi,E)
    #Plot values
    label="Phi"
    #(fig,ax)=f.plot(E,(Fnorm*E),ax,0,label,fig,Ylabel,'log','log',[""],LineWidth)
    (fig,ax)=plot(E,FnormPhi,ax,0,label,fig,Ylabel,'log','log',[""],LineWidth)


    Fchi=chi(E*10**-6,C,a,b) #Plug values into function
    #Perform the integral for Flux(E)
    Phi_int=integrate.trapz(Fchi,E*10**-6) #Chi
    #Normalized phi
    Fnormchi=Fchi/Phi_int
    #This should be one
    Phi_int2=integrate.trapz(Fnormchi,E*10**-6)
    #Plot values
    label="Chi"
    #(fig,ax)=f.plot(E,(Fnorm*E),ax,1,label,fig,Ylabel,'log','log',[""],LineWidth)
    (fig,ax)=plot(E,Fnormchi,ax,1,label,fig,Ylabel,'log','log',[""],LineWidth)

                
    ax=Legend(ax)
    plt.savefig("Figures/Phi_N_Chi.pdf")

    
    
def loop_values(list1,index):
    """                                                                                               
    This function will loop through values in list even if outside range 
    (in the positive sense not negative)                                                                                                  
    """
    while True:
        try:
            list1[index]
            break
        except IndexError:
            index=index-len(list1)
    return(list1[index])                                    

def plot(x,y,ax,Check,label,fig,Ylabel,XScale,YScale,Markers,LineWidth):
    Color=loop_values(Colors,Check)
    Marker=loop_values(Markers,Check)
    LineW=loop_values(LineWidth,Check)
    #Plot X and Y
    ax.plot(x,y,
            linestyle="solid", #"solid","dashed","dash_dot","dotted","."
            #marker="", # "*" "H" "h" "d" "^" ">"
            marker=Marker,
# good ones http://matplotlib.org/1.4.1/api/markers_api.html for more
            color=Color,
            markersize=8,
            lw=LineW,
            alpha=1,
            label=label)
    	
    #Log or linear scale?
    ax.set_xscale(XScale)
    ax.set_yscale(YScale)
    #Set Title
    fig.suptitle(Title,fontsize=TitleFontSize,
        	 fontweight=TitleFontWeight,fontdict=font,
                 ha='center')
    #Set X and y labels
    ax.set_xlabel(Xlabel,
        	  fontsize=XFontSize,fontweight=XFontWeight,
        	  fontdict=font)
    ax.set_ylabel(Ylabel,
        	  fontsize=YFontSize,
        	  fontweight=YFontWeight,
        	  fontdict=font)


    #Earlier Era
    # ax[j].plot(a[start:end+1,var.XValues],a[start:end+1,var.YValues2],
    #                                   linestyle=var.loop_values(var.LineStyles,check),
    #                                   marker=var.loop_values(var.MarkerType,check),
    #                                   color=var.loop_values(var.Colors,check),
    #                                   markersize=var.loop_values(var.MarkSize,check)*0.5,
    #                                   alpha=var.loop_values(var.Alpha_Value,check)*0.5,
    #                                   label=unique_names[check],linewidth=var.Linewidth)
            


    
    #Sets up limits of graphs
    if Xlimits:
        ax.set_xlim(XLim[0],XLim[1])
    if Ylimits:
        ax.set_ylim(YLim[0],YLim[1])
                                        
    
    return(fig,ax)

def Legend(ax):
	handles,labels=ax.get_legend_handles_labels()
	ax.legend(handles,labels,loc='best',
			fontsize=LegendFontSize,prop=font)
	return(ax)

def GETcsvFiles(directory):
    """
    This function gathers all files
    ending with ".csv" in a certain directory
    Note...NOT ".CSV" capitalization matters!
    """
    Filelist=[]
    for file in os.listdir(directory):
        if ".csv" in file:
            Filelist.append(file)
    return(Filelist)

def LoopTAPE(Protons,Isotope,Reaction):
    ZAID=Protons+Isotope+"0"
    with open('../Origen2/TAPE9_BANK.inp') as f:
        content=f.readlines()
    if 'a' in Reaction:
        Place=2
    elif 'f' in Reaction:
        Place=5
    
    for i in content:
        hold=i.split()
        if len(hold)>2:
                                   #want 600 libs not 1,2 or 3
            if ZAID in hold[1] and len(hold[0])>1:
                X_Section=hold[Place]
                break
    return(X_Section)

def MinIndex(List):
    Count=0
    List2=[]
    for item in List:
        List2.append(item[0])
    Min=min(List2)
    for item in List:
        if item[0]==Min:
            break
        Count=Count+1
    return(Count)

def DetermineAverages(List):
    Count=0
    Averages=np.zeros(4)
    if len(List[0])==5:
        for item in List:
            for i in range(0,4):
                Averages[i]=Averages[i]+float(item[i])
                Count=Count+1
        Averages=Averages/Count
    else:
        for i in range(0,4):
            Averages[i]=float(List[i])
    return(Averages)
