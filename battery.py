# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 23:27:11 2017

@author: joscha
"""
import numpy as np

#################################################################
# Battery #######################################################
    
## Temperature model
  #Temperature 
def TB(Ta, Tb, Pbatloss, h, A, m, cp, ts):
    return Tb+((np.abs(Pbatloss)-h*A*(Tb-Ta))/(cp*m/ts))
        
    ## Stationary battery model
        #Losses
def PBAT(Pt,Tb, Pnom):
    
        
        #ohmic losses for charge or discharge   
    if Pt>0.: #charge 
        etap=-0.02115*(Pt/Pnom)+1
            #losses by temperature
        deltaT= -25+Tb 
        etat= 1-1/3.6*deltaT*(-1.26e-07*Tb**3 - 1.315e-06*Tb**2 + 0.0003748*Tb - 0.006209) # hier sollte nicht 20 stehen sondern der tats. Nennstrom
    
    elif Pt==0.: #idle
        etap=0
        etat=0
            
    elif Pt<0.: #discharge
        etap= -0.02803*(Pt/Pnom)+1
            #losses by temperature
        deltaT= -25+Tb 
        etat= 1-1/3.6*deltaT*(-1.26e-07*Tb**3 - 1.315e-06*Tb**2 + 0.0003748*Tb - 0.006209) # hier sollte nicht 20 stehen sondern der tats. Nennstrom
    
            
    return Pt*etap*etat
        	
        #State of Charge - model    
def SOC(Pbat, ts, SOC, Emax, Pself):
    Pselfabs=Pself*ts*Emax
    return ((Pbat-Pselfabs)*ts/3600+SOC*Emax)/Emax
    
        #effective Capacity
def soceff(T,P,a1,b1,c1,a2,b2, Pref):
    if (c1)>0.:
        socref=b2*Pref+a2
        socT=a1*np.log(b1*(80+T))+c1
        socP=b2*P+a2
    else:
        socref=0
        socT=0
        socP=b2*P+a2
    return socT+socP-socref # note: socref is only needed in dch mode because of temperature AND power dependency
    
    ## Aging model 
    
def QLOSS(P0, P1, Tb0, Tb1, ts, SoC0, SoC1, i, Vnenn=3.6, Cnom=40, Cmax=40):
        
    if P1==0.: #idle
        Qloss = QLOSSCAL(SoC1,Tb1,ts,i, Cnom, Cmax)+Cnom-Cmax 
    else:     
        Qlossrel = QLOSSCYC(P1, Tb1, ts, SoC1, i, Vnenn, Cnom)\
    +((1-Cmax/Cnom)*100)-QLOSSCYC(P1, Tb1, ts, SoC1, (i-1), Vnenn, Cnom) # cumulated cyclic aging
        Qloss=Qlossrel*Cnom/100    
    return Qloss
        
        #Cyclic Aging 
        #gas constante R= 8.3144598 J mol-1 K-1 
        #http://physics.nist.gov/cgi-bin/cuu/Value?r
        #Gefittete Parameter in agingmodelbattery
def QLOSSCYC(P, Tb, ts, SoC, i, Vnenn=3.6, Cnom=2.2):
        
    CR=P/(Vnenn*Cnom)
    Ah=CR*ts*i/3600*2 # SOC dependency excluded because DoD doesnt contribute to Ah-throughput here
    Qlosscyc=((-47.84*CR**3)+(1215*CR**2)-(9419*CR)+(36040))*\
    np.exp((-31700+(370.3*CR))/(8.3144598*(Tb+273.15)))*(Ah**0.55) 
    #Note: This model works with kelvin temperature, the output is in percent capacity loss 
    return Qlosscyc    
    
        #Calendric Aging       
def QLOSSCAL(SoC, Tb,ts,i, Cnom=2.3, Cmax=40):
    A=4.39*(10**(-5))*np.exp(-(182/8.3144598)*((1/(Tb+273.15))-(1/(25+273.15))))
    B=1.01*(10**(-3))*np.exp(-(52.1/8.3144598)*((1/(Tb+273.15))-(1/(25+273.15))))
    # note: the model for the coeffiences A  and B is based on Kelvin Temp
    k=A*SoC+B
    Qlosscal=Cnom-Cmax
    alpha=np.exp(-(0.2*60-np.log(4))+(Tb*0.2))+3
    # note: the alpha fit curve requires Celcius Temp 
    Qlosscalit=k*(1+(Qlosscal/Cnom))**(-alpha)*ts/24/3600
    return Qlosscalit   
        
    #################################################################
 