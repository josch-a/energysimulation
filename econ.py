# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:45:29 2017
Economic Parameters
@author: josch
"""
import numpy as np

def CRF(i,n):
    return (i*(1+i)**n)/((1+i)**n-1)

def CELF(fCRF,i,n,r):
    k=(1+r)/(1+i)
    CRF=fCRF(i,n)
    return (k*(1-k**n))/(1-k)*CRF

def ARC(fCC,fCRF,trep, ny, iy, Enom,rfac):
    Nrep=len(trep)-1
    k=1
    RCtot=0
    #iff=(1+iy)**(1/(3600*24*365))-1
    for k in range(1,Nrep):
        RC=Enom*rfac*fCC(trep[k])/(1+iy)**(trep[k]/(3600*24*365))
        RCtot=RCtot+RC
        
    return fCRF(iy,ny)*RCtot

def ARV(fCRF,fCC,SoD,i,n, ny, iy, Enom):
    return ((1-SoD)*fCC(n)*Enom)/((1+iy)**ny)*fCRF(iy,ny)#*0.

def ACC(fCRF,fCC,iy,ny,Enom):
    return fCRF(iy,ny)*fCC(0)*Enom

def AOMC(fCELF,fCRF,OMC,i,n,rOM):
    return fCELF(fCRF,i,n,rOM)*OMC

def ATLCC(fARC,fCC,SoD,trep,fCRF,fCELF,i,n,rOMy, OMC, ny, iy, Nom,rfac):
    Aomc=fCELF(fCRF,iy,ny,rOMy)*OMC
    Acc=fCRF(iy,ny)*fCC(0)*Nom*rfac
    Arv=((1-SoD)*fCC(n)*Nom*rfac)/((1+iy)**ny)*fCRF(iy,ny)#*0.
    Arc=fARC(fCC,fCRF,trep, ny, iy, Nom,rfac)
    return Aomc+Acc-Arv+Arc

def LCoE(ATLCC,Etot):
    return ATLCC/Etot

# cost functions
def CCBAT(t):       # in [€/kWh]
    CC= 240+np.exp(-0.245620298*(t/365/24/3600+2017)+500)
    return CC
    
def CCMPPT(t):      # in [%/WpPV]
    CC=0.6*0.06
    return CC
    
def CCBMS(t):       # in [%/WhBattery]
    CC=0.275
    return CC
    
def CCINV(t):       # in [€/W]
    CC=0.225
    return CC
    
def CCPV(t):        # in [€/Wp]
    CC=1.1
    return CC
    