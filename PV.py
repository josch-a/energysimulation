# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 23:38:08 2017

@author: joscha
"""
import numpy as np

#################################################################
# PV #######################################################

        # The effective irradiance effective beam irradiance, 
        # effective sky diffuse irradiance,effective ground-reflected irradiance
def PVG(Gb,Gd,Gr):
    return Gb+Gd+Gr
    
        # The cell temperature in degrees Celsius (◦C):  
def PVTc(Ta, G, dT, a, b,WS):
        #The module back temperature in degrees Celsius (◦C), 
        #WS=Wind speed at 10 m height
    Tm=G*np.exp(a+b*WS)+Ta 
    return Tm + dT*(G/1000)   
        
    ## Output-Power Vergleichswerte der Genauigkeit sind in (Notton, Lazarov, & Stoyanov, 2010)
    
def PPV(Tc, G, gamma=-0.5, Gref=1000, Pmpref=100, Tref=25):
    return (G/Gref)*(Pmpref+(gamma*(Tc-Tref)))
    