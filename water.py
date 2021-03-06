# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 23:43:39 2017

@author: joscha
"""

    # Water purification
    
def PPUMP(Vdot,p): # p in [Pa] and Vdot in [m³/h]
    ap=0.2201*(p/9806.65)**3-9.928*(p/9806.65)**2+108.4*(p/9806.65)-214.4
    bp=-0.339*(p/9806.65)**3+15.04*(p/9806.65)**2-157.2*(p/9806.65)+470.5
    cp=0.1499*(p/9806.65)**3-6.547*(p/9806.65)**2+72.37*(p/9806.65)-152.8
    dp=-0.01584*(p/9806.65)**3+0.7072*(p/9806.65)**2-2.814*(p/9806.65)+16.79
    P=ap*Vdot**3+bp*Vdot**2+cp*Vdot+dp
    return P
    
    