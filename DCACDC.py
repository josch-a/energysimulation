# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 23:24:21 2017

@author: josch
"""

    #################################################################
    #DC-DC/AC, mod 0 is dependent on Pout, else is dependent on Pin #
    
def PDC(P, vloss, rloss, Pself, mod, Pmax, Pnom):
    if P>Pmax/Pnom:
        Pi=Pmax/Pnom
        if mod==0:
            eta=Pi/(Pi+Pself+(vloss*Pi)+rloss*Pi**2)
            Pi=(Pi/eta)*Pnom
            LPS=(-Pmax+P*Pnom)*eta
                #print('mod0,Pnom>Pmax',Pi)
        elif mod==1:
            eta=-((1+vloss)/(2*rloss*Pi))+(((1+vloss)/(2*rloss*Pi))**2+\
                          ((Pi-Pself)/(rloss*Pi**2)))**0.5
            Pi=eta*Pi*Pnom  
            LPS=-Pmax+P*Pnom
    
                #print('mod1,Pnom>Pmax',Pi)
    elif P<=Pself:
        LPS=0.
        if mod==0:
            Pi=+Pself*Pnom+P*Pnom
            eta=1
        elif mod==1:
            Pi=-Pself*Pnom+P*Pnom
            eta=1
            #print('P=0',Pi)
    else:
        Pi=P
        LPS=0.        
        if mod==0:
            eta=Pi/(Pi+Pself+(vloss*Pi)+rloss*Pi**2)
            Pi=(Pi/eta)*Pnom
                #print('mod0,P=Pi',Pi)
        elif mod==1:
            eta=-((1+vloss)/(2*rloss*Pi))+(((1+vloss)/(2*rloss*Pi))**2+\
                          ((Pi-Pself)/(rloss*Pi**2)))**0.5
            Pi=eta*Pi*Pnom
                #print('mod1,P=Pi',Pi)
    return Pi, eta, LPS
    #################################################################