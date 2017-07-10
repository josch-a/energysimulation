# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 20:11:46 2017

@author: joscha
"""
import numpy as np

########################################################################
# import Demand curve from external file
def PLOAD (n,ts,tauy):
    import pandas
    # add the 0th timestep
    Pload=[0]
    
    # import data
    data = pandas.read_csv('loadcurvesim.csv', header=None, usecols=[1], sep=';', skiprows=0,nrows=((n)) )
    Pload=np.append(Pload, data.values[:,0])
    
    # adapt input data if the time steps are longer than 1 min                    
    if ts>60.:
        m=int(ts/60)
        jj=0
        lst=[0]
        while jj<=len(Pload)-m:
            ss=0
            ii=0
            while ii<=m-1:
                ss=Pload[ii+jj]+ss
                ii=ii+1
            lst=np.append(lst,ss/m)
            jj=jj+m
        Pload=lst
        
    # adapt the input data if simulation timespan is longer than 1 year
    if tauy>1.:
        kk=1
        Pz=Pload[1:]*1.
        while kk<tauy:
            Pz=np.append(Pz,Pload[1:])
            kk=kk+1
        Pz=np.append([0],Pz)   
        Pload=Pz
    # adapt the input data if simulation timespan is shorter than 1 year    
    if tauy<1.:
        del Pload[n-1:]
    return Pload


#########################################################################   
# import ambient temperature data from external file
def TA (n,ts,tauy):
    import pandas
    #add the 0th timestep
    TA=[25+273.15]
    # import data
    data = pandas.read_csv('BerlinSoDa_MERRA2.csv', header=None, usecols=[4], sep=';', skiprows=22,nrows=(n), decimal=','  )
    TA=np.append(TA, data.values[:,0])
    
    # adapt input data if the time steps are longer than 1 min 
    if ts>60.:
        m=int(ts/60)
        jj=0
        lst=[25+273.15]
        while jj<=len(TA)-m:
            ss=0
            ii=0
            while ii<=m-1:
                ss=TA[ii+jj]+ss
                ii=ii+1
            lst=np.append(lst,ss/m)
            jj=jj+m
        TA=lst
        
    # adapt the input data if simulation timespan is longer than 1 year   
    if tauy>1.:
        kk=1            
        Tz=TA[1:]*1.
        while kk<tauy:
            Tz=np.append(Tz,TA[1:])
            kk=kk+1
        Tz=np.append([25+273.15],Tz)
        TA=Tz
        
        # adapt the input data if simulation timespan is shorter than 1 year    
    if tauy<1.:
        del TA[n-1:]
    return TA-273.15


###########################################################################  
# import weather data from external file
def WSP (n,ts,tauy):
    import pandas
    # add the 0th timestep
    WSp=[0]    
    # import data
    data = pandas.read_csv('BerlinSoDa_MERRA2.csv', header=None, usecols=[4,7], sep=';', skiprows=22,nrows=(n), decimal=','  )
    WSp=np.append(WSp, data.values[:,1])
    # adapt input data if the time steps are longer than 1 min 
    if ts>60.:
        m=int(ts/60)
        jj=0
        lst=[0]
        while jj<=len(WSp)-m:
            ss=0
            ii=0
            while ii<=m-1:
                if WSp[ii+jj]!='float64':
                    WSp[ii+jj]=0
                ss=WSp[ii+jj]+ss
                ii=ii+1
            lst=np.append(lst,ss/m)
            jj=jj+m
        WSp=lst
    # adapt the input data if simulation timespan is longer than 1 year    
    if tauy>1.:
        kk=1
        Wz=WSp[1:]*1.
        while kk<tauy:
            Wz=np.append(Wz,WSp[1:])
            kk=kk+1
        Wz=np.append([0],Wz)
        WSp=Wz
    # adapt the input data if simulation timespan is shorter than 1 year    
    if tauy<1.:
        del WSp[n-1:]
    return WSp

#########################################################################    
# import irradiation data from external file
def GTI (n,ts,tauy):
    import pandas
    # add the 0th timestep
    G=[0] 
    # import data
    data = pandas.read_csv('BerlinSODA.csv', header=None, usecols=[4], sep=';', skiprows=22,nrows=(n), decimal=',' )
    G=np.append(G, data.values[:,0])
    # adapt input data if the time steps are longer than 1 min 
    if ts>60.:
        m=int(ts/60)
        jj=0
        lst=[0]
        while jj<=len(G)-m:
            ss=0
            ii=0
            while ii<=m-1:
                ss=G[ii+jj]+ss
                ii=ii+1
            lst=np.append(lst,ss/m)
            jj=jj+m
        G=lst
        
    # adapt the input data if simulation timespan is longer than 1 year
    if tauy>1.:
        kk=1
        Gz=G[1:]*1.
        while kk<tauy:
            Gz=np.append(Gz,G[1:])
            kk=kk+1
        Gz=np.append([0],Gz)
        G=Gz
        
    # adapt the input data if simulation timespan is shorter than 1 year    
    if tauy<1.:
        del G[n-1:]
    return G
    #################################################################                    





    
    
    
    
    





