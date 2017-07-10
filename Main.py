
"""
Created on Wed Feb 22 13:56:59 2017

@author: joscha
"""


#######################################
#        Main Programm                #
#######################################

###################
#Import libraries 
###################
import time
starttime=time.clock()

#from scipy.optimize import *
import numpy as np
import matplotlib.pyplot as plt

import econ as eco
from input import *
from DCACDC import *
from battery import *
from PV import *
from water import *

#################
#initialization
#################
#length of time step [s]
ts=60*60
# [a]/ [s]  time horizon 
tauy=1     # in years
tau=tauy*365*24*3600    # in seconds without consideration of leap years NEEDS ADJUSTMENT
# number of time steps 
n=int(tau/ts) 

# initial randomization factor if no Monte-Carlo-Sim. is applied
rCCpv=1
rCCbat=1
rCCinv=1
rieffy=1
rrnomy=1

###############################################
#choose simulation mode: 
    #   optimization=0 ; 
    #   Monte-Carlo-Simulation: 1
    #   simple simulation: 2
mode=2
###############################################

# data for mode 1 & 2
Psol=1271.50368    # [W] - installed solar peak power
Cbatt=5685.13273   # [Wh] - installed battery capacity
Pwr=2456.61694     # [W] - installed inverter power
dv=(Psol,Cbatt,Pwr)
mc=5             # numbers of simulation rounds for Monte-Carlo-Simualation

###############################
#import data                    
##############################                    


     
# initial values weather and load
Pload=PLOAD(365*24*60,ts,tauy)
G=GTI(365*24*60,ts,tauy)
Ta=(TA(365*24*60,ts,tauy))
WS=WSP(365*24*60,ts,tauy)

    
          
#######################################################
def MAIN(dv):
    #start point
    i=1
    ####################################################################
    #Components (boundary conditions)
    
    ####################################################################
    #PV
    Gref= 1000. #[W/m²] reference irradiation (STC)
    gamma= -0.5 #[W/K] ±0.05 http://www.beamled.com/100w-semi-flexible-monocrystalline-solar-panel.html
    Pmpnom=dv[0] #[Wp] PV peak power - first decision variable 
    Trefpv=25. #[°C] reference temperature of the PV module
    PeolPV=0.8*Pmpnom   # End-of-Life condition of battery
    dgyPV=0.005         #yearly degradation of PV peak power
    dgPV=(1+dgyPV)**(1/(3600*24*365/ts))-1 # degradation of PV peak power per timestep
         
         ## PV Temperature model Gilman, P 2015 - SAM, Sandia Model (King,2004)
    a0=-3.56  #Glass/Cell Polymer Sheet -Open Rack 
    b0=-0.075
    dT0=3
    '''
    a1=-3.47  #Glass/Cell/Glass - Open Rack 
    b1=-0.0594
    dT1=3
    a2=-3.58  #Polymer/Thin Film/Steel - Open Rack 
    b2=-0.113
    dT2=3
    a3= -2.81 #Glass/Cell/Polymer Sheet Insulated Rack 
    b3=-0.0455
    dT3=0
    a4= -2.98 #Glass/Cell/Glass - Close Roof Mount 
    b4=-0.0471
    dT4=1
    '''
    ####################################################################
    #Battery
    hbat=2.       # [W/m²K] - heat transfer coefficient of the battery -
                 # http://www.schweizer-fn.de/waerme/waermeuebergang/waerme_uebergang.php
    cpbat=741    # [J/kgK] - specific heat capacity, battery [Nanda et al., 2014]
    mbatcell= 1.6   # +-0.05 [kg] - mass of a battery cell                                
    EnomBAT=dv[1] # [Wh] nominal Capacity at nominal C-Rate and nominal voltage - second decision variable
    mbat= mbatcell*EnomBAT/40/3.6 # [kg]  mass of the battery

    Vnombat=3.6 # [V] - nominal battery voltage
    Cnom=EnomBAT/Vnombat     # [Ah] - nominal battery capacity at nominal C-Rate
    CRnom=0.5    # [h^-1] - nominal C-Rate 
    PnomBAT=CRnom*EnomBAT # [W] Power at nominal C-Rate
    Abat=2*(0.115*0.047*EnomBAT/40/3.6 + 0.047*0.183*EnomBAT/40/3.6 + 0.115*0.183) # [m²] - area of a battery                               
    PselfBAT=0.25/(320*24*60*60)  # [%/s] Self-discharge of the battery    
    TopBATmin=-45. # [°C] - minimal operational temperature by manufactor's data sheet
    TopBATmax=85. # [°C] - maximal operational temperature by manufactor's data sheet    
    CeolBAT=0.4*Cnom   # End-of-Life condition of battery 
     
        # End-of-Discharge Parameters
    Vedch=2.8     # [V] - nominal end of discharge voltage
    Vech=4.0     # [V] - nominal end of charge voltage
    Prefedch=Vedch*40*CRnom  # [W] - reference power at end-of-discharge 
    Prefech=Vech*40*CRnom  # [W] - reference power at end-of-discharge
    a1batedch=-0.25359 # cut-off by Temperature, Parameter 1
    b1batedch=0.153657 # cut-off by Temperature, Parameter 2
    c1batedch=0.73176 # cut-off by Temperature, Parameter 3
    a2batedch=-0.021053 # cut-off by Discharging Power, Parameter 1
    b2batedch=0.000352111 # cut-off by Discharging Power, Parameter 2
    
        # End-of-charge Parameters 
    a1batech=0. # cut-off by Temperature, Parameter 1
    b1batech=1. # cut-off by Temperature, Parameter 2
    c1batech=0. # cut-off by Temperature, Parameter 3
    a2batech=1.14101 # cut-off by Charging Power, Parameter 1
    b2batech=-0.000225644 # cut-off by Charging Power, Parameter 2
    
    ###############################################################
    #DC-AC-inverter
    PnomINV=dv[2] # nominal inverter output power - 3rd decision variable
    PmaxINV=PnomINV*2 # maximum 
        # efficiency model parameter
    etanominv=0.98 # [-] nominal inverter efficiency
            # for output-power driven model
    vlossinv= 0.
    rlossinv=0.0375
    Pselfinv=0.0072
            # for input-power driven model
    vlossinvst= vlossinv
    rlossinvst = rlossinv*etanominv
    Pselfinvst = Pselfinv/etanominv
        # End-of-life conditions
    teolINVy=10.                          # [a] lifetime of inverter in years
    teolINV=teolINVy*365*24*3600         # [s] lifetime of inverter in sec - leap year ADOPTION IS REQUIRED
    
    ################################################################
    #MPPT-Charge Controller
    PnomMPPT=Pmpnom       # [W]  nominal charge-controller power
    PmaxMPPT=PnomMPPT*1   # [W]  maximum charge-controller power
        # efficiency model parameter
    etanomMPPT=0.975      # [-] nominal charge-controller efficiency
    vlossMPPT=0.00607463
    rlossMPPT=0.02281065
    PselfMPPT=0.00408267
    vlossMPPTst= vlossMPPT
    rlossMPPTst = rlossMPPT*etanomMPPT
    PselfMPPTst = PselfMPPT/etanomMPPT
        # End-of-life conditions
    teolMPPTy=5.                          # [a] lifetime of charge-controller
    teolMPPT=teolMPPTy*365*24*3600         # [s] lifetime of charge-controller - leap year ADOPTION IS REQUIRED
    
    ####################################################################
    #BMS
    PnomBMS=max(PnomMPPT,PmaxINV) # [W]  nominal BMS power
    PmaxBMS=PnomBMS*1             # [W]  maximum BMS power
        # efficiency model parameter
    etanomBMS=0.778 # [-] nominal BMS efficiency
    vlossBMS=-1.077e-16
    rlossBMS=0.2222
    PselfBMS=5e-07
    vlossBMSst= vlossBMS
    rlossBMSst = rlossBMS*etanomBMS
    PselfBMSst = PselfBMS/etanomBMS
        # End-of-life conditions
    teolBMSy=5.                          # [a]
    teolBMS=teolBMSy*365*24*3600         # [s] leap year ADOPTION IS REQUIRED
    
                                           
    ##############
    #Input values
    ############
        
    # electrical load
    
    # Water demand
    '''
    Vdemd=0.250 #[m³/d]
    deltap= 50000 #[Pa]
    Vdotn=Vdemd*6 # [m³/h]
    spump=0       # Signal if pump is included in the simulation
    '''
    # Economic Data
    ieffy=0.05*rieffy  # [-]
    
    rnomy=0.005*rrnomy # [-]
    ieff=(1+ieffy)**(1/(3600*24*365/ts))-1
    '''
    rnom=(1+rnomy)**(1/(3600*24*365/ts))-1
    
    ystart=2017  #  The  starting year of the simulation
    '''
    #########
    #models
    #########
    
    #################################################################
    # Economic Analysis
    # Cost functions

    OMCBAT=eco.CCBAT(0)*0.0
    OMCPV=eco.CCPV(0)*0.02
    OMCBMS=eco.CCBMS(0)*0.0
    OMCMPPT=eco.CCMPPT(0)*0.0
    OMCINV=eco.CCINV(0)*0.0
    CCREST=(0.6/0.4*Pmpnom*eco.CCPV(0)-eco.CCMPPT(0)*Pmpnom*eco.CCPV(0)-PnomINV*eco.CCINV(0)) # 10% Installation costs+ other costs of BoS
    
    
 
    #################################
    # constraints/ system description
    #################################
    
    #################
    # initial values
    ################

    LPS=np.zeros(n+1)     # Loss of Power Supply
    #Pexcess=np.zeros(n+1) # Loss by excess power which cant be used by the system
    P1=np.zeros(n+1)      # Output power PV-panels
    P2=np.zeros(n+1)      # Output power MPPT
    P3=np.zeros(n+1)      # Input power inverter
    #P3test=np.zeros(n+1)
    P4=np.zeros(n+1)      # Input power BMS
    #P5=np.zeros(n+1)      # Power loss due to over production by PV + full batteries
    P6=np.zeros(n+1)      # Terminal Power Battery 
    P7=Pload*1.           # Electrical load w/o water purification
    P8=np.zeros(n+1)      # Electrical load of water purification
    Pmpact=np.ones(n+1)*Pmpnom # actual peak power of PV-array at the corresponding time step
    SoC=np.zeros(n+1)     # State-of-Charge of the battery
    SoC[0]=0.8            # initial value of SoC
    Tb=np.zeros(n+1)      # Battery temperature
    Tb[0]=Ta[1]           # initial value of battery temperature
    '''
    Plosspv=np.zeros(n+1) # Energy loss due to inefficencies of the PV-panel
    PlossBMS=np.zeros(n+1) # Energy loss due to inefficencies of the BMS
    PlossMPPT=np.zeros(n+1) # Energy loss due to inefficencies of the MPPT-Charge controller
    Plossinv=np.zeros(n+1) # Energy loss due to inefficencies of the inverter
    '''
    PlossBAT=np.zeros(n+1) # Energy loss due to inefficencies of the Battery
    Pbat=np.zeros(n+1)    # Effective charge/discharge power of the battery after Losses
    Cmax=np.zeros(n+1) # residual nominal capacity at nominal C-Rate
    Cmax[0]=Cnom            # initial nominal capacity 
    Emaxbat=np.zeros(n+1)+EnomBAT 
    Emaxbat[0]=Cnom*Vnombat
    Eloadtot=np.sum(Pload)*ts/3600
    '''
    Qlo=np.zeros(n+1)
    Qlo0=np.zeros(n+1)
    Qlocy=np.zeros(n+1)
    Ah=np.zeros(n+1)
    B=np.zeros(n+1)
    CR=np.zeros(n+1)
    E=np.zeros(n+1)
    PexcessMPPT=np.zeros(n+1)
    LPSBMS=np.zeros(n+1)
    LPSINV=np.zeros(n+1)
    PlossBATdch=np.zeros(n+1)
    PlossBATch=np.zeros(n+1)
    '''
    PselfBATabs=np.zeros(n+1)
    SoDMPPT=np.ones(n+1)
    SoDINV=np.ones(n+1)
    SoDPV=np.ones(n+1)
    SoDBAT=np.ones(n+1)
    SoDBMS=np.ones(n+1)
    
    # control values
    Pdiff=np.zeros(n+1)     # Central Energy knot Power difference
    PTth=np.zeros(n+1)     # theoretical max Terminal power per time step
    #SOCeff=np.zeros(n+1)   # effective end of charge/discharge
    
    # efficiencies
    '''
    etaPmaxMPPT= PDC(PmaxMPPT/PnomMPPT, vlossMPPTst, rlossMPPTst, PselfMPPTst,\
         1, PmaxMPPT, PnomMPPT)[0]/PmaxMPPT         #efficiency at Pmax point
    '''
    etaPmaxBMS=PDC(PmaxBMS/PnomBMS, vloss=vlossBMS, \
                rloss=rlossBMS, Pself=PselfBMS, mod=1, Pmax=PmaxBMS, \
                Pnom=PnomBMS)[0]/PmaxBMS     	   #efficiency at Pmax point
    
    etaPmaxINV= PDC(PmaxINV/PnomINV, vloss=vlossinv, \
          rloss=rlossinv, Pself=Pselfinv, mod=1, Pmax=PmaxINV, Pnom=PnomINV)[0]/PmaxINV        #efficancy at Pmax point
    
    #etaBATdch=np.ones(n+1)
    #etaBATch=np.ones(n+1)
    #etaBAT=np.ones(n+1)
    #etaINV=np.ones(n+1)
    #etaMPPT=np.ones(n+1)
    #etaBMS=np.ones(n+1)
    '''
    # Volume flows
    Vdot=np.zeros(n+1)
    Vtank=np.zeros(n+1)
    Vtank[0]=Vdemd
       '''  
    # initializing the replacement time of every component
    trepPV=[0]
    trepMPPT=[0]
    trepBAT=[0]
    trepBMS=[0]
    trepINV=[0]     
    trepiBMS=0
    trepiMPPT=0
    trepiINV=0
    
    ###################
    # start simulation
    #################
    while i<=n:
        
        # PV Panel 
        if WS[i]!='float64':
            WS[i]=0
        P1[i]=PPV(PVTc(Ta[i], G[i], dT0, a0, b0, WS[i]), G[i], gamma, Gref, Pmpact[i-1],\
          Trefpv)
        
       # MPPT
        P2[i] =PDC(P1[i]/PnomMPPT, vlossMPPTst, rlossMPPTst, PselfMPPTst,\
         1, PmaxMPPT, PnomMPPT)[0]
        '''
        etaMPPT[i]=PDC(P1[i]/PnomMPPT, vlossMPPTst, rlossMPPTst, PselfMPPTst,\
         1, PmaxMPPT, PnomMPPT)[1]
        PexcessMPPT[i]=PDC(P1[i]/PnomMPPT, vlossMPPTst, rlossMPPTst, PselfMPPTst,\
         1, PmaxMPPT, PnomMPPT)[2]
       '''
    
        '''
        # water purification/ DSM #############################################
        if spump==1.:
            if (SoC[i-1]>0.97)&(PDC((P7[i]+PPUMP(Vdot[i],deltap))/PnomINV*etanominv, vloss=vlossinv, \
              rloss=rlossinv, Pself=Pselfinv, mod=0, Pmax=PmaxINV*etaPmaxINV, \
              Pnom=PnomINV)[0]-P2[i]<0):
                Vdot[i]=Vdotn
                P8[i]= PPUMP(Vdot[i],deltap)
                Vdot[i]=Vdotn
            elif Vtank[i-1]<Vdemd: 
                Vdot[i]=Vdotn              # one day ahead water supply is required
                P8[i]=PPUMP(Vdot[i],deltap)
                
            else:
                P8[i]=0
                Vdot[i]=0
                    
            Vtank[i]=Vdot[i]*(ts/3600)+Vtank[i-1]-Vdemd*(ts/(3600*24))
        elif spump==0.:
            P8[i]=0.
              '''  
        # DC-AC ##########################################################
        P3[i]=PDC((P7[i]+P8[i])/PnomINV*etanominv, vloss=vlossinv, \
          rloss=rlossinv, Pself=Pselfinv, mod=0, Pmax=PmaxINV*etaPmaxINV, \
          Pnom=PnomINV)[0]
        '''
        etaINV[i]=PDC((P7[i]+P8[i])/PnomINV*etanominv, vloss=vlossinv, \
          rloss=rlossinv, Pself=Pselfinv, mod=0, Pmax=PmaxINV*etaPmaxINV, \
          Pnom=PnomINV)[1] 
        LPSINV[i]= PDC((P7[i]+P8[i])/PnomINV*etanominv, vloss=vlossinv, \
          rloss=rlossinv, Pself=Pselfinv, mod=0, Pmax=PmaxINV*etaPmaxINV, \
          Pnom=PnomINV)[2]   
        P3test[i]=P3[i]*1.
        #print(P3test[i])
        #print(P7[i])
        #print(P8[i])
        '''
        # battery temperature ##############################################
        Tb[i]=TB(Ta[i], Tb[i-1], PlossBAT[i-1], hbat, Abat, mbat, cpbat,ts)
    
        # advanced energy balance central knot including BMS #################
        #   0= P2[i]-P3[i] - LPS[i]/Pexcess - P4[i] 
        #   SoCEch > SoC > SoCEdch
        
            # available Solar Power or requested Demand
        Pdiff[i]= P2[i]-P3[i] 
            #############
            #discharging
            # extra energy is required and higher than BMS losses
        if (-PDC(np.abs(Pdiff[i])/PnomBMS*etanomBMS, vloss=vlossBMS, \
                rloss=rlossBMS, Pself=PselfBMS, mod=0, Pmax=PmaxBMS*etaPmaxBMS, \
                Pnom=PnomBMS)[0]<=-PselfBMS*PnomBMS)&(Pdiff[i]<0.):        
                # BMS
            PTth[i]=-PDC(np.abs(Pdiff[i])/PnomBMS*etanomBMS, vloss=vlossBMS, \
                rloss=rlossBMS, Pself=PselfBMS, mod=0, Pmax=PmaxBMS*etaPmaxBMS, \
                Pnom=PnomBMS)[0]
            '''
            etaBMS[i]=PDC(np.abs(Pdiff[i])/PnomBMS*etanomBMS, vloss=vlossBMS, \
                rloss=rlossBMS, Pself=PselfBMS, mod=0, Pmax=PmaxBMS*etaPmaxBMS, \
                Pnom=PnomBMS)[1]
            LPSBMS[i]=PDC(np.abs(Pdiff[i])/PnomBMS*etanomBMS, vloss=vlossBMS, \
                rloss=rlossBMS, Pself=PselfBMS, mod=0, Pmax=PmaxBMS*etaPmaxBMS, \
                Pnom=PnomBMS)[2]
            '''
            # battery is not empty and battery temperature lies in between the operation bounds
            if (SOC(PTth[i], ts, SoC[i-1], Emaxbat[i-1], 0)>soceff(Tb[i],-PTth[i]/EnomBAT*40*3.6,a1=a1batedch,b1=b1batedch,\
                c1=c1batedch,a2=a2batedch,b2=b2batedch, Pref=Prefedch))&\
                (TopBATmin<=Tb[i]<=TopBATmax) :
                
                #SOCeff[i]=soceff(Tb[i],-PTth[i],a1=a1batedch,b1=b1batedch,\
                 #     c1=c1batedch,a2=a2batedch,b2=b2batedch,Pref=Prefedch) 
                
                P6[i]=PTth[i] # terminal power at Battery
                P4[i]=Pdiff[i] # terminal power at BMS
                #PlossBMS[i]=np.abs(P4[i]-P6[i]) # Energy losses at BMS
                if P7[i]>PmaxINV*etaPmaxINV:
                    P7[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                          rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                          Pnom=PnomINV)[0]-P8[i]
            # battery is empty or temperature is out of range
            else:
                #SOCeff[i]=soceff(Tb[i],-PTth[i],a1=a1batedch,b1=b1batedch,\
    #c1=c1batedch,a2=a2batedch,b2=b2batedch, Pref=Prefedch) 
                
                
                P6[i]=0
                #etaBMS[i]=1. 
                      
                # solar power is sufficient to power BMS
                if P2[i]>=PselfBMS*PnomBMS:
                    P4[i]=-PselfBMS*PnomBMS
                    P3[i]=P2[i]+P4[i]
                    #PlossBMS[i]=PselfBMS*PnomBMS
                    P7[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                      rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                      Pnom=PnomINV)[0]  
                    '''
                    etaINV[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                      rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                      Pnom=PnomINV)[1] '''
                # solar power is NOT sufficiant to power BMS
                else:
                    P3[i]=0.
                    P4[i]=0
                    P2[i]=0.
                    P1[i]=0.   # is that possible doesnt the MPPT always need supply
                    #PlossBMS[i]=1.  # same here isnt the BMS const. using power
                    P7[i]=0.
                    
              
                P8[i]=0. # if battery is empty, water is not purified
                  
            ##########    
            #charging
            # excess power available & bigger than BMS-self-consumption
        elif (PDC(np.abs(Pdiff[i])/PnomBMS,vloss=vlossBMSst,rloss=rlossBMSst,\
    Pself=PselfBMSst, mod=1, Pmax=PmaxBMS, Pnom=PnomBMS)[0]>=PselfBMS*PnomBMS)&\
    (Pdiff[i]>0.):        
                # BMS
            PTth[i]=PDC(np.abs(Pdiff[i])/PnomBMS,vloss=vlossBMSst,rloss=rlossBMSst,\
                Pself=PselfBMSst, mod=1, Pmax=PmaxBMS, Pnom=PnomBMS)[0]
            '''
            etaBMS[i]=PDC(np.abs(Pdiff[i])/PnomBMS,vloss=vlossBMSst,rloss=rlossBMSst,\
                  Pself=PselfBMSst, mod=1, Pmax=PmaxBMS, Pnom=PnomBMS)[1]
            LPSBMS[i]=PDC(np.abs(Pdiff[i])/PnomBMS,vloss=vlossBMSst,rloss=rlossBMSst,\
                  Pself=PselfBMSst, mod=1, Pmax=PmaxBMS, Pnom=PnomBMS)[2]
            '''
            # battery is NOT full and battery temperature is in the range of the manual data 
            if (SOC(PTth[i], ts, SoC[i-1], Emaxbat[i-1], 0)<soceff(Tb[i],PTth[i]/EnomBAT*40*3.6,a1=a1batech,b1=b1batech,c1=c1batech,\
                a2=a2batech,b2=b2batech,Pref=Prefech))&(TopBATmin<=Tb[i]<=TopBATmax):
                
                #SOCeff[i]=soceff(Tb[i],PTth[i],a1=a1batech,b1=b1batech,\
                #      c1=c1batech,a2=a2batech,b2=b2batech,Pref=Prefech) 
                
                P6[i]=PTth[i]
                P4[i]=Pdiff[i]
                #PlossBMS[i]=Pdiff[i]-PTth[i]
                P7[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                      rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                      Pnom=PnomINV)[0]-P8[i]
            
            # battery is full or temperature is out of range
            else:
                #SOCeff[i]=soceff(Tb[i],PTth[i],a1=a1batech,b1=b1batech,\
   # c1=c1batech,a2=a2batech,b2=b2batech,\
            #Pref=Prefech) 
                
                # Adoption of Power flows in the system due to full battery
                P6[i]=0
                P4[i]=-PselfBMS*PnomBMS
                #Pexcess[i]=P2[i]+P4[i]-P3[i]
                P2[i]=P3[i]+P4[i]            
                P1[i]=PDC(P2[i]/PnomMPPT*etanomMPPT, vlossMPPT, rlossMPPT, PselfMPPT, \
                  0, PmaxMPPT,PnomMPPT)[0]
                 
                #etaBMS[i]=1.
                #PlossBMS[i]=PselfBMS*PnomBMS
                P7[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                      rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                      Pnom=PnomINV)[0]-P8[i]
            #idle    
        else:
            
            P4[i]=0.
            #PlossBMS[i]=PselfBMS*PnomBMS 
            #etaBMS[i]=1.
            P6[i]=-PselfBMS*PnomBMS #####NEED ADJ Battery empty and 
            P7[i]=PDC((P3[i])/PnomINV, vloss=vlossinvst, \
                      rloss=rlossinvst, Pself=Pselfinvst, mod=1, Pmax=PmaxINV, \
                      Pnom=PnomINV)[0]-P8[i]
            
        # Battery  ##########################################################
        # effective Charging/ Discharging Power
        Pbat[i]= PBAT(P6[i], Tb[i], PnomBAT)
        #SoC Charge
        SoC[i]=SOC(Pbat[i], ts, SoC[i-1], Emaxbat[i-1], 0) 
        
    
        # actual capacity of the battery [Ah]
        Cmax[i]=Cnom - QLOSS(abs(P6[i-1]), abs(P6[i]), Tb[i-1], Tb[i], ts, SoC[i-1],\
     SoC[i], i, Vnombat, Cnom, \
    Cmax[i-1]) # note: the cyc aging function cant handle negative Power
    
    
        # actual energy capacity of the battery [Wh]
        Emaxbat[i]=Vnombat*Cmax[i]
        # actual PV peak power [Wp]
        Pmpact[i]=(1-dgPV)*Pmpact[i-1]
        # Power losses  
        # charging
        if P6[i]>0.:
            #PlossBATch[i]=abs(P6[i]-Pbat[i])
            PlossBAT[i]=abs(P6[i]-Pbat[i])
            #etaBATch[i]=abs(Pbat[i]/P6[i])
        #idle
        #elif P6[i]==0.:
            #etaBAT[i]=1
        # discharging
        else: 
            #PlossBATdch[i]=abs(P6[i]-Pbat[i])
            PlossBAT[i]=abs(P6[i]-Pbat[i])
            #etaBATdch[i]=abs(P6[i]/Pbat[i])
            #etaBAT[i]=PlossBAT[i]
            
        PselfBATabs[i]=PselfBAT*ts*Emaxbat[i] 
        
        PlossBAT[i]=abs(P6[i]-Pbat[i])+PselfBAT*ts*Emaxbat[i]
         
        #etaBAT[i]=PlossBAT[i]
        #Plosspv[i]= G[i]-P1[i]
        #PlossMPPT[i]=P1[i]-P2[i]
        #Plossinv[i]=P3[i]-P7[i]-P8[i] 
        
        
        # Aging conditions of the components
        SoDBMS[i]=(i*ts-trepiBMS)/teolBMS
        SoDMPPT[i]=(i*ts-trepiMPPT)/teolMPPT
        SoDPV[i]=(Pmpnom-Pmpact[i])/(Pmpnom-PeolPV)
        SoDBAT[i]=(Cnom-Cmax[i])/(Cnom-CeolBAT)
        SoDINV[i]=(i*ts-trepiINV)/teolINV     
              
        # find the time of replacement
        if Cmax[i]<CeolBAT:
            trepiBAT=i*ts
            trepBAT=np.append(trepBAT,trepiBAT) # the vector with all replacement times
            Cmax[i]=Cnom
            
        if Pmpact[i]<PeolPV:
            trepiPV=i*ts
            trepPV=np.append(trepPV,trepiPV) # the vector with all replacement times
            Pmpact[i]=Pmpnom
        
        if teolMPPT<=i*ts:
            trepiMPPT=i*ts
            trepMPPT=np.append(trepMPPT,trepiMPPT) # the vector with all replacement times
        
        if teolBMS<=i*ts:
            trepiBMS=i*ts
            trepBMS=np.append(trepBMS,trepiBMS) # the vector with all replacement times
            
        if teolINV<=i*ts:
            trepiINV=i*ts
            trepINV=np.append(trepINV,trepiINV) # the vector with all replacement times
        
        #print('LPSINV',LPSINV[i])
        #if Pload[i]<P7[i]:
         #   print(i, Pload[i]-P7[i], Pexcess[i], P6[i])
        # next step ####################################################
        i=i+1
    
    LPS=1.*Pload-P7
    '''
    ARCBAT=eco.ARC(CCBAT,eco.CRF,trepBAT,tauy,ieffy,Emaxbat[0]/1000)
    ARVBAT=eco.ARV(eco.CRF,CCBAT,SoDBAT[n],ieff,tau,tauy,ieffy,Emaxbat[0]/1000)
    AOMCBAT=eco.AOMC(eco.CELF,eco.CRF,OMCBAT,ieffy,tauy,rnomy)
    ACCBAT=eco.ACC(eco.CRF,CCBAT,ieffy,tauy,Emaxbat[0]/1000)
    '''
    ATLCCBAT=eco.ATLCC(eco.ARC,eco.CCBAT,SoDBAT[n],trepBAT,eco.CRF,eco.CELF,ieff,tau,rnomy,OMCBAT,tauy,ieffy,Emaxbat[0]/1000,rCCbat)
    ATLCCPV=eco.ATLCC(eco.ARC,eco.CCPV,SoDPV[n],trepPV,eco.CRF,eco.CELF,ieff,tau,rnomy,OMCPV,tauy,ieffy,Pmpnom,rCCpv)
    ATLCCMPPT=eco.ATLCC(eco.ARC,eco.CCMPPT,SoDMPPT[n],trepMPPT,eco.CRF,eco.CELF,ieff,tau,rnomy,OMCMPPT,tauy,ieffy,Pmpnom*eco.CCPV(0),rCCpv)
    ATLCCBMS=eco.ATLCC(eco.ARC,eco.CCBMS,SoDBMS[n],trepBMS,eco.CRF,eco.CELF,ieff,tau,rnomy,OMCBMS,tauy,ieffy,Pmpnom*eco.CCPV(0),rCCpv)
    ATLCCINV=eco.ATLCC(eco.ARC,eco.CCINV,SoDINV[n],trepINV,eco.CRF,eco.CELF,ieff,tau,rnomy,OMCINV,tauy,ieffy,PnomINV,rCCinv)
    
    # Ignoring the installation and other costs for short-term simulations
    if tauy>=10.:
        ATLCCREST=CCREST*eco.CRF(ieffy,tauy)
    else:
        ATLCCREST=0.
    ATLCCtot=ATLCCBAT+ATLCCPV+ATLCCMPPT+ATLCCBMS+ATLCCINV+ATLCCREST
    
    
    LLP=np.sum(abs(LPS))/(Eloadtot)*ts/3600 # Loss-of-Load-Probability
    LCoE=ATLCCtot/((1.-LLP)*Eloadtot/1000/tauy)
    
    ################################################################
    # output for simple simulations
    ##########################################################
    
    if mode==2:
        #################################
       
        plt.plot(np.linspace(0,n*ts/3600,n+1), G, '-y', label='G')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P1, '-r', label='P1')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P2, '--m', label='P2')
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), P7+P8, '-b', label='P7+P8')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P3, '-g', label='P3')
        #plt.plot(np.linspace(0,n/ts,n+1), P3test, '-k', label='P3test')
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        
        plt.show()
        print('---------------------------------------------------------')
        print('real Energy flows around battery')
        print('---------------------------------------------------------')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P4, '-b', label='P4')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P6, '--g', label='P6')
        plt.plot(np.linspace(0,n*ts/3600,n+1), Pbat, '-k', label='Pbat')
        plt.plot(np.linspace(0,n*ts/3600,n+1), PlossBAT, '--c', label='PlossBAT')
        plt.xlabel('Time [m]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        
        plt.show()
        
        print('---------------------------------------------------------')
        print('real Energy flows around inv')
        print('---------------------------------------------------------')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P3, ':b', label='P3')
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), P7, ':k', label='P7')
        plt.plot(np.linspace(0,n*ts/3600,n+1), P8, '-g', label='P8')
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        
        plt.show()
        
        print('---------------------------------------------------------')
        print('theor. Energy flows around battery')
        print('---------------------------------------------------------')
        plt.plot(np.linspace(0,n*ts/3600,n+1), Pload, '-r', label='Pload')
        plt.plot(np.linspace(0,n*ts/3600,n+1), PTth, '--b', label='PTth')
        plt.plot(np.linspace(0,n*ts/3600,n+1), Pdiff, '--y', label='Pdiff')
        
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        
        plt.show()
        
        print('---------------------------------------------------------')
        print('excess and unmet power')
        print('---------------------------------------------------------')
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), LPS, '-r', label='LPS')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), Pexcess, '-m', label='Pexcess')
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        
        plt.show()
        print('---------------------------------------------------------')
        print('abs. power losses')
        print('---------------------------------------------------------')
        
        #plt.plot(np.linspace(0,n*ts/3600,n+1), Plosspv, '-c', label='PlossPV')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), PlossBMS, '--y', label='PlossBMS')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), PlossMPPT/abs(P2), '-m', label='PlossMPPT')
        plt.plot(np.linspace(0,n*ts/3600,n+1), PlossBAT, '--c', label='PlossBAT')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), PlossBATch, ':b', label='PlossBATch')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), PlossBATdch, ':m', label='PlossBATdch')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), Plossinv, '--k', label='PlossINV')
        plt.plot(np.linspace(0,n*ts/3600,n+1), PselfBATabs, '--g', label='PselfBAT')
        plt.xlabel('Time [h]')
        plt.ylabel('P')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        '''
        print('---------------------------------------------------------')
        print('rel. power losses')
        print('---------------------------------------------------------')
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), Plosspv/G, '-c', label='PlossPV')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaBMS, '--y', label='etalossBMS')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaMPPT, '-m', label='PlossMPPT')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaBAT, '--c', label='PlossBAT')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaBATch, '--b', label='etalossBATch')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaBATdch, '--m', label='etalossBATdch')
        plt.plot(np.linspace(0,n*ts/3600,n+1), etaINV, '--k', label='etalossINV')
        
        plt.xlabel('Time [h]')
        plt.ylabel('eta')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        '''
        print('---------------------------------------------------------')
        print('LPS ')
        print('---------------------------------------------------------')
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), LPS, '-k', label='LPStot')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), LPSBMS, ':y', label='LPSBMS')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), LPSINV, ':m', label='LPSINV')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), LPS-LPSINV-LPSBMS, ':g', label='LPSBAT')
        plt.xlabel('Time [h]')
        plt.ylabel('P [W]')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        '''
        print('---------------------------------------------------------')
        print('excess power ')
        print('---------------------------------------------------------')
        
        #plt.plot(np.linspace(0,n*ts/3600,n+1), Pexcess, '-c', label='Pexcess,BAT')
        #plt.plot(np.linspace(0,n*ts/3600,n+1), PexcessMPPT, '--y', label='PexcessMPPT')
        
        
        plt.xlabel('Time [h]')
        plt.ylabel('P [W]')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        
        #plt.plot(np.linspace(0,n*ts/3600,n+1), SOCeff, '.r', label='SoCeff')
        plt.plot(np.linspace(0,n*ts/3600,n+1), SoC, '.k', label='SoC', alpha=0.5)
        plt.xlabel('Time [h]')
        plt.ylabel('SoC')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        #fig.tight_layout()
        
        #pyplot.savefig("outputname.eps", dpi = 100)
        plt.show()
        
        plt.plot(np.linspace(0,n*ts/3600,n+1), Emaxbat, '-k', label='Emax,bat')
        plt.xlabel('Time [h]')
        plt.ylabel('Emax')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        plt.plot(np.linspace(0,n*ts/3600,n+1), Ta, '.r', label='Ta', alpha=0.5)
        plt.plot(np.linspace(0,n*ts/3600,n+1), Tb, '.g', label='Tb', alpha=0.5)
        
        plt.xlabel('Time [h]')
        plt.ylabel('T[°C]')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        
        
        print(LLP)
        print(np.sum(LPS))
        print(np.sum((LPS**2)**0.5))
        
        print('---------------------------------------------------------')
        print('Water volume in tank')
        print('---------------------------------------------------------')
        
        plt.plot(np.linspace(0,n/ts,n+1), Vtank, '-r', label='Vtank')
        plt.xlabel('Time [m]')
        plt.ylabel('Vtank [m³]')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
        plt.show()
        '''
        
        print('---------------------------------------------------------')
        print('Text output')
        print('---------------------------------------------------------')
        
        print('LCoE',LCoE)
        print('LLP',LLP)
        np.savetxt('1jmBER2.csv', ([LCoE,LLP,ATLCCBAT,ATLCCPV,ATLCCMPPT,\
              ATLCCBMS,ATLCCINV, ATLCCREST,ATLCCtot]), delimiter = ';', fmt='%10.5f')
        '''
        np.savetxt('1jmPBER2.csv', (P1,P2,P3,P4,P5,P6,P7,P8,SoC,SOCeff,Pexcess,\
                                 PexcessMPPT,Plosspv,PlossBMS,PlossBAT,\
                                 PlossBATch,PlossBATdch,Plossinv,PselfBATabs,\
                                 Pdiff,LPS, LPS-LPSINV-LPSBMS,LPSBMS,LPSINV,PTth,Pbat,\
                                 Pload,G, P7+P8,Ta,Tb,Emaxbat,etaBMS,etaMPPT,etaINV,\
                                 etaBATch,etaBATdch),delimiter = ';', fmt='%10.5f')
            
        '''
    #####################################################################
    
    print(time.clock()-starttime) 
    if mode==0:
        return [LCoE, LLP], [LCoE, LCoE-3,LLP, LLP-0.4] 
    else:
        return [LCoE, LLP]

if mode==0:
    ################################################################    
    # Optimization algorithm
    from platypus import NSGAII, Problem, Real
    
    problem = Problem(3, 2, 4)  # initialize Problem with two decision variables, two objectives, and two constraints 
    problem.types[:] = [Real(400, 30000), Real(400, 30000), Real(400, 30000)]
    problem.constraints[0] = ">=0"
    problem.constraints[1] = "<=0"
    problem.constraints[2] = ">=0"
    problem.constraints[3] = "<=0"
    problem.function = MAIN
    problem.directions[0] = Problem.MINIMIZE
    problem.directions[1] = Problem.MINIMIZE
                      
    # instantiate the optimization algorithm
    algorithm = NSGAII(problem)
    
    # optimize the problem using 10,000 function evaluations
    algorithm.run(1300)
    
    ########################################################################
    # optimization output
    #######################################################################
    
    # save results in .csv
    np.savetxt('optimizationresBERoneh.csv', ([s.objectives[0] for s in algorithm.result],\
                                       [s.objectives[1] for s in algorithm.result],\
                                       [s.variables[0] for s in algorithm.result],\
                                       [s.variables[1] for s in algorithm.result],\
                                       [s.variables[2] for s in algorithm.result]),\
     delimiter = ';', fmt='%10.5f')
    
    
    
    # plot the results 
    plt.scatter([s.objectives[0] for s in algorithm.result],
                [s.objectives[1] for s in algorithm.result],s=75, \
                c=[s.variables[0] for s in algorithm.result], alpha=.5)
    
    #plt.xlim([0, 1.1])
    #plt.ylim([0, 1.1])
    plt.xlabel("$LCoE(Pmpnom,Cnom)$")
    plt.ylabel("$LLP(Pmpnom,Cnom)$")
    cb=plt.colorbar()
    cb.set_label('installierte Solarleistung [$W_p$]')
    plt.show()
    
    '''
    plt.scatter([s.objectives[0] for s in algorithm.result],
                [s.objectives[1] for s in algorithm.result],s=75, \
                c=[s.variables[1] for s in algorithm.result], alpha=.5)
    
    #plt.xlim([0, 1.1])
    #plt.ylim([0, 1.1])
    plt.xlabel("$LCoE(Pmpnom,Cnom)$")
    plt.ylabel("$LLP(Pmpnom,Cnom)$")
    cb=plt.colorbar()
    cb.set_label('installierte Batteriekapazität [$Wh$]')
    plt.show()
    
    
    plt.scatter([s.objectives[0] for s in algorithm.result],
                [s.objectives[1] for s in algorithm.result],s=75, \
                c=[s.variables[2] for s in algorithm.result], alpha=.5)
    
    #plt.xlim([0, 1.1])
    #plt.ylim([0, 1.1])
    plt.xlabel("$LCoE(Pmpnom,Cnom)$")
    plt.ylabel("$LLP(Pmpnom,Cnom)$")
    cb=plt.colorbar()
    cb.set_label('installierte Wechselrichterleistung [$W_nom$]')
    plt.show()
    
    plt.scatter([s.variables[0] for s in algorithm.result],
                [s.variables[1] for s in algorithm.result])
    
    #plt.xlim([0, 1.1])
    #plt.ylim([0, 1.1])
    plt.xlabel("$Pmpnom$")
    plt.ylabel("$Cnom$")
    plt.show()
    '''
    
    '''
    for solution in algorithm.result:    
        print(solution)
    '''


elif mode==1:
    #######################################################################
    # monte-carlo-cycle
    i=0
    RES=np.zeros(shape=(mc,2))
    for i in range(0,mc):    
        # Randomisation factors
        rload=abs((1+1.51*np.random.randn(len(Pload))))
        rsun=(1+0.12*np.random.randn(len(Pload)))
        rT=(1+0.20*np.random.randn(len(Pload)))         
        
        #random numbers - normal distr. 
        rCCpv=(1+0.51*np.random.randn())
        rCCbat=(1+0.51*np.random.randn())
        rCCinv=(1+0.51*np.random.randn())
        rieffy=(1+0.51*np.random.randn())
        rrnomy=(1+0.51*np.random.randn()) 
        
        Pload=Pload*rload*(np.sum(Pload)/np.sum(Pload*rload))
        Ta=Ta*rT*(np.sum(Ta)/np.sum(Ta*rT))
        if np.sum(G*rsun)==0.:
            G=G
        else:    
            G=G*rsun*(np.sum(G)/np.sum(G*rsun))
            
        # simulation
        RES[i]=MAIN(dv)
        print(-starttime+time.clock())
        #print('RES',RES[i])    
    
    #######################################################################
    # Monte-Carlo-output
    ########################################################################
    
    np.savetxt('montecres27test.csv', RES, delimiter = ';', fmt='%10.5f')
    RES=np.array(RES)
    plt.plot(RES[:,0], RES[:,1], '.r')

    plt.xlabel(r"LCoE [$Euro/kWh_{el}$]")
    plt.ylabel("LLP [-]")
    plt.legend(loc='best')
    plt.show()
    
    plt.hist(RES[:,0], bins='auto')  # plt.hist passes it's arguments to np.histogram
    plt.title("Histogram with 'auto' bins")
    plt.show()
    plt.hist(RES[:,1], bins='auto')  # plt.hist passes it's arguments to np.histogram
    plt.title("Histogram with 'auto' bins")
    plt.show()

elif mode==2:
    ####################################################################
    #simple simulation
    RES=MAIN(dv)
    
print('total simulation time')
print(time.clock()-starttime)