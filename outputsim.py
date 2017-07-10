# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 20:42:12 2017

@author: josch
"""

#################################
# output for simple simulations
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
print('---------------------------------------------------------')
print('rel. power losses')
print('---------------------------------------------------------')

#plt.plot(np.linspace(0,n*ts/3600,n+1), Plosspv/G, '-c', label='PlossPV')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaBMS, '--y', label='etalossBMS')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaMPPT, '-m', label='PlossMPPT')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaBAT, '--c', label='PlossBAT')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaBATch, '--b', label='etalossBATch')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaBATdch, '--m', label='etalossBATdch')
#plt.plot(np.linspace(0,n*ts/3600,n+1), etaINV, '--k', label='etalossINV')

plt.xlabel('Time [h]')
plt.ylabel('eta')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
plt.show()

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
'''
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