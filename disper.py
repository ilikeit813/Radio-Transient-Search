"""
Calculate the DM error of Signal noise ratio 
"""

import numpy as np
import math
import sys
import math
import re
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Searches for radio transient equation (13) by Cordes&Mclaughflin
def kersci(deltaDM = 1., BandwidthMHz = 1., width_ms = 1., nuGHz = 1.):
    return 6.91*10**(-3)*deltaDM*BandwidthMHz/width_ms/nuGHz**3

#Searches for radio transient equation (12) by Cordes&Mclaughflin
def S_delatDM_vs_S_ratio(kersci):
    return np.sqrt(np.pi)/2.*kersci**(-1)*math.erf(kersci)

#return dispersion lag in second
def dispersion_lag_second(DM, nuLowMHz, nuHighMHz):
    return DM*4.148808*10**3*((1./nuLowMHz)**2-(1./nuHighMHz)**2)

#return S/S ratio
def snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz):
    k = kersci(dDM, Bandwidth_MHz, W_ms, freq_centeral_GHz)
    return np.sqrt(np.pi)/2./k*math.erf(k)

#return dDM
def dDMi(ssratio, W_ms, freq_centeral_GHz, Bandwidth_MHz):
    dDM = 1000.
    while snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz)<ssratio:
        dDM -= 0.001
        if dDM<0.001:return dDM
    return dDM
    

def main(args):
    dDM               = 30. 
    ssratio         = 0.5
    pulse_width_sec = 0.089328385024/512*40
    DM              = 690                   #pc/cc
    centeralfreq_MHz= 835.5957031           
    bandwidth_MHz   = 31.25

    W_ms              = 1.* pulse_width_sec*1000
    freq_centeral_GHz = 1.* centeralfreq_MHz/1000
    Bandwidth_MHz     = 1.* bandwidth_MHz
    print "pulse width       =", '%.2f' %W_ms, "ms"
    print "freq_centeral_GHz =",'%.2f' % freq_centeral_GHz, 'GHz'
    print "Bandwidth_MHz    =", '%.2f' % Bandwidth_MHz, 'MHz'
    print 'S(error)/S', snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz)

    ssratio         = snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz)

    print 'dDM =', dDMi(ssratio, W_ms, freq_centeral_GHz, Bandwidth_MHz), 'S(error)/S', ssratio

    sys.exit()
    print snrratio(dDM, 1, 1, 1)
    dDM = np.arange(-2300,2300,1.1)
    sso = []
    for i in dDM:
        sso.append(snrratio(i, 1., 1., 1.))

    plt.plot(dDM, sso)
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])