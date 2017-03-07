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

#return the dispersed time (sec) accross two frequencies in MHz
def dispersion_t_sec(DM, nuLowMHz, nuHighMHz):
    return DM*4.148808*10**3*((1./nuLowMHz)**2-(1./nuHighMHz)**2)

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
    return S_delatDM_vs_S_ratio(k)


#return dDM
def dDMi(ssratio, W_ms, freq_centeral_GHz, Bandwidth_MHz):
    #dDM = 1000.
    for dDM in (0.01, 0.1, 1, 10, 100, 1000):
        #print dDM, snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz), '1', ssratio , '2'
        if snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz) > ssratio: continue
        else:
            while snrratio(dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz) < ssratio:
                dDM -= 0.001
                #if dDM<0.001:return dDM
            return dDM
    return 1001
    
def cal_snrratio(DM, dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec):
    w_ms = np.array([W_ms, 10.**3*dispersion_t_sec(DM, freq_centeral_GHz*10**3-Bandwidth_MHz/channels, freq_centeral_GHz*10**3+Bandwidth_MHz/channels), tInt_sec]).max()
    #print  '\n', 'un_disperse_able_width_ms', w_ms#, 10.**3*dispersion_t_sec(DM, freq_centeral_GHz*10**3-0.5*Bandwidth_MHz, freq_centeral_GHz*10**3-0.5*Bandwidth_MHz+Bandwidth_MHz/channels),'\n'
    k = kersci(dDM, Bandwidth_MHz, w_ms, freq_centeral_GHz)
    return S_delatDM_vs_S_ratio(k)

def cal_dDMi(ssratio, DM, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec):
    w_ms = np.array([W_ms, 10.**3*dispersion_t_sec(DM, freq_centeral_GHz*10**3-Bandwidth_MHz/channels, freq_centeral_GHz*10**3+Bandwidth_MHz/channels), tInt_sec]).max()
    #print  '\n', 'un_disperse_able_width_ms', w_ms#, 10.**3*dispersion_t_sec(DM, freq_centeral_GHz*10**3-0.5*Bandwidth_MHz, freq_centeral_GHz*10**3-0.
    return dDMi(ssratio, w_ms, freq_centeral_GHz, Bandwidth_MHz)


def main(args):

    ssratio         = 0.8
    DM              = 271#69.6937 #826.421
    dDM             = 33.043#869.705-862.421
    pulse_width_sec = 0.089328385024/512*40

    centeralfreq_MHz= 835.5957031           
    Bandwidth_MHz   = 31.25
    channels        = 320.
    tInt_sec        = 0.0003276796875
    pulse_width_sec = tInt_sec

    W_ms              = 1.* pulse_width_sec*1000
    freq_centeral_GHz = 1.* centeralfreq_MHz/1000

    print "pulse width       =", '%.2f' %W_ms, "ms"
    #print "freq_centeral_GHz =",'%.2f' % freq_centeral_GHz, 'GHz'
    #print "Bandwidth_MHz    =", '%.2f' % Bandwidth_MHz, 'MHz'
    print 'SSratio =', cal_snrratio(DM, dDM, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec)

    print 'dDM =', cal_dDMi(ssratio, DM, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec)


    """
    #ssratio = snrratio(DM, W_ms, freq_centeral_GHz, Bandwidth_MHz)

    #print dDMi(ssratio, W_ms, freq_centeral_GHz, Bandwidth_MHz), 'S(error)/S', ssratio

    #sys.exit()
    #print snrratio(dDM, 1, 1, 1)

    dDM = np.arange(-50,50,1.01)
    sso = []
    for i in dDM:
        #sso.append(snrratio(i, 1., 1., 1.))
        sso.append( cal_snrratio(DM, i, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec) )
    plt.plot(dDM, sso)
    """

    sso=0.5
    DMtrial = 1.
    DM =[]
    for i in range(300):
        dDM = cal_dDMi(sso, DMtrial, W_ms, freq_centeral_GHz, Bandwidth_MHz, channels, tInt_sec)
        DM.append(DMtrial +  dDM)
        DMtrial += dDM
    plt.plot(np.arange(i+1), DM)
    plt.xlabel('DM trial index')
    plt.ylabel(r'DM (pc cm$^{-3}$)')

    plt.show()
    #"""

if __name__ == '__main__':
    main(sys.argv[1:])