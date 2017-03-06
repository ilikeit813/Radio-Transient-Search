import numpy as np
import math
import sys


#Searches for radio transient equation (13) by Cordes&Mclaughflin
def kersci(deltaDM = 1., BandwidthMHz = 1., width_ms = 1., nuGHz = 1.):
    return 6.91*10**(-3)*deltaDM*BandwidthMHz/width_ms/nuGHz**3

#Searches for radio transient equation (12) by Cordes&Mclaughflin
def S_delatDM_vs_S_ratio(kersci):
    return np.sqrt(np.pi)/2.*kersci**(-1)*math.erf(kersci)

#return dispersion lag in second
def dispersion_lag_second(DM, nuLowMHz, nuHighMHz):
    return DM*4.148808*10**3*((1./nuLowMHz)**2-(1./nuHighMHz)**2)

def dDMi(DMtrial = 1., nuCenteralMHz = 2., channelMHz = 0.1, BMHz = 1., SSratio = 0.95, temporal_resol = 0.00021, pulsewidth = 0.01):
    dDM = 0.001
    w_sec = dispersion_lag_second(DM = DMtrial, nuLowMHz = nuCenteralMHz - 0*BMHz/2 - channelMHz/2, nuHighMHz =  nuCenteralMHz - 0*BMHz/2 + channelMHz/2)
    w_sec = np.array([pulsewidth, temporal_resol, w_sec]).max()
    print 'Un-de-dispersed-able width', w_sec, 'sec', '(temporal resol. =', temporal_resol, 'sec)'
    k = kersci(deltaDM = 1.*dDM, BandwidthMHz = 1.*BMHz, width_ms = 1.*w_sec*1000, nuGHz = 1.*nuCenteralMHz/1000)
    counter = 0
    if S_delatDM_vs_S_ratio(k) <= SSratio: return dDM
    while S_delatDM_vs_S_ratio(k) > SSratio: # > 0.01*SSratio:
        dDM += 0.001
        k = kersci(deltaDM = 1.*dDM, BandwidthMHz = 1.*BMHz, width_ms = 1.*w_sec*1000, nuGHz = 1.*nuCenteralMHz/1000)
        counter += 1
        if counter>10000:
            break
    return(dDM)

def main(args):
    ssratio         = 0.5
    pulse_width     = 0.089328385024/512*40
    DM              = 69.                   #pc/cc
    centeralfreq    = 835.5957031           #MHz
    bandwidth       = 31.25                 #MHz
    tInt            = 0.0003276798693344222 #sec
    numberOfchannel = 320.                  #bins
    print "pulse width = ", pulse_width, "sec"
    print 'channel width', 1.*bandwidth/numberOfchannel*1000, 'kHz'
    ddm = dDMi(1.*DM, 1.*centeralfreq, 1.*bandwidth/numberOfchannel, 1.*bandwidth, 1.*ssratio, 1.*tInt, 1.*pulse_width)
    print 'DM error = ',ddm, 'pc/cc, S/R ratio =', ssratio

if __name__ == '__main__':
    main(sys.argv[1:])