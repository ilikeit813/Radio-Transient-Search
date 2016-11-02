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

def dDMi(DMtrial = 1., nuCenteralMHz = 2., channelMHz = 0.1, BMHz = 1., SSratio = 0.95, temporal_resol = 0.00021):
    dDM = 0.0001
    w_sec = dispersion_lag_second(DM = DMtrial, nuLowMHz = nuCenteralMHz - 0*BMHz/2 - channelMHz/2, nuHighMHz =  nuCenteralMHz - 0*BMHz/2 + channelMHz/2)
    if w_sec<temporal_resol: w_sec = temporal_resol
    k = kersci(deltaDM = 1.*dDM, BandwidthMHz = 1.*BMHz, width_ms = 1.*w_sec*1000, nuGHz = 1.*nuCenteralMHz/1000)
    #print 'Hi, is ',np.abs(S_delatDM_vs_S_ratio(k) - SSratio),"< 0.001 ?"
    while S_delatDM_vs_S_ratio(k) - SSratio > 0.000001:
        dDM += 0.001
        w_sec = dispersion_lag_second(DM = DMtrial, nuLowMHz = nuCenteralMHz - 0*BMHz/2 - channelMHz/2, nuHighMHz =  nuCenteralMHz - 0*BMHz/2 + channelMHz/2)
        if w_sec<temporal_resol: w_sec = temporal_resol
        k = kersci(deltaDM = 1.*dDM, BandwidthMHz = 1.*BMHz, width_ms = 1.*w_sec*1000, nuGHz = 1.*nuCenteralMHz/1000)
    return(dDM)
