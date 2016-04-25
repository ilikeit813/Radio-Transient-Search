#script of Python script to dedisperse a file, save each DM trials time series
#Codes from Sean, modified by Jamie
#V3, super cluster technique

from mpi4py import MPI
import sys
import numpy as np
import glob
import os
import time
import sys

def DMs(DMstart,DMend,dDM):
    """
    Calculate the number of DMs searched between DMstart and DMend, with spacing dDM * DM.

    Required:

    DMstart   - Starting Dispersion Measure in pc cm-3
    DMend     - Ending Dispersion Measure in pc cm-3
    dDM       - DM spacing in pc cm-3
    """

    #NDMs = np.log10(float(DMend)/float(DMstart))/np.log10(1.0+dDM)
    NDMs = (DMend-DMstart)/dDM

    return int(np.round(NDMs))


def delay2(freq, dm):
    """
    Calculate the relative delay due to dispersion over a given frequency
    range in Hz for a particular dispersion measure in pc cm^-3.  Return
    the dispersive delay in seconds.  Same as delay, but w.r.t to highest frequency.
    ***Used to simulate a dispersed pulse.***

    Required:

    freq - 1-D array of frequencies in MHz
    dm   - Dispersion Measure in pc cm-3
    """
    # Dispersion constant in MHz^2 s / pc cm^-3
    _D = 4.148808e3
    # Delay in s
    tDelay = dm*_D*((1/freq)**2 - (1/freq.max())**2)

    return tDelay


def Threshold(ts, thresh, clip=3, niter=1):
    """
    Wrapper to scipy threshold a given time series using Scipy's threshold function (in 
    scipy.stats.stats).  First it calculates the mean and rms of the given time series.  It then 
    makes the time series in terms of SNR.  If a given SNR value is less than the threshold, it is 
    set to "-1".  Returns a SNR array with values less than thresh = -1, all other values = SNR.
    Also returns the mean and rms of the timeseries.

    Required:

    ts   -  input time series.

    Options:

    thresh  -  Time series signal-to-noise ratio threshold.  default = 5.
    clip    -  Clipping SNR threshold for values to leave out of a mean/rms calculation.  default = 3.
    niter   -  Number of iterations in mean/rms calculation.  default = 1.

    Usage:
    >>sn, mean, rms = Threshold(ts, *options*)

    """
    #  Calculate, robustly, the mean and rms of the time series.  Any values greater than 3sigma are left
    #  out of the calculation.  This keeps the mean and rms free from sturation due to large deviations.

    mean = np.mean(ts) 
    std  = np.std(ts)  
    #print mean,std

    if niter > 0:
        for i in range(niter):
            ones = np.where((ts-mean)/std < clip)[0]  # only keep the values less than 3sigma
            mean = np.mean(ts[ones])
            std  = np.std(ts[ones])
    SNR = (ts-mean)/std
    SNR[SNR<thresh]=-1
    #SNR = st.threshold((ts-mean)/std, threshmin=thresh, newval=-1)

    return SNR, mean, std

def Decimate_ts(ts, ndown=2):
    """
    Takes a 1-D timeseries and decimates it by a factore of ndown, default = 2.  
    Code adapted from analysis.binarray module: 
      http://www.astro.ucla.edu/~ianc/python/_modules/analysis.html#binarray 
    from Ian's Python Code (http://www.astro.ucla.edu/~ianc/python/index.html)
    
    Optimized for time series' with length = multiple of 2.  Will handle others, though.

    Required:
    
    ts  -  input time series

    Options:
    
    ndown  -  Factor by which to decimate time series. Default = 2.
              if ndown = 1, returns ts       
    """

    if ndown==1:
       return ts

    ncols = len(ts)
    n_rep = ncols / ndown
    ts_ds = np.array([ts[i::ndown][0:n_rep] for i in range(ndown)]).mean(0)

    return ts_ds


class OutputSource():

      pulse = None  # Pulse Number
      SNR   = None  # SNR of pulse
      DM    = None  # DM (pc/cm3) of pulse
      time  = None  # Time at which pulse ocurred
      dtau  = None  # Temporal resolution of time series
      dnu   = None  # Spectral resolution
      nu    = None  # Central Observing Frequency
      mean  = None  # Mean in the time series
      rms   = None  # RMS in the time series

      formatter = "{0.pulse:07d}    {0.SNR:10.6f}     {0.DM:10.4f}     {0.time:10.6f} "+\
                 "     {0.dtau:10.6f}     {0.dnu:.4f}     {0.nu:.4f}    {0.mean:.5f}"+\
                 "    {0.rms:0.5f}\n " 

      def __str__(self):
          return self.formatter.format(self)

def savitzky_golay(y, window_size, order, deriv=0):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter

    This implementation is based on [1]_.

    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.

    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute
        (default = 0 means only smoothing)

    Returns
    -------
    y_smooth : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).

    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    Examples
    --------
    >>> t = np.linspace(-4, 4, 500)
    >>> y = np.exp(-t ** 2)
    >>> np.random.seed(0)
    >>> y_noisy = y + np.random.normal(0, 0.05, t.shape)
    >>> y_smooth = savitzky_golay(y, window_size=31, order=4)
    >>> print np.rms(y_noisy - y)
    >>> print np.rms(y_smooth - y)

    References
    ----------
    .. [1] http://www.scipy.org/Cookbook/SavitzkyGolay
    .. [2] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [3] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")

    if window_size % 2 != 1 or window_size < 1:
        #raise TypeError("window_size size must be a positive odd number")
        window_size += 1

    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")

    order_range = range(order + 1)

    half_window = (window_size - 1) // 2

    # precompute coefficients
    b = np.mat([[k ** i for i in order_range]
                for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv]

    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])

    y = np.concatenate((firstvals, y, lastvals))

    return np.convolve(y, m, mode='valid')

def RFI(sp,std):
    sp[:, np.where( np.abs(sp[:,:].mean(0)-np.median(sp[:,:].mean(0)) ) > std/np.sqrt(sp.shape[0]) )   ] = np.median(sp)
    sp[   np.where( np.abs(sp[:,:].mean(1)-np.median(sp[:,:].mean(1)) ) > std/np.sqrt(sp.shape[1]) ), :] = np.median(sp)
    return sp

if __name__ == '__main__':
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    fpp   =  44 #spectrogram per processer you want, limited mainly by 64GB memory per node
    nodes =  2 #the number of node requensted in sh
    pps   =  16 #processer per node requensted in sh
    numberofFiles=fpp*nodes*pps #totalnumberofspec = 6895.

    maxpw = 10 #Maximum pulse width to search in seconds. default = 1 s.
    thresh= 5.0 #SNR cut off

    fn   = sorted(glob.glob('05*.npy')) 
    tInt = np.load('tInt.npy')
    std = np.zeros((fpp, nodes*pps))

    pol = 1

    DMstart= 1000. #1.0 #initial DM trial
    DMend  = 5000. #90.0 #finial  DM trial

    fcl = 6000+7000#low frequency cut off
    fch = fcl+3343#high frequency cut off
    bp = np.zeros((fpp, nodes*pps, fch-fcl))

    npws = int(np.round(np.log2(maxpw/tInt)))+1 # +1 Due to in range(y) it goes to y-1 only

    spect=np.load(fn[0],mmap_mode='r')
    spectarray = np.zeros((fpp,spect.shape[0],spect.shape[1])) # X and Y are merged already after bandpass

    #cobimed spectrogram and remove background
    for i in range(fpp):
        spectarray[i,:,:] = np.load(fn[rank*fpp+i])
        bp[i,rank,:]= np.median(spectarray[i,:,:], 0)
 
    bpall=bp*0. #initiate a 4 hour blank std
    comm.Allreduce(bp,bpall,op=MPI.SUM) #merge the 4 hour std from all processor

    if rank == 0:
        np.save('bpallp%.1i' % pol, bpall)

    bp = np.median(np.median(bpall[:,:,:],0),0)

    bp = savitzky_golay(bp,1510,2)

    for i in range(fpp):
        spectarray[i,:,:] /= bp
        #spectarray[i,:,:] /= np.median(spectarray[i,:,:])
        std[i,rank]=spectarray[i,:,:].std()

    stdall=std*0. #initiate a 4 hour blank std
    comm.Allreduce(std,stdall,op=MPI.SUM) #merge the 4 hour std from all processor

    std = np.median(stdall)

    if rank ==0:
        np.save('stdall%.1i' % pol, stdall)

    #'''
    #RFI removal
    for i in range(fpp):
        spectarray[i] = RFI(spectarray[i],5.*std)
        spectarray[i] /= spectarray[i].mean()
    #'''

    if  pol < 4:
        if pol==0:
            freq=np.load('freq1.npy')
        else: 
            freq=np.load('freq2.npy')
        freq /= 10**6
        cent_freq = np.median(freq)
        BW   = freq.max()-freq.min()
        DM=DMstart
    
        txtsize=np.zeros((npws,2),dtype=np.int32) #fileno = txtsize[ranki,0], pulse number = txtsize[ranki,1],ranki is the decimated order of 2
        txtsize[:,0]=1 #fileno star from 1

        tbmax=0 #repeat control, if dedispersion time series are identical, skip dedispersion calculation
        while DM < DMend:
            if DM >=1000.: dDM = 1.
            else: dDM = 0.1
            tb=np.round((delay2(freq,DM)/tInt)).astype(np.int32)
            if tb.max()-tbmax==0:#identical dedispersion time series checker
                tbmax=tb.max()
                if DM > 0.:
                    DM+=dDM
                continue
            tbmax=tb.max()

            ts=np.zeros((tb.max()+numberofFiles*np.load(fn[0],mmap_mode='r').shape[0]))
            for freqbin in range(len(freq)): 
                for i in range(fpp):
                    ts[tb.max()-tb[freqbin] + (rank*fpp+i)*spect.shape[0] :tb.max()-tb[freqbin] + (rank*fpp+i+1)*spect.shape[0] ] += spectarray[i,:,freqbin]

            tstotal=ts*0#initiate a 4 hour blank time series
            comm.Allreduce(ts,tstotal,op=MPI.SUM)#merge the 4 hour timeseries from all processor
            tstotal = tstotal[tb.max():len(tstotal)-tb.max()]#cut the dispersed time lag


            #'''
            # save the time series around the Pulsar's DM
            if rank == 0:
                if np.abs(DM - 9.1950) <= dDM:
                    print 'DM=',DM
                    np.save('ts_pol%.1i_DMx100_%.6i' % (pol,DM*100),tstotal)
            #'''

            #"""#search for signal with decimated timeseries
            if rank<npws:#timeseries is ready for signal search
                ranki=rank
                filename = "ppc_SNR_pol_%.1i_td_%.2i_no_%.05i.txt" % (pol,ranki,txtsize[ranki,0])
                outfile = open(filename,'a')
                ndown = 2**ranki #decimate the time series
                sn,mean,rms = Threshold(Decimate_ts(tstotal,ndown),thresh,niter=0)
                ones = np.where(sn!=-1)[0]
                for one in ones:# Now record all pulses above threshold
                    pulse = OutputSource()
                    txtsize[ranki,1] += 1
                    pulse.pulse = txtsize[ranki,1]
                    pulse.SNR = sn[one]
                    pulse.DM = DM
                    pulse.time = one*tInt*ndown
                    pulse.dtau = tInt*ndown
                    pulse.dnu = freq[1]-freq[0]
                    pulse.nu = cent_freq
                    pulse.mean = mean
                    pulse.rms = rms
                    outfile.write(pulse.formatter.format(pulse)[:-1]) 
                    if txtsize[ranki,1] >200000*txtsize[ranki,0]:
                        outfile.close()
                        txtsize[ranki,0]+=1
                        filename = "ppc_SNR_pol_%.1i_td_%.2i_no_%.05d.txt" % (pol,ranki,txtsize[ranki,0])
                        outfile = open(filename,'a')

            DM+=dDM # End of DM loop
