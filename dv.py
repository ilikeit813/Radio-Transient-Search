from mpi4py import MPI
import disper
import sys
import numpy as np
import glob
import os
import time
import sys
import getopt
from apputils import forceIntValue

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

def snr(x):return (x-x.mean())/x.std()

def bpf(x, windows = 40):
    bp = savitzky_golay(x,windows,1)
    x2 = x / bp
    mask = np.where(snr(x2)>1)[0]
    mask2= np.zeros(x.shape[0])
    mask2[mask] = 1.
    y = np.ma.array(x, mask = mask2)
    bp = savitzky_golay(y,windows,1)
    fit = np.ma.polyfit(np.arange(len(y)),y,4)
    p = np.poly1d(fit)(np.arange(len(y)))[mask]
    bp = x
    bp[mask] = np.poly1d(fit)(np.arange(len(y)))[mask]
    bp = savitzky_golay(bp,windows,2)
    return bp

def fold(t,period,T0=0):
    time = np.arange(len(t))
    epoch = np.floor( 1.*(time - T0)/period )
    phase = 1.*(time - T0)/period - epoch
    foldt = t[np.argsort(phase)]
    return Decimate_ts(foldt, 1.*len(t)/period )


def RFImask(spr):#sp has shape[:,:]
    x = np.where(abs(spr.mean(1))>np.sort(spr.mean(1))[spr.shape[0]/2]+np.sort(spr.mean(1))[spr.shape[0]/2]-np.sort(spr.mean(1))[1])
    y = np.where(abs(spr.mean(0))>np.sort(spr.mean(0))[spr.shape[1]/2]+np.sort(spr.mean(0))[spr.shape[1]/2]-np.sort(spr.mean(0))[1])
    return [x[0],y[0]]

def massagesp(spectrometer, windows_x=43,windows_y=100):
    bp = bpf(spectrometer.mean(0),windows_x)
    spectrometer /= bp
    bl = bpf(spectrometer.mean(1),windows_y)
    spectrometer = (spectrometer.T - bl).T
    mask  = np.array ( RFImask(spectrometer) )
    mask2 = np.zeros((spectrometer.shape))
    mask2[mask[0],:] = 1.
    mask2[:,mask[1]] = 1.
    temp_spec = np.ma.array(spectrometer, mask = mask2 )
    mean = temp_spec.mean()
    spectrometer[mask[0],:] = mean
    spectrometer[:,mask[1]] = mean
    spectrometer -= mean
    return spectrometer



if __name__ == '__main__':
    #fcl = 360/4
    #fch = 3700/4
    fcl = 0
    fch = 4095
    pol = 0
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    fpp   =  264/12 #spectrogram per processer you want, limited mainly by 64GB memory per node (32GB Hokieone)
    nodes =  2 #the number of node requensted in sh
    pps   =  6 #processer per node requensted in sh
    numberofFiles=fpp*nodes*pps #totalnumberofspec = 6895.

    maxpw = 600 #Maximum pulse width to search in seconds. default = 1 s.
    thresh= 5.0 #SNR cut off

    # CCY - This modification allows specifying the original time series data filepath and fcl and fch
    # parameters from the command line.  This gives more flexibility in that it allows dv.py to find any
    # pattern of frame files and allows specifying the frequency band without having to make hard-coded
    # modifications to dv.py.
    #
    szShortOpts = 'lut'
    szLongOpts = ['lower', 'upper', 'tuning']
    filepath = getopt.getopt(sys.argv[1:], szShortOpts, szLongOpts)[1][0]
    # Make sure that a path to the original data file has been provided.
    if len(filepath) == 0:
        print('Path to the original data file must be provided')
        exit(1)
    # end if
    filename = os.path.basename(os.path.splitext(filepath)[0])
    (cmdLnOpts, cmdLnParams) = getopt.getopt(sys.argv[1:], szShortOpts, szLongOpts)[0]
    for index in range(len(cmdLnOpts)):
        # Parse the command-line options.
        if cmdLnOpts[index] in [szShortOpts[0], szLongOpts[0]]:
            fcl = forceIntValue(cmdLnParams[index], 0, 4095)
        elif cmdLnOpts[index] in [szShortOpts[1], szLongOpts[1]]:
            fch = forceIntValue(cmdLnParams[index], 0, 4095)
        elif cmdLnOpts[index] in [szShortOpts[2], szLongOpts[2]]:
            if cmdLnParams[index] in ['low', 'high']:
               if cmdLnParams[index] == 'low':
                  pol = 0
               else
                  pol = 1
               # end if
            else
               print('Allowed tuning specifications: \'low\', \'high\'')
               exit(1)
            # end if
        else
            print('UNKNOWN OPTION: {opt}'.format(opt=cmdLnOpts[index]))
            exit(1)
        # end if
    # end for index in range(len(cmdLnOpts))


    fn   = sorted(glob.glob('{filename}*.npy'.format(filename=filename)) 
    tInt = np.load('tInt.npy')

    DMstart =  0 #1.0 #initial DM trial
    DMend   =  5000 #90.0 #finial  DM trial
    npws = int(np.round(np.log2(maxpw/tInt)))+1 # +1 Due to in range(y) it goes to y-1 only

    spect=np.load(fn[0],mmap_mode='r')[:,:,fcl:fch]
    spectarray = np.zeros((fpp,spect.shape[0],spect.shape[2])) # X and Y are merged already after bandpass

    #cobimed spectrogram and remove background
    for i in range(fpp):
        print '1',(np.load(fn[rank*fpp+i])[:,pol,fcl:fch]).shape
        print '2',massagesp( np.load(fn[rank*fpp+i])[:,pol,fcl:fch] ).shape
        spectarray[i,:,:] = massagesp( np.load(fn[rank*fpp+i])[:,pol,fcl:fch], 10, 50 )

    np.save('spectarray%.2i' % rank, spectarray)
    #sys.exit()

    if  pol < 4:
        if pol==0:
            freq=np.load('freq1.npy')[fcl:fch]
        else: 
            freq=np.load('freq2.npy')[fcl:fch]
        freq /= 10**6
        cent_freq = np.median(freq)
        BW   = freq.max()-freq.min()
        DM = DMstart
    
        txtsize=np.zeros((npws,2),dtype=np.int32) #fileno = txtsize[ranki,0], pulse number = txtsize[ranki,1],ranki is the decimated order of 2
        txtsize[:,0]=1 #fileno star from 1


        DMtrials = DMstart # 0
        if rank == 0:
            while DM < DMend:
                #dDM = disper.dDMi(DMtrial = 1.*DM, nuCenteralMHz = 1.*cent_freq, channelMHz = freq[1]-freq[0], BMHz = freq[-1]-freq[0], SSratio = 0.8, temporal_resol = 1.*tInt)
                #print dDM
                if DM < 1000:
                    dDM = 0.1
                elif DM >= 1000:
                    dDM = 1.
                DM += dDM
                DMtrials = np.append(DMtrials,DM)

        DMtrials = comm.bcast(DMtrials,root =0)

        for DM in DMtrials:
            #print 'rank = ',rank, ' DM trial = ',DM
            tb=np.round((delay2(freq,DM)/tInt)).astype(np.int32)

            ts=np.zeros((tb.max()+numberofFiles*np.load(fn[0],mmap_mode='r').shape[0]))
            for freqbin in range(len(freq)): 
                for i in range(fpp):
                    ts[tb.max()-tb[freqbin] + (rank*fpp+i)*spect.shape[0] :tb.max()-tb[freqbin] + (rank*fpp+i+
                        1)*spect.shape[0] ] += spectarray[i,:,freqbin]

            tstotal=ts*0#initiate a 4 hour blank time series
            comm.Allreduce(ts,tstotal,op=MPI.SUM)#merge the 4 hour timeseries from all processor
            tstotal = tstotal[tb.max():len(tstotal)-tb.max()]#cut the dispersed time lag


            '''
            # save the time series around the Pulsar's DM
            if rank == 0:
                if np.abs(DM - 10.922) <= dDM:
                    print 'DM=',DM
                    np.save('ts_pol%.1i_DMx100_%.6i' % (pol,DM*100),tstotal)
            sys.exit()
            '''

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
                    # end if
                # end for
            # end if
        # end for
    # end if
