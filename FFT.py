#script of Python script to dedisperse a file, save each DM trials time series
#Codes from Sean, modified by Jamie
#V3, super cluster technique

from mpi4py import MPI
import sys
import numpy
import numpy as np
import glob
import os
import time
import sys
import getopt
import drx
import errors
import dp
def DMs(DMstart,DMend,dDM):
    """
    Calculate the number of DMs searched between DMstart and DMend, with spacing dDM * DM.

    Required:

    DMstart   - Starting Dispersion Measure in pc cm-3
    DMend    - Ending Dispersion Measure in pc cm-3
    dDM    - DM spacing in pc cm-3
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
    DM  = None  # DM (pc/cm3) of pulse
    time  = None  # Time at which pulse ocurred
    dtau  = None  # Temporal resolution of time series
    dnu   = None  # Spectral resolution
    nu  = None  # Central Observing Frequency
    mean  = None  # Mean in the time series
    rms   = None  # RMS in the time series

    formatter = "{0.pulse:07d}  {0.SNR:10.6f}    {0.DM:10.4f}    {0.time:10.6f} "+\
                "    {0.dtau:10.6f}  {0.dnu:.4f}     {0.nu:.4f} {0.mean:.5f}"+\
                "   {0.rms:0.5f}\n " 

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


def main(args):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    tunes = int(getopt.getopt(args,':')[1][1]) # 0: low, 1: hig, read in from your command line
    nodes =  1 #89 #the number of node requensted in sh
    pps   =  16 #processer per node requensted in sh

    fcl =  6000+7000 #low frequency cut off
    fch =  fcl+3343 #high frequency cut off

    nChunks = 3000 #the temporal shape of a file.
    LFFT = 4096*16 #Length of the FFT.4096 is the size of a frame readed.
    nFramesAvg = 1*4*LFFT/4096 # the intergration time under LFFT, 4 = beampols = 2X + 2Y (high and low tunes)

    offset_i = rank # range to 4309


    fn = sorted(glob.glob('05*.npy'))
    j = np.zeros((len(fn)))
    for i in range(len(fn)):
        j[i] = fn[i][30:39]

    k=[]
    for i in range(len(j)-1):
        if j[i+1] != j[i] + nChunks*nFramesAvg:
            k.append(j[i]+nChunks*nFramesAvg)

    l = len(k)/(nodes*pps)+1

    if rank == 0 :
        fh = open(getopt.getopt(args,':')[1][0], "rb")
        nFramesFile = os.path.getsize(getopt.getopt(args,':')[1][0]) / drx.FrameSize #drx.FrameSize = 4128
        junkFrame = drx.readFrame(fh)
        srate = junkFrame.getSampleRate()
        fh.seek(-drx.FrameSize, 1)
        beam,tune,pol = junkFrame.parseID()
        beams = drx.getBeamCount(fh)
        tunepols = drx.getFramesPerObs(fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        nFrames = nFramesAvg*nChunks
        centralFreq1 = 0.0
        centralFreq2 = 0.0
        for i in xrange(4):
            junkFrame = drx.readFrame(fh)
            b,t,p = junkFrame.parseID()
            if p == 0 and t == 0:
                try:
                    centralFreq1 = junkFrame.getCentralFreq()
                except AttributeError:
                    centralFreq1 = dp.fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
            elif p == 0 and t == 2:
                try:
                    centralFreq2 = junkFrame.getCentralFreq()
                except AttributeError:
                    centralFreq2 = dp.fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
            else:
                pass
        fh.seek(-4*drx.FrameSize, 1)
        freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
        tInt = 1.0*LFFT/srate
        freq1 = (freq+centralFreq1)[fcl:fch]
        freq2 = (freq+centralFreq2)[fcl:fch]
        np.save('freq1',freq1)
        np.save('freq2',freq2)
        np.save('tInt',tInt)

if __name__ == "__main__":
        main(sys.argv[1:])
