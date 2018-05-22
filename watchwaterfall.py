import sys
import numpy as np
import glob
import os
import time
import sys
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from apputils import forceIntValue

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

    # Create the filter kernel.  We only want to create this _ONCE_ for any given set of parameters.
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")

    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")

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
    #sp[:, np.where( np.abs(sp[:,:].mean(0)-np.median(sp[:,:].mean(0)) ) > std/np.sqrt(sp.shape[0]) )   ] = np.median(sp)
    #sp[   np.where( np.abs(sp[:,:].mean(1)-np.median(sp[:,:].mean(1)) ) > std/np.sqrt(sp.shape[1]) ), :] = np.median(sp)

    spMean0 = sp.mean(0)
    spMean1 = sp.mean(1)
    spMedian = np.median(sp)
    sp[:, np.where( np.abs(spMean0 - np.median(spMean0)) > std/np.sqrt(sp.shape[0]) )] = spMedian
    sp[np.where( np.abs(spMean0 - np.median(spMean0)) > std/np.sqrt(sp.shape[0]) ), :] = spMedian

    return sp

def snr(a):
    return (a-a.mean() )/a.std()



# == Main routine ==
#
# CCY - The following modification to the command-line adds options and flexibility to bandpasscheck.py
# to allow it to be more scriptable in the radio transient workflow.
#
usage="Usage: %prog [options] <waterfall file>"
cmdlnParser = OptionParser(usage=usage)
cmdlnParser.add_option('-a', '--low-tuning-lower',dest='lowTuneLower', default=0, type='int',
                        action='store',
                        help='Lower FFT index (between 0 and 4095) for the low tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-b', '--low-tuning-upper',dest='lowTuneUpper', default=4095, type='int',
                        action='store',
                        help='Upper FFT index (between 0 and 4095) for the low tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-y', '--high-tuning-lower',dest='highTuneLower', default=0, type='int',
                        action='store',
                        help='Lower FFT index (between 0 and 4095) for the high tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-z', '--high-tuning-upper',dest='highTuneUpper', default=4095, type='int',
                        action='store',
                        help='Upper FFT index (between 0 and 4095) for the high tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                       action="store",
                       help="Working directory path.", metavar="PATH")
(cmdlnOpts, cmdlnArgs) = cmdlnParser.parse_args()
workDir = cmdlnOpts.workDir
# Check that a combined waterfall file has been provided and that it exists.
if len(cmdlnArgs[0]) > 0:
   if not os.path.exists(cmdlnArgs[0]):
      print('Cannot find the combined waterfall file \'{path}\''.format(path=cmdlnArgs[0]))
      exit(1)
   # endif
else:
   print('A path to the combined waterfall file must be given.')
   exit(1)
# endif

# Force the FFT index values to be between 0 and 4095.
Lfcl = forceIntValue(cmdlnOpts.lowTuneLower, 0, 4095)
Lfch = forceIntValue(cmdlnOpts.lowTuneUpper, 0, 4095)
Hfcl = forceIntValue(cmdlnOpts.highTuneLower, 0, 4095)
Hfch = forceIntValue(cmdlnOpts.highTuneUpper, 0, 4095)
# Check that the lower FFT indice are less than the upper FFT indices.
if not Lfch > Lfcl:
   print('Low tuning lower FFT index must be less than the upper FFT index.')
   exit(1)
# endif
if not Hfch > Hfcl:
   print('High tuning lower FFT index must be less than the upper FFT index.')
   exit(1)
# endif

# CCY - Because we are now allowing the spectrogram to have a bandpass selection for the different
# tunings, we must explicitly pull out the bandpass for each tuning separately.
#
rawData = np.load(cmdlnArgs[0])
lowSpectrogram = rawData[:,0,Lfcl:(Lfch + 1)]
highSpectrogram = rawData[:,1,Hfcl:(Hfch + 1)]
# Correct bandpass and baseline for the low tuning
lowBandpass = savitzky_golay( np.median(lowSpectrogram, 0), 151, 2 )
lowSpectrogram = lowSpectrogram - lowBandpass
lowBaseline = savitzky_golay( np.median(lowSpectrogram, 1), 151, 2 )
lowSpectrogram = (lowSpectrogram.T - lowBaseline).T

# Correct bandpass and baseline for the high tuning.
highBandpass = savitzky_golay( np.median(highSpectrogram, 0), 151, 2 )
highSpectrogram[:,:] = highSpectrogram[:,:] - highBandpass
highBaseline = savitzky_golay( np.median(highSpectrogram, 1), 151, 2 )
highSpectrogram = (highSpectrogram.T - highBaseline).T


lowSTD = lowSpectrogram.std()
lowSpectrogram = RFI(lowSpectrogram, 5.0*lowSTD)
lowSpectrogram = snr(lowSpectrogram)
lowSTD = lowSpectrogram.std()
lowSpectrogram[ np.where( (abs(lowSpectrogram) > 3.0*lowSTD) ) ] = lowSpectrogram.mean()

highSTD = highSpectrogram.std()
highSpectrogram = RFI(highSpectrogram, 5.0*highSTD)
highSpectrogram = snr(highSpectrogram)
highSTD = highSpectrogram.std()
highSpectrogram[ np.where( (abs(highSpectrogram) > 3.0*highSTD) ) ] = highSpectrogram.mean()

#bandpass
#bp = 0.*sp[0,:,:]
#baseline
#bl = 0.*sp[:,:,0]
#bp[:,:]= np.median(sp, 0)
#bp[0,] = savitzky_golay(bp[0,:],151,2)
#bp[1,] = savitzky_golay(bp[1,:],111,2)
#correct the bandpass
#for tuning in (0,1):
#    sp[:,tuning,:] = sp[:,tuning,:]-bp[tuning,]
#bl[:,:]= np.median(sp, 2)
#bl[:,0] = savitzky_golay(bl[:,0],151,2)
#bl[:,1] = savitzky_golay(bl[:,1],151,2)

#correct the baseline
#for tuning in (0,1):
#    sp[:,tuning,:] = (sp[:,tuning,:].T - bl[:,tuning].T).T

#for tuning in (0,1):
#    sp[:,tuning,:] = RFI(sp[:,tuning,:],5.*sp[:,tuning,:].std())

#for tuning in (0,1):
#    sp[:,tuning,:] = snr(sp[:,tuning,:])

#for tuning in (0,1):
#    sp[:,tuning,:][ np.where( ( abs(sp[:,tuning,:]) > 3.*sp[:,tuning,:].std()) )] = sp[:,tuning,:].mean()


plt.imshow(lowSpectrogram.T, cmap='Greys_r', origin = 'low', aspect = 'auto')
plt.suptitle('Spectrogram Low Tuning', fontsize = 30)
plt.xlabel('Time (14 sec)',fontdict={'fontsize':16})
plt.ylabel('Frequency (4.79 kHz)',fontdict={'fontsize':14})
plt.colorbar().set_label('STDs About Mean',size=18)
plt.savefig('{dir}/lowspectrogram'.format(dir=workDir))
plt.clf()

plt.imshow(highSpectrogram.T, cmap='Greys_r', origin = 'low', aspect = 'auto')
plt.suptitle('Spectrogram High Tuning', fontsize = 30)
plt.xlabel('Time (14 sec)',fontdict={'fontsize':16})
plt.ylabel('Frequency (4.79 kHz)',fontdict={'fontsize':14})
plt.colorbar().set_label('STDs About Mean',size=18)
plt.savefig('{dir}/highspectrogram'.format(dir=workDir))
plt.clf()

exit(0)
