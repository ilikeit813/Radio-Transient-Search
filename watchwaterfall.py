import sys
import numpy as np
import glob
import os
import time
import sys
import matplotlib.pyplot as plt

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
    sp[:, np.where( np.abs(sp[:,:].mean(0)-np.median(sp[:,:].mean(0)) ) > std/np.sqrt(sp.shape[0]) )   ] = np.median(sp)
    sp[   np.where( np.abs(sp[:,:].mean(1)-np.median(sp[:,:].mean(1)) ) > std/np.sqrt(sp.shape[1]) ), :] = np.median(sp)
    return sp

def snr(a):
    return (a-a.mean() )/a.std()


sp = np.load('waterfall.npy')
#bandpass
bp = 0.*sp[0,:,:]
#baseline
bl = 0.*sp[:,:,0]
bp[:,:]= np.median(sp, 0)
bp[0,] = savitzky_golay(bp[0,],151,2)
bp[1,] = savitzky_golay(bp[1,],111,2)
#correct the bandpass
for tuning in (0,1):
    sp[:,tuning,:] = sp[:,tuning,:]-bp[tuning,]
bl[:,:]= np.median(sp, 2)
bl[:,0] = savitzky_golay(bl[:,0],151,2)
bl[:,1] = savitzky_golay(bl[:,1],151,2)
#correct the baseline
for tuning in (0,1):
    sp[:,tuning,:] = (sp[:,tuning,:].T - bl[:,tuning].T).T

for tuning in (0,1):
    sp[:,tuning,:] = RFI(sp[:,tuning,:],5.*sp[:,tuning,:].std())

for tuning in (0,1):
    sp[:,tuning,:] = snr(sp[:,tuning,:])

for tuning in (0,1):
    sp[:,tuning,:][ np.where( ( abs(sp[:,tuning,:]) > 3.*sp[:,tuning,:].std()) )] = sp[:,tuning,:].mean()

#'''
plt.imshow(sp[:,0,:].T, cmap='Greys_r', origin = 'low', aspect = 'auto')
plt.suptitle('Spectrogram Low Tuning', fontsize = 30)
plt.xlabel('Time (14 sec)',fontdict={'fontsize':16})
plt.ylabel('Frequency (4.7 kHz)',fontdict={'fontsize':14})
plt.colorbar().set_label('std',size=18)
#plt.show()
plt.savefig('Low')
#'''
plt.clf()

plt.imshow(sp[:,1,:].T, cmap='Greys_r', origin = 'low', aspect = 'auto')
plt.suptitle('Spectrogram High Tuning', fontsize = 30)
plt.xlabel('Time (14 sec)',fontdict={'fontsize':16})
plt.ylabel('Frequency (4.7 kHz)',fontdict={'fontsize':14})
plt.colorbar().set_label('std',size=18)
#plt.show()
plt.savefig('High')
