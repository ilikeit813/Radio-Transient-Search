from mpi4py import MPI
import disper
import sys
import numpy as np
import glob
import os
import time
import sys
import getopt
from optparse import OptionParser
from apputils import forceIntValue, clipValue, procMessage, atProcMessage

SGkernels = None  # Dictionary of kernels to use for the Savitsky-Golay smoothing.


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
   global SGkernels

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

   # Only build the kernel if the kernel for the specified window_size does not exist in out dictionary.
   half_window = (window_size - 1) // 2
   if (window_size, order, deriv) not in SGkernels: 
      order_range = range(order + 1)
      # precompute coefficients
      b = np.mat([[k ** i for i in order_range]
               for k in range(-half_window, half_window + 1)])
      m = np.linalg.pinv(b).A[deriv]
      SGkernels[window_size, order, deriv] = m
   else:
      m = SGkernels[window_size, order, deriv]
   # endif

   # pad the signal at the extremes with
   # values taken from the signal itself
   firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
   lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])

   y = np.concatenate((firstvals, y, lastvals))

   return np.convolve(y, m, mode='valid')
# end savitzky_golay()


def snr(x): 
   return (x - x.mean())/x.std()
# end snr()


def bpf(x, windows = 40):
   bp = savitzky_golay(x,windows,1)
   mask = np.where(snr(x / bp) > 1)[0]
   mask2= np.zeros(x.shape[0])
   mask2[mask] = 1.
   y = np.ma.array(x, mask = mask2)
   s = np.arange(len(y))
   fit = np.ma.polyfit(s,y,4)
   bp = x
   bp[mask] = np.poly1d(fit)(s)[mask]
   return savitzky_golay(bp,windows,2)
# end bpf()


def RFImask(spr):#sp has shape[:,:]
   #x = np.where(abs(spr.mean(1))>np.sort(spr.mean(1))[spr.shape[0]/2]+np.sort(spr.mean(1))[spr.shape[0]/2]-np.sort(spr.mean(1))[1])
   #y = np.where(abs(spr.mean(0))>np.sort(spr.mean(0))[spr.shape[1]/2]+np.sort(spr.mean(0))[spr.shape[1]/2]-np.sort(spr.mean(0))[1])

   sprMean1 = spr.mean(1)
   sprMean1.sort()
   sprMean0 = spr.mean(0)
   sprMean0.sort()
   sprMedian1 = sprMean1[sprMean1.shape[0]/2]
   sprMedian0 = sprMean0[sprMean0.shape[0]/2]
   #x = np.where(abs(sprMean1)>2*sprMedian1 - np.sort(sprMean1)[1])
   #y = np.where(abs(sprMean0)>2*sprMedian0 - np.sort(sprMean0)[1])
   x = np.nonzero( np.fabs(sprMean1) > (2*sprMedian1 - sprMean1[1]) )[0]
   y = np.nonzero( np.fabs(sprMean0) > (2*sprMedian0 - sprMean0[1]) )[0]
 
   return [x, y]
# end RFImask()

def massagesp(spectrometer, windows_x=43,windows_y=100):
   bp = bpf(spectrometer.mean(0),windows_x)
   spectrometer /= bp
   bl = bpf(spectrometer.mean(1),windows_y)
   spectrometer = (spectrometer.T - bl).T
   indices = RFImask(spectrometer)
   mask = np.zeros((spectrometer.shape), dtype=np.int32)
   mask[indices[0],:] = 1
   mask[:,indices[1]] = 1
   mean = np.ma.array(spectrometer, mask = mask ).mean()
   spectrometer[indices[0],:] = mean
   spectrometer[:,indices[1]] = mean
   spectrometer -= mean
   return spectrometer
# end massagesp()


# main routine
#
def main_routine():
   global SGkernels
   SGkernels = dict()

   comm  = MPI.COMM_WORLD
   rank  = comm.Get_rank()
   nProcs = comm.Get_size()


   # Setup command-line parsing.
   usage="USAGE: %prog [options] <radio data filepath>" 
   cmdlnParser = OptionParser(usage=usage)
   cmdlnParser.add_option("-a","--low-tuning-lower",dest="Lfcl", default=0, type="int",
                           action="store",
                           help="Lower FFT index (between 0 and 4095) for the low tuning.", 
                           metavar="INDEX")
   cmdlnParser.add_option("-b","--low-tuning-upper",dest="Lfch", default=4095, type="int",
                           action="store",
                           help="Upper FFT index (between 0 and 4095) for the low tuning.", metavar="INDEX")
   cmdlnParser.add_option("-y","--high-tuning-lower",dest="Hfcl", default=0, type="int",
                           action="store",
                           help="Lower FFT index (between 0 and 4095) for the high tuning.", 
                           metavar="INDEX")
   cmdlnParser.add_option("-z","--high-tuning-upper",dest="Hfch", default=4095, type="int",
                           action="store",
                           help="Upper FFT (between 0 and 4095) index for the high tuning.", 
                           metavar="INDEX")
   cmdlnParser.add_option("-p","--bandpass-window",dest="BPwindow", default=10, type="int",
                           action="store",
                           help="Window size for bandpass smoothing.", 
                           metavar="SIZE")
   cmdlnParser.add_option("-l","--baseline-window",dest="BLwindow", default=50, type="int",
                           action="store",
                           help="Window size for baseline smoothing.", 
                           metavar="SIZE")
   cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                          action="store",
                          help="Working directory path.", metavar="PATH")
   
   # Parse command-line for FFT bandpass indices and the radio data file path.
   (cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()
   filepath = cmdlnParams[0]
   Lfcl = forceIntValue(cmdlnOpts.Lfcl, 0, 4095)
   Lfch = forceIntValue(cmdlnOpts.Lfch, 0, 4095)
   Hfcl = forceIntValue(cmdlnOpts.Hfcl, 0, 4095)
   Hfch = forceIntValue(cmdlnOpts.Hfch, 0, 4095)
   BPwindow = abs(cmdlnOpts.BPwindow)
   BLwindow = abs(cmdlnOpts.BLwindow)
   workDir = cmdlnOpts.workDir

   # Validate command-line arguments.
   if len(filepath) == 0:
      print('Path to the original data file must be provided')
      exit(1)
   # endif
   if Lfch <= Lfcl:
      print('Low tuning lower FFT index must be less than the upper FFT index.')
      exit(1)
   # endif
   if Hfch <= Hfcl:
      print('High tuning lower FFT index must be less than the upper FFT index.')
      exit(1)
   # endif

   # Collect the detailed spectrogram filenames.
   filename = os.path.basename(os.path.splitext(filepath)[0])
   filenamePattern = '{dir}/master{name}*.npy'.format(dir=workDir, name=filename)
   tagSplit = '{dir}/master{name}'.format(dir=workDir, name=filename)
   spectFilenames = None
   if rank == 0:
      spectFilenames = sorted(glob.glob(filenamePattern))
   # endif
   spectFilenames = comm.bcast(spectFilenames, root=0)

   # If we have files, process them.
   numFiles = len(spectFilenames)
   if numFiles > 0:

      # Preallocate the low and high tuning spectrogram segments.
      sampleSpect = np.load(spectFilenames[0])
      lowSpect = np.zeros((sampleSpect.shape[0], Lfch - Lfcl + 1))
      highSpect = np.zeros((sampleSpect.shape[0], Hfch - Hfcl + 1))

      # Determine which files will be handled by this process.
      #numFiles = min(numFiles, 100) # Profiling short-circuit.
      filesPerProc = int(numFiles/nProcs)
      fileIndices = np.arange(filesPerProc) + filesPerProc*rank
      if rank == 0:
         fileIndices = np.concatenate((fileIndices, np.arange(filesPerProc*nProcs, numFiles)), axis=0)
      # endif

      # Iterate through the files assigned to this process.
      fileIndex = filesPerProc*rank
      filesCount = 0
      for fileIndex in fileIndices:
         # Perform the data smoothing, RFI cleaning, and bandpass filtering on the spectrograms for the low
         # and high tunings, separately.
         procMessage(''.join(['Smoothing, RFI cleaning, and bandpass filtering '.format(rank=rank),
                              'file {index} of {total}'.format(index=fileIndex + 1, total=numFiles)]) )
         spectArray = np.load(spectFilenames[fileIndex])
         lowSpect[:,:] = massagesp( spectArray[:,0,:], BPwindow, BLwindow )[:, Lfcl:(Lfch + 1) ]
         highSpect[:,:] = massagesp( spectArray[:,1,:], BPwindow, BLwindow )[:, Hfcl:(Hfch + 1) ]

         # Save the transpose of the new smoothed, cleaned, and filtered spectrograms.
         # CCY - the transpose is need for memory efficiency during the de-dispersion phase, which comes
         # after the phase that uses this module.
         (root, ext) = os.path.splitext(spectFilenames[fileIndex])
         tag = root.split(tagSplit)[-1]
         np.save('{dir}/master0_{name}{tag}'.format(dir=workDir, name=filename, tag=tag), lowSpect.T)
         np.save('{dir}/master1_{name}{tag}'.format(dir=workDir, name=filename, tag=tag), highSpect.T) 

         # Delete the old spectrogram file.
         os.remove(spectFilenames[fileIndex])

      # end for fileIndex in fileIndices:
   else:
      procMessage('No files to process.')
   # endif
# end main_routine()



if __name__ == '__main__' :
   main_routine()
   exit(0)
# endif
