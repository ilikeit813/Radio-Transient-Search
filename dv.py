from mpi4py import MPI
import disper
import sys
import numpy as np
import scipy.sparse as spsparse
import glob
import os
import time
import sys
import getopt
from optparse import OptionParser
from apputils import forceIntValue, clipValue, procMessage, procMessage
import drx


def DMs(DMstart,DMend,dDM):
   """
   Calculate the number of DMs searched between DMstart and DMend, with spacing dDM * DM.
   Required:
   DMstart   - Starting Dispersion Measure in pc cm-3
   DMend    - Ending Dispersion Measure in pc cm-3
   dDM      - DM spacing in pc cm-3
   """

   #NDMs = np.log10(float(DMend)/float(DMstart))/np.log10(1.0+dDM)
   NDMs = (DMend-DMstart)/dDM

   return int(np.round(NDMs))
# end DMs()


def scaledDelay(freq1, freq2):
   # Dispersion constant in MHz^2 s / pc cm^-3 = 4.148808e3

   return (4.148808e3)*(1.0/(freq1*freq1) - 1.0/(freq2*freq2))
   
# end scaledDelay()


def delay2(freq1, freq2, dm):
   """
   Calculate the relative delay due to dispersion over a given frequency
   range in Hz for a particular dispersion measure in pc cm^-3.  Return
   the dispersive delay in seconds.  Same as delay, but w.r.t to highest frequency.
   ***Used to simulate a dispersed pulse.***
   Required:
   freq - 1-D array of frequencies in MHz
   dm   - Dispersion Measure in pc cm-3
   """
   # Delay in s
   tDelay = dm*scaledDelay(freq1, freq2)

   return tDelay
# end delay2()


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
   clip   -  Clipping SNR threshold for values to leave out of a mean/rms calculation.  default = 3.
   niter   -  Number of iterations in mean/rms calculation.  default = 1.
   Usage:
   >>sn, mean, rms = Threshold(ts, *options*)
   """
   #  Calculate, robustly, the mean and rms of the time series.  Any values greater than 3sigma are left
   #  out of the calculation.  This keeps the mean and rms free from sturation due to large deviations.

   mean = np.mean(ts) 
   std  = np.std(ts)  
   pulseIndices = []

   if niter > 0:
      for i in range(niter):
         ones = np.where((ts-mean)/std < clip)[0]  # only keep the values less than 3sigma
         mean = np.mean(ts[ones])
         std  = np.std(ts[ones])
   if std > 0.0:
      SNR = (ts-mean)/std
      pulseIndices = np.nonzero( SNR < thresh )[0]
   else:
      SNR = np.full(shape=ts.shape, dtype=np.float64, fill_value=-1.0)

   return (SNR, mean, std, pulseIndices)
# end Threshold()

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
# end Decimate_ts()


class OutputSource():

   def __init__(self):
      self.pulse = None  # Pulse Number
      self.SNR   = None  # SNR of pulse
      self.DM   = None  # DM (pc/cm3) of pulse
      self.time  = None  # Time at which pulse ocurred
      self.dtau  = None  # Temporal resolution of time series
      self.dnu   = None  # Spectral resolution
      self.nu   = None  # Central Observing Frequency
      self.mean  = None  # Mean in the time series
      self.std   = None  # STD in the time series
   # end __init__() constructor.

   def str(self):
      str1 = '{pulse:10.4f}   {SNR:10.6f}    '.format(pulse=float(self.pulse), SNR=float(self.SNR))
      str2 = '{DM:10.4f}    {time:10.6f}     '.format(DM=float(self.DM), time=float(self.time))
      str3 = '{dtau:10.6f}    {dnu:.4f}    '.format(dtau=float(self.dtau), dnu=float(self.dnu))
      str4 = '{nu:.4f}   {mean:.5f}   '.format(nu=float(self.nu), mean=float(self.mean))
      str5 = '{std:0.5f}\n'.format(std=float(self.std))
      return ''.join([str1, str2, str3, str4, str5])
   # end str()

   def __str__(self):
      return self.str()
   # end __str__()
# end class OutputSource


def computeFrameIndexDelay(tInt, freq1, freq2, dm):
   return int(delay2(freq1, freq2, dm)/tInt + 0.5)
# end computeFrameIndexDelay()

def computeDelayedFrameIndex(startIndex, tInt, freq1, freq2, dm):
   return startIndex + computeFrameIndexDelay(tInt, freq1, freq2, dm)
# end computeDelayedFrameIndex()


def main_routine():
   comm  = MPI.COMM_WORLD
   rank  = comm.Get_rank()
   nProcs = comm.Get_size()

   maxpw = 600 #Maximum pulse width to search in seconds. default = 1 s.
   thresh= 5.0 #SNR cut off


   # Setup command-line parsing.
   usage="USAGE: %prog [options] <radio data filepath>" 
   cmdlnParser = OptionParser(usage=usage)
   cmdlnParser.add_option("-l", "--lower-FFT", dest="fcl", default=0, type="int", action="store",
                          help="Lower FFT index (between 0 and 4095) for the bandpass filtering.", 
                          metavar="INDEX")
   cmdlnParser.add_option("-u", "--upper-FFT", dest="fch", default=4095, type="int", action="store",
                          help="Upper FFT index (between 0 and 4095) for the bandpass filtering.", 
                          metavar="INDEX")
   cmdlnParser.add_option("-t", "--tuning", dest="tuning", default=0, type="int", action="store",
                          help="Tuning selection, 0 for low tuning, 1 for high tuning.",
                          metavar="TUNING")
   cmdlnParser.add_option("-s","--dm-start", dest="DMStart", default=0, type="float", action="store",
                          help="Starting DM (between 0 and 5000) for DM trials",
                          metavar="DM")
   cmdlnParser.add_option("-e","--dm-end", dest="DMEnd", default=5000, type="float", action="store",
                          help="Ending DM (between 0 and 5000) for DM trials",
                          metavar="DM")
   cmdlnParser.add_option("-f", "--frequency-file", dest="freqFilepath", default="freq1.npy", 
                          type="string", action="store",
                          help="Path to the bandpass frequency file.", metavar="PATH")
   cmdlnParser.add_option("-i", "--integration-time", dest="spectIntegTime", default=0.001,
                          type="float", action="store",
                          help="Spectral integration time in seconds.", metavar="SECS")
   cmdlnParser.add_option("-m", "--memory-limit", dest="memLimit", default=16384, type="int",
                          action="store",
                          help="Memory limit in MB for loading spectral files.  Default is 16384 MB.",
                          metavar="MB")
   cmdlnParser.add_option("-p", "--samples-per-sec", dest="nSamplesSec", default=0, type="int",
                          action="store",
                          help="Number of samples per second extracted from the data file.  " +
                          "This will override the total number of samples extracted.",
                          metavar="NUM")
   cmdlnParser.add_option("-n", "--samples", dest="nSamples", default=10000, type="int",
                          action="store",
                          help="Total number of samples extracted from the data file.",
                          metavar="NUM")
   cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                          action="store",
                          help="Working directory path.", metavar="PATH")

   # Parse the command-line and extract parameters.
   (cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()
   workDir = cmdlnOpts.workDir
   filepath = cmdlnParams[0]
   fcl = clipValue(cmdlnOpts.fcl, 0, 4095, valueType=int)
   fch = clipValue(cmdlnOpts.fch, 0, 4095, valueType=int)
   pol = clipValue(cmdlnOpts.tuning, 0, 1, valueType=int)
   DMstart = clipValue(cmdlnOpts.DMStart, 0.0, 5000.0, valueType=float)
   DMend = clipValue(cmdlnOpts.DMEnd, 0.0, 5000.0, valueType=float)
   freqFilepath = cmdlnOpts.freqFilepath
   memLimit = cmdlnOpts.memLimit
   spectIntegTime = cmdlnOpts.spectIntegTime
   nSamples = max(abs(cmdlnOpts.nSamples), nProcs)
   nSamplesSec = min(abs(cmdlnOpts.nSamplesSec), 1/spectIntegTime)
   nSamplesProc = nSamples/nProcs

   # Validate command-line inputs.
   if len(filepath) == 0:
      print('Path to the original data file must be provided')
      exit(1)
   # endif
   if fch <= fcl:
      print('Lower FFT index must be less than upper FFT index')
      exit(1)
   # endif
   if DMstart >= DMend:
      print('Starting DM trial value must be less than ending DM trial value.')
      exit(1)
   # endif

   # Extract the basename of the radio data file and use that as a pattern for loading the waterfall
   # files.
   filename = os.path.basename(os.path.splitext(filepath)[0])
   fn = sorted( glob.glob('{dir}/master{tune}_{filename}*.npy'.format(dir=workDir, tune=pol, 
                           filename=filename)) )

   # Load the time integration and frequency bandpass information.
   tInt = np.load('{dir}/tInt.npy'.format(dir=workDir))
   fullFreq=np.load(freqFilepath)
   fullFreq /= 1e6
   cent_freq = np.median(fullFreq)
   freq = fullFreq[fcl:(fch + 1)]
   freqDelta = freq[1] - freq[0]
   freqHigh = fullFreq[-1]
   numFreqs = len(freq)
   freqIndices = np.arange(numFreqs)
   BW = freq.max()-freq.min()

   npws = int(np.round(np.log2(maxpw/tInt)))+1 # +1 Due to in range(y) it goes to y-1 only


   try:
      nBytesFile = os.path.getsize(filepath)  # Radio data file size in bytes.
   except:
      print filepath,'not found'
      sys.exit(1)
   # end try
   #
   nFramesBeam = int(4)                            # Number of frames per beam data sample.
   nFramesFile = int(nBytesFile/drx.FrameSize)     # Number of frames in the radio data file.
   nFFTsTotal = int(nFramesFile/nFramesBeam)       # Number of FFTs in the total spectrogram.

   spectFileSize = os.path.getsize(fn[0])
   spect = np.load(fn[0])
   nFFTsPerFile = spect.shape[1]                   # Number of FFTs in a spectrogram sample file.
   spectSample = np.zeros((numFreqs, spect.shape[1]),
                           dtype=np.float32)       # Preallocated spectrogram sample.

   secsChunk = tInt   # Number of seconds span by a single beam data chunk (4 frames).
   nChunksSample = int(spectIntegTime/secsChunk + 0.5) # Number of beam data chunks per spectral sample.
   nFramesSample = nFramesBeam*nChunksSample       # Frames per spectral sample.
   nSamplesFile = int(nFramesFile/nFramesSample)   # Number of spectral samples in the data file.

   # If the user specified a number of samples per second of data for the waterfall, then compute 
   # the number of samples that need to be generated.  This is overriding any explicit setting of 
   # number of samples.
   if nSamplesSec > 0:
      nSecsFile = nFramesFile*secsChunk/nFramesBeam
      # Determine the number of samples to extract given the rate of samples per second of data.
      # Make sure this is not more than the number of samples that can extracted given the spectral
      # integration time.
      nSamples = int(min(nSecsFile*nSamplesSec, nSecsFile/spectIntegTime))
   # endif
   nSamples = min(nSamples, nSamplesFile)

   # The following values are used to control how we step through the radio data file.  We are not
   # trying to read every frame in the file.  Only the frames needed to build the coarse waterfall.
   #
   nSamplesProc = int(nSamples/nProcs)          # Number of samples each MPI process extracts.

   nDecimation = int(max(nSamplesFile/nSamples, 1 ))    # Decimation.


   # Compute the initial range of file indices that need to be used for the first decimation region of
   # the de-dispersed time-series.
   # CCY - TODO: Pull this segment of code into a separate utility in the workflow to allow setting the
   # DM range from within the workflow.
#  nFramesDecimation = nDecimation*nFFTsPerFile 
#  DMtrials = np.array([10, 20, 50, 100, 200, 500, 1000, 1500, 1600, 1700, 1800, 1900, 2000, 
#                       3000, 4000, 5000], dtype=np.float)
#  for DM in DMtrials:
#     startFileIndex = computeDelayedFrameIndex(0, tInt, freq[-1], freqHigh, DM)/nFramesDecimation
#     endFileIndex = computeDelayedFrameIndex(0, tInt, freq[0], freqHigh, DM)/nFramesDecimation
#     endFileIndex = min(endFileIndex, len(fn))
#     memoryRequired = (endFileIndex - startFileIndex + 1)*spectFileSize/float(1024**3)

#     print 'For DM =', DM
#     print '    Start file Index = ', startFileIndex
#     print '    End file Index = ', endFileIndex
#     print '    Memory required = {mem:.3f}'.format(mem=memoryRequired), 'GB'

#  exit(1)


   # Compute the maximum range of files that can be loaded given the memory constraints.
   nFilesPerSpan = int(memLimit*(1024**2)/(spectFileSize/2.0)) # Files are float64, but we are using
                                                               # float32 arrays.
   numFiles = len(fn)
   
   # Compute the maximum dispersion measure supported by the maximum file load.
   invScaledDelay = 1.0/scaledDelay(freq[0], freq[-1])
   DM_max = ((nFilesPerSpan*nDecimation + 1)*nFFTsPerFile - 1)*tInt*invScaledDelay
   # Warn user of the maximum and request if computation should proceed, regardless.
   procMessage('WARNING: Maximum supported dispersion measure ' + \
               'given memory constraints = {dm:.3f} pc cm^-3'.format(dm=DM_max), root=0)
   procMessage('DM trials start and end will be truncated to this maximum, if necessary.', root=0)
   procMessage('Current requested DM range is {start} to {end}.'.format(start=DMstart, end=DMend), root=0)
   procMessage('Proceed anyway [y/n]?', root=0)
   userAnswer = None
   if rank == 0:
      userAnswer = raw_input()
   # endif
   userAnswer = comm.bcast(userAnswer, root=0)
   if userAnswer.lower() not in ['y', 'yes', 'yeah', 'yea']:
      procMessage('De-dispersion search cancelled.', root=0)
      exit(1)
   # endif

   # Clip the start and end dispersion measures to the maximum.  We simply can't do better than this due
   # to memory limitations.
   DMstart = min(DMstart, DM_max)
   DMend = min(DMend, DM_max)
   # Build the array of dispersion-measure trials.
   DM = DMstart
   procMessage('Generating DM trial values'.format(rank=rank), root=0)
   numLowerTrials = 0
   numUpperTrials = 0
   numTotalTrials = 0
   if DMend > 1000:
      if DMstart < 1000:
         numLowerTrials = int((1000 - DMstart)*10.0)
         numUpperTrials = int(DMend - (numLowerTrials*0.1 + DMstart))
      else:
         numUpperTrials = int(DMend - DMstart)
      # endif
   else:
      numLowerTrials = int((DMend - DMstart)*10)
   # endif
   numTotalTrials = numLowerTrials + numUpperTrials
   DMtrials = np.arange(numTotalTrials, dtype=np.float64)
   DMtrials[:numLowerTrials] = DMtrials[:numLowerTrials]*0.1 + DMstart
   # DMtrials[numLowerTrials:numTotalTrials] += numLowerTrials*0.1 + DMstart - numLowerTrials
   DMtrials[numLowerTrials:numTotalTrials] += DMstart - 0.9*numLowerTrials
   # Get the DM trials to be executed by this process.
   nDMTrialsPerProc = int(numTotalTrials/nProcs)
   DMtrialStartIndex = rank*nDMTrialsPerProc
   DMtrialEndIndex = DMtrialStartIndex + nDMTrialsPerProc
   procDMTrials = DMtrials[ DMtrialStartIndex : DMtrialEndIndex ]
   # Rank 0 process will deal with the remainder at the end.
   if rank == 0 and nProcs*nDMTrialsPerProc < numTotalTrials:
      remainStartIndex = nProcs*nDMTrialsPerProc
      procDMTrials = np.concatenate((procDMTrials, DMtrials[remainStartIndex : numTotalTrials]),
                                    axis=0)
   # endif

   # Compute the largest segment of the de-dispersed time series we will see and allocate it.
   tsSize = ((nFilesPerSpan - 1)*nDecimation + 1)*nFFTsPerFile - \
               computeFrameIndexDelay(tInt, freq[-1], freqHigh, DMstart)
   ts = np.zeros(tsSize, dtype=np.float32)
   snr = np.zeros(tsSize, dtype=np.float32)
   #tsIndices = np.arange(tsSize, dtype=np.int64)

   # Open a single output file that is shared amongst all the processes. 
   filename = "{dir}/ppc_SNR_pol_{tune}.txt".format(dir=workDir, tune=pol)
   outfile = MPI.File.Open(comm, filename, MPI.MODE_CREATE | MPI.MODE_WRONLY)
   # Preallocate the pulse object that will be written to file.
   pulse = OutputSource()
   pulseID = 0.0


   fileIndex = 0
   nFilesLoaded = 0
   fileRingStart = 0
   fileRing = None
   datumSize = np.dtype(np.float32).itemsize
   # Prepare memory in process 0 to hold the loaded spectral files.  This memory will be shared to the
   # other processes.
   membuff = None
   if rank == 0:
      # Allocate the memory buffer and reshape it into the needed file ring buffer.
      memsize = nFilesPerSpan*nFFTsPerFile*numFreqs*datumSize
      membuff = MPI.Alloc_mem(memsize)
      fileRing = np.array(membuff, copy=False).view(dtype=np.float32).reshape(nFilesPerSpan, 
                           numFreqs, nFFTsPerFile)
      # Load initial files into the file ring.
      procMessage('Loading initial set of spectral sample files into memory shared file ring.', root=0)
      while fileIndex < numFiles and nFilesLoaded < nFilesPerSpan:
         ringIndex = (fileIndex + fileRingStart) % nFilesPerSpan
         fileRing[ringIndex, : , : ] = np.load(fn[fileIndex])
         fileIndex = fileIndex + 1
         nFilesLoaded = nFilesLoaded + 1
      # end while
      procMessage('{num} spectral sample files loaded.'.format(num=nFilesLoaded), root=0)
      fileIndex = fileIndex - 1
   # endif
   rmaWin = MPI.Win.Create(membuff, datumSize, comm=comm)
   rmaWin.Fence()
   # Broadcast the file ring information to all processes and reconstruct the ring in the current
   # process.
   fileIndex = comm.bcast(fileIndex, root=0)
   fileRingStart = comm.bcast(fileRingStart, root=0)
   nFilesLoaded = comm.bcast(nFilesLoaded, root=0)

   # == Precalculation and preallocation ==
   procMessage('Performing preallocation and precalculation.', root=0)
   fftIndices = np.arange(nFFTsPerFile, dtype=np.int32)
   freqIndices = np.arange(numFreqs, dtype=np.int32)
   # Preallocate the spectral sample cache.  Preallocate the spectral sample indices.  Preallocate the
   # DM indices.
   spectSample = np.zeros((numFreqs, nFFTsPerFile), dtype=np.float32)
   spectIndices = np.arange(nFilesLoaded, dtype=np.int32)
   dmIndices = np.arange(len(procDMTrials), dtype=np.int32)

   # Preallocate the indexed-time delays for spectral samples at all relevant DMs and frequencies.
   scaleDelay = scaledDelay(freq, freq[-1])/tInt
   spectScaleDelay = scaledDelay(freq[-1], freqHigh)/tInt
   freqIndexDelay = np.array( map(lambda x: scaleDelay*x + 0.5, procDMTrials), dtype=np.int32, copy=False )
   freqDedispIndices = np.array( map(lambda x: freqIndexDelay[x] - freqIndexDelay[x][0], dmIndices),
                                 dtype=np.int32, copy=False)
   spectIndexDelay = np.array( spectScaleDelay*procDMTrials + 0.5, dtype=np.int32, copy=False )
   # Preallocate de-dispersion arrays for de-dispersing a spectral sample.
   dedispSize = np.array( map(lambda x: freqIndexDelay[x][0], dmIndices) ) + nFFTsPerFile
   #dedispPower = np.zeros(dedispSize[-1], dtype=np.float32)
   csrIndptr = np.arange(numFreqs + 1, dtype=np.int32)*nFFTsPerFile
   csrIndices_init = np.tile(fftIndices, numFreqs)
   csrIndices_dedisp = np.array( [ [csrIndices_init[x*nFFTsPerFile : (x + 1)*nFFTsPerFile] + \
                                       freqDedispIndices[y][x] for x in freqIndices] \
                                     for y in dmIndices ], dtype=np.int32, copy=False).reshape(
                                     (dmIndices[-1] + 1, numFreqs*nFFTsPerFile) )
   csrSpectSample = spsparse.csr_matrix((spectSample.flat, csrIndices_init, csrIndptr), 
                                          shape=(numFreqs, dedispSize[-1]), dtype=np.float32, copy=False)
   # Preallocate boundary indices for merging de-dispersed samples into the de-dispersed time-series.
   offsetIndices = spectIndices*nDecimation*nFFTsPerFile
   endIndices = np.array( map(lambda x: offsetIndices - spectIndexDelay[x], dmIndices), 
                           dtype=np.int32, copy=False) + nFFTsPerFile
   beginIndices = np.array( map(lambda x: endIndices[x][:] - freqIndexDelay[x][0] - nFFTsPerFile, 
                                 dmIndices), 
                              dtype=np.int32, copy=False)
   cutIndices = np.maximum(-beginIndices, 0)
   beginIndices = np.maximum(beginIndices, 0)
   endIndices = np.maximum(endIndices, 0)
   blockIndex = 0
   nFilesPerSpanArray = np.full(shape=spectIndices.shape, fill_value=nFilesPerSpan, dtype=np.int32)
   nDoublesPerFile = nFFTsPerFile*numFreqs
   # Loop over spectral sample files.
   while fileIndex < numFiles:
      ringIndex = np.mod(spectIndices + fileRingStart, nFilesPerSpanArray)
      # Loop over DM trials.
      for dmIndex in dmIndices:
         DM = procDMTrials[dmIndex]
         procMessage('De-dispersing spectral block {block} with DM = {dm:.3f}'.format(block=blockIndex, 
                     dm=DM))
         # Loop over loaded spectral samples.
         for spectIndex in spectIndices:
            # Copy the current spectral sample from the memory buffer of loaded samples.
            # CCY - NOTE: because we created csrSpectSample with the option copy=False, its data array
            # points to the same memory region as spectSample.  So, changing the data values in
            # spectSample automatically changes the values in csrSpectSample.  Thus, we avoid a copy
            # step here.  Instead, we can just proceed to operate directly on csrSpectSample.
            rmaWin.Get(origin=(spectSample, spectSample.size, MPI.DOUBLE), 
                       target_rank=0, 
                       target=(ringIndex[spectIndex]*nDoublesPerFile, spectSample.size, 
                                 MPI.DOUBLE))
            # De-disperse the current spectral sample and add its de-dispersed power into the
            # de-dispersed time-series.
            if endIndices[dmIndex, spectIndex] > 0:
               cutIndex = cutIndices[dmIndex, spectIndex]
               sizeIndex = dedispSize[dmIndex]
               csrSpectSample.indices = csrIndices_dedisp[dmIndex]
               dedispPower = np.array(csrSpectSample.sum(axis=0), copy=False).reshape((dedispSize[-1],))
               np.add(ts[ beginIndices[dmIndex, spectIndex] : endIndices[dmIndex, spectIndex] ], 
                      dedispPower[cutIndex : sizeIndex],
                      out=ts[ beginIndices[dmIndex, spectIndex] : endIndices[dmIndex, spectIndex] ])
            # endif
         # end for spectIndex in spectIndices:

         procMessage('Thresholding de-dispersed SNR power and searching ' + \
                     'for pulses at DM = {dm:.3f}'.format(dm=DM))
         (snr[:], mean, std, pulseTimeIndices) = Threshold(ts, thresh, niter=0)

         # If there were pulses above threshold, record them.
         if len(pulseTimeIndices) > 0:
            fileStartIndex = (fileIndex + 1) - nFilesLoaded 
            timeIndexOffset = np.minimum(fileStartIndex*nDecimation*nFFTsPerFile - \
                                          indexDelay[dmIndex, -1], 0)
            procMessage('Found {num} signals above threshold for DM = {dm:.3f}'.format(dm=DM,
                        num=len(pulseTimeIndices)))
            procMessage('Outputting signals to results file.')
            for index in pulseTimeIndices:
               pulseID += 1.0
               pulse.pulse = pulseID + DM/10000.0
               pulse.SNR = snr[index]
               pulse.DM = DM
               pulse.time = (timeIndexOffset + index)*tInt
               pulse.dtau = tInt
               pulse.dnu = freqDelta
               pulse.nu = cent_freq
               pulse.mean = mean
               pulse.std = std
               outfile.Iwrite_shared( pulse.str() )
            # end for one in ones
            procMessage('Done outputting signals to results file.')
         else:
            procMessage('No signals found above threshold for DM = {dm:.3f}'.format(dm=DM))
         # endif

         # Reset the de-dispersed time-series.
         ts.fill(0.0)
         snr.fill(0.0)
      # end for DM in procTrials:
      
      rmaWin.Fence()
      #comm.Barrier()
      # Have process 0 load the next spectral sample file into the ring, if there are more such files.
      if rank == 0:
         fileIndex = fileIndex + 1
         if fileIndex < numFiles:
            procMessage('Loading file {index} into data ring.'.format(index=fileIndex), root=0)
            fileRingStart = (fileRingStart + 1) % nFilesPerSpan
            ringIndex = (nFilesLoaded + fileRingStart - 1) % nFilesPerSpan
            fileRing[ringIndex, : , : ] = np.load(fn[fileIndex])
         # endif
      # endif
      rmaWin.Fence()
      #comm.Barrier()
      fileIndex = comm.bcast(fileIndex, root=0)
      fileRingStart = comm.bcast(fileRingStart, root=0)
   # end while fileIndex < numFiles
   #
   rmaWin.Free()
   if rank == 0:
      MPI.Free_mem(membuff)
   # endif
   outfile.Close()
   procMessage('Done with de-dispersion.')
# end main



if __name__ == '__main__':
   main_routine()
   exit(0)
# endif
