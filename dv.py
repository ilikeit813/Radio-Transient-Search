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
from apputils import forceIntValue
from apputils import clipValue


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
   # Dispersion constant in MHz^2 s / pc cm^-3
   _D = 4.148808e3
   # Delay in s
   tDelay = dm*_D*((1.0/(freq1*freq1)) - (1.0/(freq2*freq2)))

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
   clip   -  Clipping SNR threshold for values to leave out of a mean/rms calculation.  default = 3.
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
     DM   = None  # DM (pc/cm3) of pulse
     time  = None  # Time at which pulse ocurred
     dtau  = None  # Temporal resolution of time series
     dnu   = None  # Spectral resolution
     nu   = None  # Central Observing Frequency
     mean  = None  # Mean in the time series
     rms   = None  # RMS in the time series

     formatter = "{0.pulse:07d}   {0.SNR:10.6f}    {0.DM:10.4f}    {0.time:10.6f} "+\
             "    {0.dtau:10.6f}    {0.dnu:.4f}    {0.nu:.4f}   {0.mean:.5f}"+\
             "   {0.rms:0.5f}\n " 

     def __str__(self):
        return self.formatter.format(self)




if __name__ == '__main__':
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
   cmdlnParser.add_option("-m", "--memory-limit", dest="memLimit", default=10, type="float",
                          action="store",
                          help="Memory limit per MPI process in MB.  Default is 10 MB.",
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

   # Parse the command-line and extract parameters.
   (cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()
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
   nSamplesSec = min(abs(cmdlnOpts.nSamplesSec), 1/spectIntTime)
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
   fn = sorted( glob.glob('master{tune}{filename}*.npy'.format(tune=pol, filename=filename)) )

   # Load the time integration and frequency bandpass information.
   tInt = np.load('tInt.npy')
   fullFreq=np.load(freqFilepath)
   fullFreq /= 1e6
   cent_freq = np.median(fullFreq)
   freq = fullFreq[fcl:fch]
   numFreqs = len(freq)
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
   nFFTsPerFile = spect.shape[0]                   # Number of FFTs in a spectrogram sample file.

   secsChunk = tInt   # Number of seconds span by a single beam data chunk (4 frames).
   nChunksSample = int(spectIntTime/secsChunk + 0.5) # Number of beam data chunks per spectral sample.
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
      nSamples = int(min(nSecsFile*nSamplesSec, nSecsFile/spectIntTime))
   # endif
   nSamples = min(nSamples, nSamplesFile)

   # The following values are used to control how we step through the radio data file.  We are not
   # trying to read every frame in the file.  Only the frames needed to build the coarse waterfall.
   #
   nSamplesProc = int(nSamples/nProcs)          # Number of samples each MPI process extracts.

   nDecimation = max(int(nSamplesFile/nSamples), int(1))    # Decimation.


   DM = DMstart
   DMtrials = None
   if rank == 0:
      print 'Process {rank}: Generating DM trial values'.format(rank=rank)
      numLowerTrials = 0
      numUpperTrials = 0
      numTotalTrials = 0
      if DMend > 1000:
         if DMstart < 1000:
            numLowerTrials = int((1000 - DMstart)/0.1)
            numUpperTrials = int(DMend - (numLowerTrials*0.1 + DMstart))
         else:
            numUpperTrials = int(DMend - DMstart)
         # endif
      else:
         numLowerTrials = int(DMend - DMstart)
      # endif
      numTotalTrials = numLowerTrials + numUpperTrials
      DMtrials = np.arange(numTotalTrials, dtype=np.float)
      DMtrials[:numLowerTrials] = DMtrials[:numLowerTrials]*0.1 + DMstart
      # DMtrials[numLowerTrials:numTotalTrials] += numLowerTrials*0.1 + DMstart - numLowerTrials
      DMtrials[numLowerTrials:numTotalTrials] += DMstart - 0.9*numLowerTrials
   # endif
   DMtrials = comm.bcast(DMtrials,root =0)

   numFiles = len(fn)
   filesPerProc = int(32*memLimit*tInt/spectIntegTime + 0.5)
   filesPerBlock = nProcs*filesPerProc
   spect = np.load(fn[0]) # Load sample file to get number of frames per file.
   spectSample = np.zeros((spect.shape[0], numFreqs)) # Preallocated spectrogram sample.

   # Open a single output file that is shared amongst all the processes. 
   filename = "ppc_SNR_pol_%.1i.txt" % (pol)
   outfile = MPI.File.Open(comm, filename, MPI.MODE_CREATE | MPI.MODE_WRONLY)

   # Preallocate the pulse object that will be written to file.
   pulse = OutputSource()

   # Preallocate the de-dispersed time-series.
   ts = np.zeros(nFFTsTotal)  # De-dispersed time-series built only from the spectrogram samples
                              # assigned to this process.
   tsMerge = None             # Merged de-dispersed time-series from all processes.

   # Preallocate the lag-free time-series segment that is searched for pulses.
   endLagIndex = len(tsMerge) - tbMax
   searchSize = endLagIndex - tbMax
   tsLagFree = np.zeros(searchSize)
   sn = np.zeros(searchSize)


   txtsize=np.zeros((npws,2),dtype=np.int32) # fileno = txtsize[ranki,0], pulse number = txtsize[ranki,1],
                                             # ranki is the decimated order of 2
   txtsize[:,0]=1 #fileno star from 1
   # Loop over each trial DM.
   for DM in DMtrials:
      print 'Process {rank}: Performing trial DM = {dm}'.format(rank=rank, dm=DM)

      if rank == 0:
         print 'De-dispersing spectrogram samples'

      # De-disperse the spectrogram samples being handled by this process into the de-dispersed
      # time-series.
      fileIndex = rank*filesPerProc
      while fileIndex < numFiles:
         nFileCount = 0
         while fileIndex < numFiles and nFileCount < filesPerProc:
            spectSample[:,:] = np.load(fn[fileIndex])

            for freqIndex in range(numFreqs):
               # Compute the de-dispersion shift, in terms of spectrogram FFTs, for a given frequency in
               # the bandpass.
               timeBin = np.int32(delay2(freq[freqbin], freq[-1], DM)/tInt + 0.5)
               # Compute where the de-dispersed spectrogram sample lies in the de-dispersed time-series.
               tsOffset = fileIndex*nDecimation*nFFTsPerFile - timeBin
               tsBegin = tsOffset + nFileCount*nDecimation*nFFTsPerFile
               tsEnd = beginIndex + FFTsPerFile
               # Add to the de-dispersed time-series only the portion of the de-dispersed spectrogram
               # sample that lies within the de-dispersed time-series.
               if endIndex >= 0:
                  beginCutIndex = 0
                  if beginIndex < 0:
                     beginCutIndex = -beginIndex
                     beginIndex = 0
                  # endif
                  ts[beginIndex : endIndex] += spectSample[beginCutIndex:,freqbin]
               # endif
            # end for freqIndex in range(numFreqs):
            
            nFileCount = nFileCount + 1
            fileIndex = fileIndex + 1
         # end while fileIndex < numFiles and nFileCount < filesPerProc

         #fileIndex += filesPerBlock
         fileIndex = fileIndex - filesPerProc + filesPerBlock
      # end while fileIndex < numFiles

      if rank == 0:
         print 'Merging time-series'

      # Merge the de-dispersed time-series from all processes.
      comm.Allreduce(ts, tsMerge, op=MPI.SUM)
      ts.fill(0.0)


      #"""#search for signal with decimated timeseries
      if rank < npws: # timeseries is ready for signal search
         if rank == 0:
            print 'Thresholding time-series to search for bursts.'

         ndown = 2**rank # decimate the time series
         sn[:],mean,rms = Threshold(Decimate_ts(tsMerge,ndown),thresh,niter=0)

         if rank == 0:
            print 'Outputing found signals'

         # Record all pulses above threshold in the section searched by this process.
         for one in np.where(sn != -1)[0]:
            txtsize[rank,1] += 1
            pulse.pulse = txtsize[rank,1]
            pulse.SNR = sn[one]
            pulse.DM = DM
            pulse.time = one*tInt*ndown
            pulse.dtau = tInt*ndown
            pulse.dnu = freq[1]-freq[0]
            pulse.nu = cent_freq
            pulse.mean = mean
            pulse.rms = rms
            outfile.Write_Ordered(pulse.formatter.format(pulse)[:-1]) 
            #if txtsize[ranki,1] >200000*txtsize[ranki,0]:
            #   outfile.close()
            #   txtsize[ranki,0]+=1
            #   filename = "ppc_SNR_pol_%.1i_td_%.2i_no_%.05d.txt" % (pol,ranki,txtsize[ranki,0])
            #   outfile = open(filename,'a')
            # endif
         # end for one in ones
         if rank == 0:
            print 'Proceeding to next file set'
      # endif
   # end for DM in DMtrials

   outfile.Close()
# end main
