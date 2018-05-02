#from mpi4py import MPI
import mpi4py 
import mpi4py.MPI as MPI
import os
import sys
import numpy
import getopt
import drx
import errors
import time
import math
import matplotlib.pyplot as plt
from optparse import OptionParser

from apputils import Decimate


def main(args):
   # Get multi-processing information.
   MPIComm = MPI.COMM_WORLD
   nProcs = MPIComm.Get_size()
   procRank = MPIComm.Get_rank()

   LFFT = 4096 # Length of the FFT. 4096 is the size of a frame read.

   # Initialize command-line parser and then parse the command-line
   usage = " Usage: %prog [options] <radio filepath>"
   cmdlnParser = OptionParser(usage=usage)
   cmdlnParser.add_option("-i", "--integrate-time", dest="spectIntTime", default=0.001, type="float", 
                           action="store", 
                           help="Spectral integration time in seconds.", metavar="SECS")
   cmdlnParser.add_option("-s", "--samples", dest="nSamples", default=4000, type="long",
                           action="store",
                           help="Number of spectral samples for coarse spectrogram.",
                           metavar="NUM")
   (cmdlnOpts, args) = cmdlnParser.parse_args()
   if len(args) == 0:
      print "Must supply a path to the radio data file."
      sys.exit(1)
   # endif

   # Extract command-line options and arguments.
   inFilepath = args[0]
   inFilename = os.path.basename(os.path.splitext(inFilepath)[0])
   spectIntTime = cmdlnOpts.spectIntTime
   nSamples = max(cmdlnOpts.nSamples, nProcs)
   nSamplesProc = nSamples/nProcs

   # Open the radio data file.
   try:
      inFile = open(inFilepath, "rb")
      nBytesFile = long(os.path.getsize(inFilepath))  # Radio data file size in bytes.
   except:
      print inFilepath,' not found'
      sys.exit(1)
   # end try
   #
   nFramesBeam = int(4)                            # Number of frames per beam data sample.
   nFramesFile = nBytesFile/int(drx.FrameSize)     # Number of frames in the radio data file.

   # Parse the first few frames to extract metadata associated with the radio data stream.
   try:
      junkFrame = drx.readFrame(inFile)
      # CCY - The information obtained from the 4 lines below is not used because the main loop that
      # reads the radio data file only assumes a single beam for the entire radio data file.  However,
      # we keep these lines, just in case.
      beams = drx.getBeamCount(inFile)
      tunepols = drx.getFramesPerObs(inFile)
      tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
      beampols = tunepol
      try:
         srate = junkFrame.getSampleRate()
         pass
      except ZeroDivisionError:
         print 'zero division error computing sample rate.'
         inFile.close()
         sys.exit(1)
      # end try

      # The following values are used to control how we step through the radio data file.  We are not
      # trying to read every frame in the file.  Only the frames needed to build the coarse waterfall.
      #
      nChunksSample = long(srate*spectIntTime)     # Number of beam data chunks per spectral sample.
      nFramesSample = nFramesBeam*nChunksSample    # Frames per spectral sample.
      nSamplesFile = nFramesFile/nFramesSample     # Number of spectral samples in the data file.
      nFileSamplesProc = nSamplesFile/nProcs       # Number of spectral samples in the file covered by 
                                                   # each MPI process.

      nDecimateFactor = max(nSamplesFile/nSamples, long(1))    # Decimation factor.

      sampleFileOffset = nFramesSample*drx.FrameSize  # File offset of each spectral sample.
      procFileOffset = nFileSamplesProc*sampleFileOffset  # File offset of each MPI process region.
      deciFileOffset = nDecimateFactor*sampleFileOffset  # File offset of each decimated sample that we
                                                         # actually process.


      # Scan the first 4 frames to determine the low and high tuning frequencies.
      # CCY - NOTE: there is an implicit assumption that a given radio data file contains only 1 beam.
      #
      inFile.seek(-drx.FrameSize, 1)
      centralFreq1 = 0.0
      centralFreq2 = 0.0
      for i in xrange(nFramesBeam):
         junkFrame = drx.readFrame(inFile)
         beam,tune,pol = junkFrame.parseID()
         if pol == 0:
            if tune == 0:
                  centralFreq1 = junkFrame.getCentralFreq()
            else:
                  centralFreq2 = junkFrame.getCentralFreq()
            # endif
         # endif
      # end for i in xrange(nFramesBeam)
   except errors.syncError:
      print 'Could not read radio data file for metadata.'
      inFile.close()
      sys.exit(1)
   # end try
   #
   # Jump to the point in the file that we need to read for the current MPI process.
   fileOffset = procRank*procFileOffset
   inFile.seek(fileOffset - nFramesBeam*drx.FrameSize, 1)

   
   framesRemaining = nFramesFile
   framesWork = nFramesBeam
   masterSpectra = numpy.zeros((2, LFFT-1))  # CCY - NOTE: The DC component is removed from the FFT
                                             # before saving the spectrogram.  Not sure why this is
                                             # done, but it essentially removes the component with the
                                             # central frequency.
   data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
   freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
   tInt = 1.0*LFFT/srate
   channelWidth = 1e-6/tInt
   centralFreq1MHz = centralFreq1*1e-6
   centralFreq2MHz = centralFreq2*1e-6
   normFactor = 2*LFFT*nChunksSample
   # Master loop over all of the spectral samples to be processed by the current MPI process.
   for nProcSpecial in range(2):
      nProcStartFrame = nFileSamplesProc*procRank*nFramesSample
      sampleIndex = nSamplesProc*procRank + 1
      # We want the rank 0 process, and only the rank 0 process, to handle the trailing extra samples 
      # remaining once it is done with the normal segment of samples.  This requires us to change the 
      # starting frame number and position in the file to just the samples that remain.  Everything else
      # proceeds in the same manner.
      if nProcSpecial == 1:
         if procRank == 0:
            samplesDone = nSamplesProc*nProcs
            if samplesDone < nSamples:
               sampleIndex = samplesDone
               nProcStartFrame = nFileSamplesProc*nProcs*nFramesSample
               # Determine how many samples remain in the file and jump to the start of that remainder.
               nFileSamplesProc = max(nSamplesFile - nProcs*nFileSamplesProc, 0)
               inFile.seek(-nFramesSample*nFileSamplesProc*drx.FrameSize, 2)
            else:
               # We have all the samples, so don't do the loop to process samples.
               nFileSamplesProc = 0
            # endif
         else:
            # This is not the rank 0 process, so don't process anymore samples.
            nFileSamplesProc = 0
         # endif
      # endif
      nCount = 0
      while nCount < nFileSamplesProc:
         nStartFrame = nProcStartFrame + nCount*nFramesSample
         nEndFrame = nStartFrame + nFramesSample - 1
         print 'Spectral integration of sample {sample} of {totsamples}'.format(sample=sampleIndex,
               totsamples=nSamples), \
               '(frames {beg} - {end} of {frames}):'.format(beg=nStartFrame, end=nEndFrame, 
                                                            frames=nFramesFile)
         print '   Temporal resl = {value} secs'.format(value=tInt)
         print '   Channel width = {value:7.6g} MHz'.format(value=channelWidth)
         print '   Low Tuning = {value:6.5g} MHz'.format(value=centralFreq1MHz)
         print '   High Tuning = {value:6.5g} MHz'.format(value=centralFreq2MHz)
         for i in xrange(nChunksSample):
            # If the number of remaining frames is less than a complete set for the number of beams, then
            # just work with however many frames remain.
            if framesRemaining < framesWork:
               framesWork = framesRemaining
               data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
            # endif

            count = {0:0, 1:0, 2:0, 3:0}
            # If there are fewer samples than we need to fill an FFT, skip this chunk
            if data.shape[1] < LFFT:
               print 'data.shape[1]< LFFT, break'
               break
            # Inner loop that actually reads the frames into the data array
            for j in xrange(framesWork):
               # Read in the next frame and anticipate any problems that could occur
               try:
                  cFrame = drx.readFrame(inFile, Verbose=False)
               except errors.eofError:
                  print "EOF Error"
                  inFile.close()
                  sys.exit(1)
               except errors.syncError:
                  print "Sync Error"
                  continue
               beam,tune,pol = cFrame.parseID()
               aStand = tune + pol
               data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
               count[aStand] +=  1
            # end for j in xrange(framesWork)

            # Calculate the spectra for the current data chunk.
            #
            dataFFT = numpy.fft.fftshift( numpy.fft.fft(data[:2,:].mean(0))[1:] )
            diffPower = (dataFFT.real*dataFFT.real + dataFFT.imag*dataFFT.imag)/(normFactor)
            masterSpectra[0,:] += diffPower

            dataFFT = numpy.fft.fftshift( numpy.fft.fft(data[2:,:].mean(0))[1:] )
            diffPower = (dataFFT.real*dataFFT.real + dataFFT.imag*dataFFT.imag)/(normFactor)
            masterSpectra[1,:] += diffPower


            #masterSpectra[i,0,:] = ((numpy.fft.fftshift(numpy.abs(
            #      numpy.fft.fft2(data[:2,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
            #masterSpectra[i,1,:] = ((numpy.fft.fftshift(numpy.abs(
            #      numpy.fft.fft2(data[2:,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
         # end for i in xrange(nChunksSample)

         # Build the DRX file
         outname = "%s_%i_fft_offset_%.9i_frames" % (inFilename, beam, nStartFrame)
         numpy.save('waterfall' + outname, masterSpectra )

         # Skip forward to the next decimated sample.  Then reset the spectrum to zero.
         inFile.seek(nDecimateFactor*nFramesSample*drx.FrameSize, 1)
         nCount = nCount + nDecimateFactor
         masterSpectra.fill(0.0)
         sampleIndex = sampleIndex + 1
      # end while nCount < nFileSamplesProc:

   # end while nProcSpecial < 2:

   inFile.close()
# end def main()

# Main routine
if __name__ == "__main__":
   main(sys.argv[1:])
   sys.exit(0)
