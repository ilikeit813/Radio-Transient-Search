#from mpi4py import MPI
import mpi4py 
import mpi4py.MPI as MPI
import os
import sys
import numpy as np
import getopt
import drx
import errors
import time
import math
import matplotlib.pyplot as plt
from optparse import OptionParser

from apputils import Decimate, procMessage, atProcMessage


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
   cmdlnParser.add_option("-s", "--samples", dest="nSamples", default=10000, type="int",
                           action="store",
                           help="Number of total spectral samples to produce for the waterfall.",
                           metavar="NUM")
   cmdlnParser.add_option("-p", "--samples-per-sec", dest="nSamplesSec", default=0, type="int",
                          action="store",
                          help="Number of spectral samples per second generated for the waterfall." +
                          "This will override the setting for the total number of samples",
                          metavar="NUM")
   cmdlnParser.add_option("-d", "--detailed", dest="fDetailed", default=False,
                          action="store_true",
                          help="Flag denoting to perform a detailed spectrogram.")
   cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                          action="store",
                          help="Working directory path.", metavar="PATH")
   (cmdlnOpts, args) = cmdlnParser.parse_args()
   if len(args) == 0:
      print "Must supply a path to the radio data file."
      sys.exit(1)
   # endif

   # Extract command-line options and arguments.
   inFilepath = args[0]
   inFilename = os.path.basename(os.path.splitext(inFilepath)[0])
   spectIntTime = cmdlnOpts.spectIntTime
   nSamples = max(abs(cmdlnOpts.nSamples), nProcs)
   nSamplesSec = min(abs(cmdlnOpts.nSamplesSec), 1/spectIntTime)
   nSamplesProc = nSamples/nProcs
   workDir = cmdlnOpts.workDir


   # Open the radio data file.
   try:
      inFile = open(inFilepath, "rb")
      nBytesFile = os.path.getsize(inFilepath)  # Radio data file size in bytes.
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

      secsChunk = 4096.0/srate   # Number of seconds span by a single beam data chunk (4 frames).
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

      nFileSamplesProc = nDecimation*nSamplesProc        # Number of samples in the file spanned by the
                                                         # region processed by a single MPI process.
      sampleFileOffset = nFramesSample*drx.FrameSize     # File offset of each spectral sample.
      procFileOffset = nFileSamplesProc*sampleFileOffset # File offset of each MPI process region.
      deciFileOffset = nDecimation*sampleFileOffset      # File offset of each decimated sample that we
                                                         # actually extract.


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
   # Allocate the output spectra depending whether we are doing a coarse or detailed spectrogram.
   if not cmdlnOpts.fDetailed:
      outSpectra = np.zeros((2, LFFT)) # Waterfall spectrum.
      spectraTag = 'waterfall'
   else:
      outSpectra = np.zeros((nChunksSample, 2, LFFT)) # Spectrogram of a single sample.
      spectraTag = 'master'
   # endif
   data = np.zeros((4,framesWork*4096/beampols), dtype=np.csingle)
   freq = np.fft.fftshift(np.fft.fftfreq(LFFT, d = 1.0/srate))
   tInt = 1.0*LFFT/srate
   channelWidth = 1e-6/tInt
   centralFreq1MHz = centralFreq1*1e-6
   centralFreq2MHz = centralFreq2*1e-6
   normFactor = 1.0/(4*LFFT)
   avgFactor = 1.0/nChunksSample


   # Determine which samples will be generated by this process.
   sampleIndices = np.arange(nSamplesProc) + nSamplesProc*procRank
   if procRank == 0:
      sampleIndices = np.concatenate( (sampleIndices, np.arange(nSamplesProc*nProcs, nSamples)), axis=0 )
   # endif
   chunkIndices = np.arange(nChunksSample)

   # Loop over all the samples to be generated by this process.
   inFile.seek(sampleIndices[0]*nDecimation*nFramesSample*drx.FrameSize, 0)
   for sampleIndex in sampleIndices:
      nStartFrame = sampleIndex*nDecimation*nFramesSample
      nEndFrame = nStartFrame + nFramesSample - 1
      procMessage(''.join(['Spectral integration of sample ', 
            '{sample} of {totsamples} '.format(sample=sampleIndex+1, totsamples=nSamples),  
            '(frames {beg} - {end} of {frames}):'.format(beg=nStartFrame, end=nEndFrame,  
                                                         frames=nFramesFile)]) )
      if procRank == 0:
         print '   Temporal resl = {value} secs'.format(value=tInt)
         print '   Channel width = {value:7.6g} MHz'.format(value=channelWidth)
         print '   Low Tuning = {value:6.5g} MHz'.format(value=centralFreq1MHz)
         print '   High Tuning = {value:6.5g} MHz'.format(value=centralFreq2MHz)

         # The rank 0 process is handling the remainder.  So, if it has done nSamplesProc samples,
         # already, then the next samples it will do is the remainder.  Seek to the start position of the
         # remainder.
         if sampleIndex > nSamplesProc:
            inFile.seek(sampleIndex*nDecimation*nFramesSample*drx.FrameSize, 0)
      # endif
      
      # Seek to the beginning of the sample.
      for i in chunkIndices:
         # If the number of remaining frames is less than a complete set for the number of beams, then
         # just work with however many frames remain.
         if framesRemaining < framesWork:
            framesWork = framesRemaining
            data = np.zeros((4,framesWork*4096/beampols), dtype=np.csingle)
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

         # Integrate the spectra for the current data chunk.
         #
         # Extract low-tuning FFTs.
         dataFFT = np.fft.fftshift( np.fft.fft(data[:2,:].mean(0))[:] )
         lowPowerFFT = normFactor*(dataFFT.real*dataFFT.real + dataFFT.imag*dataFFT.imag)
         # Extract high-tuning FFTs.
         dataFFT = np.fft.fftshift( np.fft.fft(data[2:,:].mean(0))[:] )
         highPowerFFT = normFactor*(dataFFT.real*dataFFT.real + dataFFT.imag*dataFFT.imag)
         if not cmdlnOpts.fDetailed:
            outSpectra[0,:] += lowPowerFFT*avgFactor
            outSpectra[1,:] += highPowerFFT*avgFactor
         else:
            outSpectra[i,0,:] = lowPowerFFT
            outSpectra[i,1,:] = highPowerFFT
         # endif

         #masterSpectra[i,0,:] = ((np.fft.fftshift(np.abs(
         #      np.fft.fft2(data[:2,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
         #masterSpectra[i,1,:] = ((np.fft.fftshift(np.abs(
         #      np.fft.fft2(data[2:,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
      # end for i in xrange(nChunksSample)

      # Build the DRX file
      outname = "%s_%i_fft_offset_%.9i_frames" % (inFilename, beam, nStartFrame)
      np.save('{dir}/{tag}{label}'.format(dir=workDir, tag=spectraTag, label=outname), 
                  outSpectra)

      # Move to the next sample.
      inFile.seek(deciFileOffset - sampleFileOffset, 1)
   # end while sampleIndex < nFileSamplesProc && sampleCount <= nSamples:

   inFile.close()
# end main()

# Main routine
if __name__ == "__main__":
   main(sys.argv[1:])
   sys.exit(0)
