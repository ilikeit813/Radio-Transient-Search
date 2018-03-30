#from mpi4py import MPI
import mpi4py 
import mpi4py.MPI as MPI
import os
import sys
import numpy
import getopt
import drx
import time
import matplotlib.pyplot as plt

def Decimate_ts(ts, ndown=2):
   """
   Takes a N dimension array and decimates it by a factore of ndown, default = 2. with axis = 0
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
# end def Decimate_ts()



def main(args):
   totalrank = 12
   comm  = MPI.COMM_WORLD
   rank  = comm.Get_rank()
   t0 = time.time()
   nChunks = 10000 #the temporal shape of a file.
   LFFT = 4096 # Length of the FFT. 4096 is the size of a frame read.
   nFramesAvg = 4*LFFT/4096 # the intergration time under LFFT, 4 = beampols = 2X + 2Y (high and low tunes)
   
   #for offset_i in range(4306, 4309):# one offset = nChunks*nFramesAvg skiped
   for offset_i in range(100, 1000 ): # one offset = nChunks*nFramesAvg skiped
      print 'Top-level loop: offset_i = ', offset_i
      offset_i = 1.*totalrank*offset_i + rank
      offset = nChunks*nFramesAvg*offset_i
      # Build the DRX file
      try:
         # CCY - the original code required the original data file to be in the same directory from
         # which waterfall.py is run.  Unfortunately, this is rather unsafe to disk space and directory
         # integrity as it required either copying (extremely large) data files into the run directory
         # or running in the original data directory, dumping files that probably shouldn't be there.
         # The change here is to extract the base filename (no extension or directory path) from the
         # original input filepath---the base filename is used later for constructing the filename
         # for the output file.  This allows specifying a full path to the data file while dumping
         # resultant files in the current working directory.
         inFilepath = getopt.getopt(args,':')[1][0]
         if len(inFilepath) == 0:
            print('Path to the original data file must be provided')
            exit(1)
         # end if
         inFilename = os.path.basename(os.path.splitext(inFilepath)[0])
         inFile = open(inFilepath, "rb")
         nFramesFile = os.path.getsize(inFilepath) / drx.FrameSize #drx.FrameSize = 4128
      except:
         print inFilepath,' not found'
         sys.exit(1)
      # end try
      try:
         junkFrame = drx.readFrame(inFile)
         try:
            srate = junkFrame.getSampleRate()
            pass
         except ZeroDivisionError:
            print 'zero division error'
            break
         # end try
      except errors.syncError:
         print 'assuming the srate is 19.6 MHz'
         inFile.seek(-drx.FrameSize+1, 1)
      # end try
      inFile.seek(-drx.FrameSize, 1)
      beam,tune,pol = junkFrame.parseID()
      beams = drx.getBeamCount(inFile)
      tunepols = drx.getFramesPerObs(inFile)
      tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
      beampols = tunepol
      if offset != 0:
         inFile.seek(offset*drx.FrameSize, 1)
      if nChunks == 0:
         nChunks = 1
      nFrames = nFramesAvg*nChunks
      centralFreq1 = 0.0
      centralFreq2 = 0.0
      for i in xrange(4):
         junkFrame = drx.readFrame(inFile)
         b,t,p = junkFrame.parseID()
         if p == 0 and t == 0:
            try:
               centralFreq1 = junkFrame.getCentralFreq()
            except AttributeError:
               from dp import fS
               centralFreq1 = fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
         elif p == 0 and t == 2:
            try:
               centralFreq2 = junkFrame.getCentralFreq()
            except AttributeError:
               from dp import fS
               centralFreq2 = fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
         else:
            pass
      # end for i in xrange(4)
      #
      inFile.seek(-4*drx.FrameSize, 1)
      # Sanity check
      if nFrames > (nFramesFile - offset):
         raise RuntimeError("Requested integration time + offset is greater than file length")
      # Master loop over all of the file chunks
      freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
      tInt = 1.0*LFFT/srate
      print 'offset_i = ', offset_i
      print 'offset = ', offset
      print 'Temporal resl = ',tInt
      print 'Channel width = ',1./tInt
      freq1 = freq+centralFreq1
      freq2 = freq+centralFreq2
      #print tInt,freq1.mean(),freq2.mean()
      masterSpectra = numpy.zeros((nChunks, 2, LFFT-1))
      for i in xrange(nChunks):
         # Find out how many frames remain in the file.  If this number is larger
         # than the maximum of frames we can work with at a time (nFramesAvg),
         # only deal with that chunk
         framesRemaining = nFrames - i*nFramesAvg
         if framesRemaining > nFramesAvg:
            framesWork = nFramesAvg
         else:
            framesWork = framesRemaining
         #if framesRemaining%(nFrames/10)==0:
         #  print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
         count = {0:0, 1:0, 2:0, 3:0}
         data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
         # If there are fewer frames than we need to fill an FFT, skip this chunk
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
               break
            except errors.syncError:
               print "Sync Error"
               continue
            beam,tune,pol = cFrame.parseID()
            if tune == 0:
               tune += 1
            aStand = 2*(tune-1) + pol
            data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
            count[aStand] +=  1
         # end for j in xrange(framesWork)
         #
         # Calculate the spectra for this block of data
         masterSpectra[i,0,:] = ((numpy.fft.fftshift(numpy.abs(numpy.fft.fft2(data[:2,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
         masterSpectra[i,1,:] = ((numpy.fft.fftshift(numpy.abs(numpy.fft.fft2(data[2:,:]))[:,1:]))**2.).mean(0)/LFFT/2. #in unit of energy
         # Save the results to the various master arrays
         #print masterSpectra.shape
         #numpy.save('data',data)
         #sys.exit()
         #if i % 100 ==1 :
         #  print i, ' / ', nChunks
      # end for i in xrange(nChunks)
      #
      outname = "%s_%i_fft_offset_%.9i_frames" % (inFilename, beam, offset)
      numpy.save('waterfall' + outname, masterSpectra.mean(0) )
   #print time.time()-t0
   #print masterSpectra.shape
   #print masterSpectra.shape
# end def main()

# Main routine
if __name__ == "__main__":
   main(sys.argv[1:])
