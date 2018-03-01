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

def main(args):
   nodes = 4 #total blades used
   pps = 6   #process per blade

   windownumber = 4 # The length of FFT = windownumber * 4096

   #Low tuning frequency range
   Lfcl = 1700 * windownumber
   Lfch = 2100 * windownumber
   #High tuning frequency range
   Hfcl =  670 * windownumber
   Hfch = 1070 * windownumber

   totalrank = nodes*pps
        comm  = MPI.COMM_WORLD
        rank  = comm.Get_rank()
   t0 = time.time()
   nChunks = 3000 #the temporal shape of a file.
   LFFT = 4096 * windownumber #Length of the FFT. 4096 is the size of a frame readed. The mini quantized window lenght is 4096
   nFramesAvg = 1*4* windownumber # the intergration time under LFFT, 4 = beampols = 2X + 2Y (high and low tunes)
   
   #for offset_i in range(4306, 4309):# one offset = nChunks*nFramesAvg skiped
   for offset_i in range(0, 1000 ):# one offset = nChunks*nFramesAvg*worker_rank skiped
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
         inFilename = os.path.basename(os.path.splitext(inFilepath)[0])
         inFile = open(inFilepath, "rb")
         nFramesFile = os.path.getsize(inFilepath) / drx.FrameSize #drx.FrameSize = 4128
      except:
         print inFilepath,' not found'
         sys.exit(1)
      try:
         junkFrame = drx.readFrame(inFile)
         try:
            srate = junkFrame.getSampleRate()
            pass
         except ZeroDivisionError:
            print 'zero division error'
            break
      except errors.syncError:
         print 'assuming the srate is 19.6 MHz'
         inFile.seek(-drx.FrameSize+1, 1)
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
      inFile.seek(-4*drx.FrameSize, 1)
      # Sanity check
      if nFrames > (nFramesFile - offset):
         raise RuntimeError("Requested integration time + offset is greater than file length")
      # Master loop over all of the file chunks
      #freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
      #tInt = 1.0*LFFT/srate
                #print 'Temporal resl = ',tInt
                #print 'Channel width = ',1./tInt
      #freq1 = freq+centralFreq1
      #freq2 = freq+centralFreq2
      #print tInt,freq1.mean(),freq2.mean()
      masterSpectra = numpy.zeros((nChunks, 2, Lfch-Lfcl))
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
         # Calculate the spectra for this block of data, in the unit of intensity
         masterSpectra[i,0,:] = ((numpy.fft.fftshift(numpy.abs(numpy.fft.fft2(data[:2,:]))[:,1:])[:,Lfcl:Lfch])**2.).mean(0)/LFFT/2.
         masterSpectra[i,1,:] = ((numpy.fft.fftshift(numpy.abs(numpy.fft.fft2(data[2:,:]))[:,1:])[:,Hfcl:Hfch])**2.).mean(0)/LFFT/2.
      # Save the results to the various master arrays
      outname = "%s_%i_fft_offset_%.9i_frames" % (inFilename, beam,offset)
      numpy.save(outname,masterSpectra)
if __name__ == "__main__":
   main(sys.argv[1:])
