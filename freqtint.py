import os
import sys
import numpy
import getopt
import drx
from apputils import forceIntValue

def main(args):
   windownumber = 4

   #Low tuning frequency range
   #Lfcl = 2700 * windownumber
   #Lfch = 2800 * windownumber
   Lfcl = 0
   Lfch = 4095 * windownumber
   #High tuning frequency range
   #Hfcl = 1500 * windownumber
   #Hfch = 1600 * windownumber
   Hfcl = 0
   Hfch = 4095 * windownumber

   nChunks = 3000 #the temporal shape of a file.
   LFFT = 4096 * windownumber #Length of the FFT.4096 is the size of a frame readed.
   nFramesAvg = 1 * 4 * windownumber # the intergration time under LFFT, 4 = beampols = 2X + 2Y (high and low tunes)

   #for offset_i in range(4306, 4309):# one offset = nChunks*nFramesAvg skiped
   for offset_i in range(0, 1):# one offset = nChunks*nFramesAvg skiped
      offset = 0
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
         #
         inFilepath = getopt.getopt(args,':')[1][0]
         if len(inFilepath) == 0:
            print('Path to the original data file must be provided')
            exit(1)
         # end if
         inFilename = os.path.basename(os.path.splitext(inFilepath)[0])
         inFile = open(inFilepath, "rb")
         nFramesFile = os.path.getsize(inFilepath) / drx.FrameSize #drx.FrameSize = 4128

         # Get the FFT indices for the low and high frequency tunings from the user.
         szLongOpts = ['low-tuning-lower', 'low-tuning-upper', 'high-tuning-lower', 'high-tuning-upper']
         szShortOpts = 'abyz'
         (cmdLnOpts, cmdLnParams) = getopt.getopt(sys.argv[1:], szShortOpts, szLongOpts)[0]
         for index in range(len(cmdLnOpts)):
             # Set the FFT index for the specified tuning and force it to be a non-negative integer less
             # than 4096.
             if cmdLnOpts[index] in [szShortOpts[0], szLongOpts[0]]:
                 Lfcl = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
             elif cmdLnOpts[index] in [szShortOpts[1], szLongOpts[1]]:
                 Lfch = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
             elif cmdLnOpts[index] in [szShortOpts[2], szLongOpts[2]]:
                 Hfcl = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
             elif cmdLnOpts[index] in [szShortOpts[3], szLongOpts[3]]:
                 Hfch = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
             else:
                 print('UNKNOWN OPTION: {opt}'.format(opt=cmdLnOpts[index]))
                 exit(1)
             # end if
         # end for index in range(len(cmdLnOpts))
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
      freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
      tInt = 1.0*LFFT/srate
      print 'Temporal resl = ',tInt
      print 'Channel width = ',1./tInt
      freq1 = freq+centralFreq1
      freq2 = freq+centralFreq2
      print 'Low  freq  = ', freq1[Lfcl],freq1[Lfch],' at', freq1[Lfcl]/2+freq1[Lfch]/2
      print 'High freq  = ', freq2[Hfcl],freq2[Hfch],' at', freq2[Hfcl]/2+freq2[Hfch]/2
      numpy.save('tInt',tInt)
      numpy.save('freq1',freq1[Lfcl:Lfch])
      numpy.save('freq2',freq2[Hfcl:Hfch])
if __name__ == "__main__":
        main(sys.argv[1:])
