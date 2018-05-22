import os
import sys
from optparse import OptionParser
import numpy
import drx
import errors
from apputils import forceIntValue

def main(args):
   LFFT = 4096 # Length of the FFT. 4096 is the size of a frame read.

   # Setup the command-line options.
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
   cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                          action="store",
                          help="Working directory path.", metavar="PATH")
   
   # Parse command-line for FFT indices and the radio data file path.
   (cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()
   Lfcl = forceIntValue(cmdlnOpts.Lfcl, 0, 4095)
   Lfch = forceIntValue(cmdlnOpts.Lfch, 0, 4095)
   Hfcl = forceIntValue(cmdlnOpts.Hfcl, 0, 4095)
   Hfch = forceIntValue(cmdlnOpts.Hfch, 0, 4095)
   inFilepath = cmdlnParams[0]
   workDir = cmdlnOpts.workDir

   # Validate command-line inputs.
   if Lfch <= Lfcl:
      print('ERROR: Low tuning lower FFT index must be less than the upper FFT index')
      exit(1)
   # endif
   if Hfch <= Hfcl:
      print('ERROR: High tuning lower FFT index must be less than the upper FFT index')
      exit(1)
   # endif
   if len(inFilepath) == 0:
      print('Path to the original data file must be provided')
      exit(1)
   # end if
   #

   # Open the radio data file.
   try:
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
         print 'zero division error computing sampling rate.'
         inFile.close()
         exit(1)
   except errors.syncError:
      print 'assuming the srate is 19.6 MHz'
      inFile.seek(-drx.FrameSize+1, 1)

   # Extract metadata for the radio data file using the first 4 frames.
   # CCY - NOTE: There is an implicit assumption that a given radio data file is associated with only a
   # single beam.
   #
   beam,tune,pol = junkFrame.parseID()
   beams = drx.getBeamCount(inFile)
   tunepols = drx.getFramesPerObs(inFile)
   tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
   beampols = tunepol

   # Use the first 4 frames to determing the high and low tuning frequencies
   #
   inFile.seek(-drx.FrameSize, 1)
   centralFreq1 = 0.0
   centralFreq2 = 0.0
   for i in xrange(4):
      junkFrame = drx.readFrame(inFile)
      b,t,p = junkFrame.parseID()
      if p == 0:
         if t == 0:
            try:
               centralFreq1 = junkFrame.getCentralFreq()
            except AttributeError:
               from dp import fS
               centralFreq1 = fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
         else:
            try:
               centralFreq2 = junkFrame.getCentralFreq()
            except AttributeError:
               from dp import fS
               centralFreq2 = fS * ((junkFrame.data.flags[0]>>32) & (2**32-1)) / 2**32
         # endif
      # endif
   # end for i in xrange(4)

   # Determine the set of frequencies on the FFT for both tunings and the temporal resolution.  Save
   # this information to disk.
   freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d = 1.0/srate))
   tInt = 1.0*LFFT/srate
   print 'Temporal resl = {time} secs'.format(time=tInt)
   print 'Channel width = {freq} Hz'.format(freq=1.0/tInt)
   freq1 = freq+centralFreq1
   freq2 = freq+centralFreq2
   print 'Low freq bandpass = {low} - {high} Hz at tuning {tuning} Hz'.format(low=freq1[Lfcl],
         high=freq1[Lfch],tuning=centralFreq1)
   print 'High freq bandpass = {low} - {high} Hz at tuning {tuning} Hz'.format(low=freq2[Hfcl],
         high=freq2[Hfch],tuning=centralFreq2)
   numpy.save('{dir}/tInt'.format(dir=workDir), tInt)
   numpy.save('{dir}/lowtunefreq'.format(dir=workDir), freq1)
   numpy.save('{dir}/hightunefreq'.format(dir=workDir), freq2)

   inFile.close()
# end main()
#
if __name__ == "__main__":
        main(sys.argv[1:])
