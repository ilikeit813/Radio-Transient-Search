import os
import sys
import numpy
import getopt
import drx

def main(args):
	windownumber = 4

	#Low tuning frequency range
	Lfcl = 2700 * windownumber
	Lfch = 2800 * windownumber
	#High tuning frequency range
	Hfcl = 1500 * windownumber
	Hfch = 1600 * windownumber

	nChunks = 3000 #the temporal shape of a file.
	LFFT = 4096 * windownumber #Length of the FFT.4096 is the size of a frame readed.
	nFramesAvg = 1 * 4 * windownumber # the intergration time under LFFT, 4 = beampols = 2X + 2Y (high and low tunes)

	#for offset_i in range(4306, 4309):# one offset = nChunks*nFramesAvg skiped
	for offset_i in range(0, 1):# one offset = nChunks*nFramesAvg skiped
		offset = 0
		# Build the DRX file
		try:
                        fh = open(getopt.getopt(args,':')[1][0], "rb")
                        nFramesFile = os.path.getsize(getopt.getopt(args,':')[1][0]) / drx.FrameSize #drx.FrameSize = 4128
		except:
			print getopt.getopt(args,':')[1][0],' not found'
			sys.exit(1)
		try:
			junkFrame = drx.readFrame(fh)
			try:
				srate = junkFrame.getSampleRate()
				pass
			except ZeroDivisionError:
				print 'zero division error'
				break
		except errors.syncError:
			print 'assuming the srate is 19.6 MHz'
			fh.seek(-drx.FrameSize+1, 1)
		fh.seek(-drx.FrameSize, 1)
		beam,tune,pol = junkFrame.parseID()
		beams = drx.getBeamCount(fh)
		tunepols = drx.getFramesPerObs(fh)
		tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
		beampols = tunepol
		if offset != 0:
			fh.seek(offset*drx.FrameSize, 1)
		if nChunks == 0:
			nChunks = 1
		nFrames = nFramesAvg*nChunks
		centralFreq1 = 0.0
		centralFreq2 = 0.0
		for i in xrange(4):
			junkFrame = drx.readFrame(fh)
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
		fh.seek(-4*drx.FrameSize, 1)
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
