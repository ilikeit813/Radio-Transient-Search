#This tool tells information about the data, need LSL library !!!

import os
import sys
import time
import ephem

from lsl.reader.ldp import LWA1DataFile
from lsl.writer import tsfits
from lsl.astro import unix_to_utcjd, DJD_OFFSET


def main(args):
    idf = LWA1DataFile(args[0])
    nFramesFile = idf.getInfo('nFrames')

    srate = idf.getInfo('sampleRate')
    beam = idf.getInfo('beam')
    beampols = idf.getInfo('beampols')

    # Date
    beginDate = ephem.Date(unix_to_utcjd(idf.getInfo('tStart')) - DJD_OFFSET)

    # File summary
    print "Filename: %s" % args[0]
    print "Date of First Frame: %s" % str(beginDate)
    print "Beam: %i" % beam
    print "Tune/Pols: %ii" % beampols
    print "Sample Rate: %i Hz" % srate
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
    print "---"


if __name__ == "__main__":
    main(sys.argv[1:])
                                                                                                                                                5,0-1         12%
