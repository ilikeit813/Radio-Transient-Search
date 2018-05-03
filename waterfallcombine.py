import os
import sys
import getopt
import glob
import numpy as np
from optparse import OptionParser

from apputils import Decimate



# CCY - This modification allows specifying the original time series data file path on the command-line.
# This is to avoid having to hard-code the pattern of files used by waterfallcombine.py and make the
# module more scriptable.
#
usage = " Usage: %prog [options] <radio filepath>"
cmdlnParser = OptionParser(usage=usage)
cmdlnParser.add_option('-o','--outfile', dest='outFilepath', default='waterfall', type='string',
                        action='store',
                        help='Path to the output file (do not include the .npy extension)',
                        metavar='FILEPATH')
(cmdlnOpts, cmdlnArgs) = cmdlnParser.parse_args()
# Make sure that a path to the original data file has been provided.
if len(cmdlnArgs[0]) == 0:
   print('Path to the original data file must be provided to determine waterfall filename pattern.')
   exit(1)
# endif
filename = os.path.basename(os.path.splitext(cmdlnArgs[0])[0])
filenamePattern = 'waterfall{name}*.npy'.format(name=filename)
waterfallFilenames = sorted(glob.glob(filenamePattern))

numFiles = len(waterfallFilenames)
if numFiles > 0:
   sample = np.load(waterfallFilenames[0])
   sp = np.zeros((numFiles, sample.shape[0], sample.shape[1]))
   for i in range(numFiles):
      sp[i,:,:]=np.load(waterfallFilenames[i])

   np.save(cmdlnOpts.outFilepath, Decimate(sp, sp.shape[0]/6000))

   # Check that the combined file was created, safely.
   if not os.path.exists('{prefix}.npy'.format(prefix=cmdlnOpts.outFilepath)):
      print('Failed to make combined waterfall file \'{path}.npy\'.'.format(path=cmdlnOpts.outFilepath))
      exit(1)
   # endif
else:
   print('No waterfall files matching pattern \'{pattern}\' found.'.format(pattern=filenamePattern))
   exit(1)
# endif

exit(0)
