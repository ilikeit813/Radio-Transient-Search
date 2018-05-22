import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser

from apputils import forceIntValue

# CCY - The following modification to the command-line adds options and flexibility to bandpasscheck.py
# to allow it to be more scriptable in the radio transient workflow.
#
usage="Usage: %prog [options] <waterfall file>"
cmdlnParser = OptionParser(usage=usage)
cmdlnParser.add_option('-a', '--low-tuning-lower',dest='lowTuneLower', default=0, type='int',
                        action='store',
                        help='Lower FFT index (between 0 and 4095) for the low tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-b', '--low-tuning-upper',dest='lowTuneUpper', default=4095, type='int',
                        action='store',
                        help='Upper FFT index (between 0 and 4095) for the low tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-y', '--high-tuning-lower',dest='highTuneLower', default=0, type='int',
                        action='store',
                        help='Lower FFT index (between 0 and 4095) for the high tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option('-z', '--high-tuning-upper',dest='highTuneUpper', default=4095, type='int',
                        action='store',
                        help='Upper FFT index (between 0 and 4095) for the high tuning frequency.',
                        metavar='INDEX')
cmdlnParser.add_option("-w", "--work-dir", dest="workDir", default=".",
                       action="store",
                       help="Working directory path.", metavar="PATH")
(cmdlnOpts, cmdlnArgs) = cmdlnParser.parse_args()
workDir = cmdlnOpts.workDir
# Check that a combined waterfall file has been provided and that it exists.
if len(cmdlnArgs[0]) > 0:
   if not os.path.exists(cmdlnArgs[0]):
      print('Cannot find the combined waterfall file \'{path}\''.format(path=cmdlnArgs[0]))
      exit(1)
   # endif
else:
   print('A path to the combined waterfall file must be given.')
   exit(1)
# endif

# Force the FFT index values to be between 0 and 4095.
Lfcl = forceIntValue(cmdlnOpts.lowTuneLower, 0, 4095)
Lfch = forceIntValue(cmdlnOpts.lowTuneUpper, 0, 4095)
Hfcl = forceIntValue(cmdlnOpts.highTuneLower, 0, 4095)
Hfch = forceIntValue(cmdlnOpts.highTuneUpper, 0, 4095)
# Check that the lower FFT indice are less than the upper FFT indices.
if Lfch <= Lfcl:
   print('Low tuning lower FFT index must be less than the upper FFT index.')
   exit(1)
# endif
if Hfch <= Hfcl:
   print('High tuning lower FFT index must be less than the upper FFT index.')
   exit(1)
# endif
#

# Create the raw, unfiltered spectrogram for each tuning.
#
# Low frequency bandpass.
sp = np.load(cmdlnArgs[0])
plt.plot(sp[:,0,Lfcl:(Lfch + 1)].mean(0))
plt.suptitle('Bandpass Low Tuning', fontsize = 30)
plt.ylabel('Intensity',fontdict={'fontsize':16})
plt.xlabel('Frequency',fontdict={'fontsize':14})
plt.savefig('{dir}/lowbandpass'.format(dir=workDir))
plt.clf()
# High frequency bandpass.
plt.plot(sp[:,1,Hfcl:(Hfch + 1)].mean(0))
plt.suptitle('Bandpass High Tuning', fontsize = 30)
plt.ylabel('Intensity',fontdict={'fontsize':16})
plt.xlabel('Frequency',fontdict={'fontsize':14})
plt.savefig('{dir}/highbandpass'.format(dir=workDir))
plt.clf()

exit(0)
