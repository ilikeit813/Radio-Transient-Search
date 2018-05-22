import os
import sys
import glob
from optparse import OptionParser

cmdlnParser = OptionParser()
cmdlnParser.add_option('-w', '--work-dir', dest='workDir', default='.',
                        type='string', action='store',
                        help='Root directory where files to be deleted are stored.',
                        metavar='PATH')
cmdlnParser.add_option('-t', '--waterfall', dest='fWaterfall', default=False,
                        action='store_true', help='Flag denoting to delete coarse waterfall files.')
cmdlnParser.add_option('-c', '--combined-waterfall', dest='fCombWaterfall', default=False,
                        action='store_true', help='Flag denoting to delete coarse combined waterfall file.')
cmdlnParser.add_option('-d', '--detail-waterfall', dest='fDetWaterfall', default=False,
                        action='store_true', help='Flag denoting to delete detailed waterfall files.')
cmdlnParser.add_option('-s', '--spectrogram-images', dest='fImages', default=False,
                        action='store_true', help='Flag denoting to delete spectrogram images files.')
cmdlnParser.add_option('-i', '--integration', dest='fIntegration', default=False,
                        action='store_true', 
                        help='Flag denoting to delete waterfall integration description files.')
cmdlnParser.add_option('-a', '--all', dest='fAll', default=False,
                        action='store_true', help='Flag denoting to delete all radio transient files.')
(cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()                        
workDir = cmdlnOpts.workDir

myfiles = []
if cmdlnOpts.fDetWaterfall or cmdlnOpts.fAll:
   myfiles = myfiles + glob.glob('{dir}/master*.npy'.format(dir=workDir))

if cmdlnOpts.fWaterfall or cmdlnOpts.fAll:
   myfiles = myfiles + glob.glob('{dir}/waterfall*.npy'.format(dir=workDir))

if cmdlnOpts.fImages or cmdlnOpts.fAll:
   myfiles = myfiles + glob.glob('{dir}/*.png'.format(dir=workDir))

if cmdlnOpts.fCombWaterfall or cmdlnOpts.fAll:
   myfiles = myfiles + ['{dir}/combined_waterfall.npy'.format(dir=workDir)] 

if cmdlnOpts.fIntegration or cmdlnOpts.fAll:
   myfiles = myfiles + ['{dir}/tInt.npy'.format(dir=workDir)] + \
               glob.glob('{dir}/*tunefreq.npy'.format(dir=workDir))

nCount = 0
numFiles = len(myfiles)
for filepath in myfiles:
   if os.path.exists(filepath):
      print 'Attempting to remove', filepath
      os.remove(filepath)
      if not os.path.exists(filepath):
         nCount = nCount + 1
         print 'Removal SUCCESS'
      else:
         print 'Removal FAIL'
      # endif
   else:
      numFiles = numFiles - 1
   # endif
# end for

print 'Successfully removed', nCount, 'files of', numFiles
