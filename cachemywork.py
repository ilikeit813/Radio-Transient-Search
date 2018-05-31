import os
import sys
import shutil
import glob
from optparse import OptionParser

cmdlnParser = OptionParser()
cmdlnParser.add_option('-w', '--work-dir', dest='workDir', default='.',
                        type='string', action='store',
                        help='Root directory where files are currently stored.',
                        metavar='PATH')
cmdlnParser.add_option('-c', '--cache-dir', dest='cacheDir', default='.',
                        type='string', action='store',
                        help='Root directory where files will be copied.',
                        metavar='PATH')
(cmdlnOpts, cmdlnParams) = cmdlnParser.parse_args()                        
workDir = cmdlnOpts.workDir
cacheDir = cmdlnOpts.cacheDir

if workDir != cacheDir:

   if not os.path.exists(cacheDir):
      print 'Cache directory does not exist.  It must be created before files can be cached.'
      sys.exit(1)
   # endif
   workFilepaths = glob.glob('{dir}/*.*'.format(dir=workDir))
   workFilenames = map(lambda x: os.path.split(x)[1], workFilepaths)
   cacheFilepaths = map(lambda x: '{dir}/{name}'.format(dir=cacheDir, name=x), workFilenames)

   nCount = 0
   numFiles = len(workFilepaths)
   for index in range(numFiles):
      workPath = workFilepaths[index]
      cachePath = cacheFilepaths[index]
      if os.path.exists(workPath):
         print 'Copying ', workPath
         shutil.copy(workPath, cachePath)
         if os.path.exists(cachePath):
            nCount = nCount + 1
            print 'Copy SUCCESS'
         else:
            print 'Copy FAIL'
         # endif
      else:
         print 'File =', workPath, 'not found.  Skipping.'
         numFiles = numFiles - 1
      # endif
   # end for

   print 'Successfully copied', nCount, 'files of', numFiles
else:
   print 'Work and cache directory are the same.  No files need to be copied.'
# endif

sys.exit(0)
