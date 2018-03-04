'''
Use this code to examine if there are missing files.
'''
import numpy
import glob
import getopt

windownumber = 2
nodes = 2
pps = 6
nChunks = 1000
nFramesAvg = 4*windownumber

#fn = sorted(glob.glob('waterfall05*.npy'))
#j = numpy.zeros((len(fn)))
#for i in range(len(fn)):
#        j[i] = fn[i][39:48]

# CCY - This modification allows specifying the original time series data file path on the command-line.
# This is to avoid having to hard-code the pattern of files used by eyexam.py.  This will also allow
# eyexam.py to be scriptable within a workflow.
#
filepath = getopt.getopt(sys.argv[1:], szShortOpts, szLongOpts)[1][0]
# Make sure that a path to the original data file has been provided.
if len(filepath) == 0:
   print('Path to the original data file must be provided')
   exit(1)
# end if
filename = os.path.basename(os.path.splitext(filepath)[0])
fn = sorted(glob.glob('{filename}*.npy'.format(filename=filename)))

# Extract the frame number from the filename.
j = numpy.zeros((len(fn)))
rootLength = len(filename)
for i in range(len(fn)):
   j[i] = fn[i][rootLength:-4]

x = numpy.arange(j[-1]/nChunks/nFramesAvg)*nChunks*nFramesAvg
offset = numpy.setdiff1d(x, j)
print 'missing ', offset.shape, 'files'
