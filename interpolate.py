import numpy
import glob

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
# end for

nChunks_nFramesAvg = j[1]-j[0]
#x = total perfect offset
x = numpy.arange(j[-1]/nChunks_nFramesAvg)*nChunks_nFramesAvg
# k = the different between perfect and real
k = numpy.setdiff1d(x, j)

for i in range(len(k)):
    a=fn[numpy.where( j == k[i]-nChunks_nFramesAvg)[0]]
    b=fn[numpy.where( j == k[i]+nChunks_nFramesAvg)[0]]
    newfn = fn[0][0:rootLength-1]+'%.9i.npy' % int(k[i]) 
    print newfn
    numpy.save(newfn,(numpy.load(a)+numpy.load(b))*.5)
