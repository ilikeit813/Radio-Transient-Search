import os
import sys
import getopt
import glob
import numpy as np

def Decimate(ts, ndown=2):
   #return a decimated array shape = x, y, z, ...) with shape =  x/ndown, y, z, ....
    if ndown==1:
       return ts
    ncols = len(ts)
    n_rep = ncols / ndown
    ts_ds = np.array([ts[i::ndown][0:n_rep] for i in range(ndown)]).mean(0)
    return ts_ds



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
fn = sorted(glob.glob('waterfall{filename}*.npy'.format(filename=filename)))
sp = np.zeros((len(fn),np.load(fn[0]).shape[0],np.load(fn[0]).shape[1]))
for i in range(len(fn)):
    sp[i,:,:]=np.load(fn[i])

np.save('waterfall',Decimate(sp, sp.shape[0]/4000))
