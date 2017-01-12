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

fn = sorted(glob.glob('waterfall05*.npy'))
sp = np.zeros((len(fn),np.load(fn[0]).shape[0],np.load(fn[0]).shape[1]))
for i in range(len(fn)):
    sp[i,:,:]=np.load(fn[i])

np.save('waterfall',Decimate(sp, sp.shape[0]/4000))
