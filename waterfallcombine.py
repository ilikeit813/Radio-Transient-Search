import glob
import numpy as np

fn = sorted(glob.glob('waterfall05*.npy'))
sp = np.zeros((len(fn),np.load(fn[0]).shape[0],np.load(fn[0]).shape[1]))
for i in range(len(fn)):
    sp[i,:,:]=np.load(fn[i])

np.save('waterfall',sp)
