'''
Use this code to examine if there are missing files.
'''
import numpy
import glob

windownumber = 2
nodes = 2
pps = 6
nChunks = 1000
nFramesAvg = 4*windownumber

#fn = sorted(glob.glob('waterfall05*.npy'))
#j = numpy.zeros((len(fn)))
#for i in range(len(fn)):
#        j[i] = fn[i][39:48]

fn = sorted(glob.glob('05*.npy'))
j = numpy.zeros((len(fn)))
for i in range(len(fn)):
        j[i] = fn[i][30:39]

x = numpy.arange(j[-1]/nChunks/nFramesAvg)*nChunks*nFramesAvg
offset = numpy.setdiff1d(x, j)
print 'missing ', offset.shape, 'files'
