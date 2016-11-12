'''
Use this code to examine if there are missing files.
'''
import numpy
import glob

nodes = 2
pps = 6
nChunks = 10000
nFramesAvg = 4

fn = sorted(glob.glob('waterfall05*.npy'))
j = numpy.zeros((len(fn)))
for i in range(len(fn)):
        j[i] = fn[i][39:48]

x = numpy.arange(j[-1]/nChunks/nFramesAvg)*nChunks*nFramesAvg
offset = numpy.setdiff1d(x, j)
