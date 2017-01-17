import numpy
import glob

fn = sorted(glob.glob('05*.npy'))
j = numpy.zeros((len(fn)))
for i in range(len(fn)):
        j[i] = fn[i][30:39]

nChunks_nFramesAvg = j[1]-j[0]
#x = total perfect offset
x = numpy.arange(j[-1]/nChunks_nFramesAvg)*nChunks_nFramesAvg
# k = the different between perfect and real
k = numpy.setdiff1d(x, j)

for i in range(len(k)):
    a=fn[numpy.where( j == k[i]-nChunks_nFramesAvg)[0]]
    b=fn[numpy.where( j == k[i]+nChunks_nFramesAvg)[0]]
    newfn = fn[0][0:30]+'%.9i' % int(k[i]) +fn[0][39:]
    print newfn
    #numpy.save(newfn,numpy.load(a)*.5+numpy.load(b)*.5)
