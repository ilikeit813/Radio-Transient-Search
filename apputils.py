# 
# apputils.py
#
# Purpose: Contain various utility functions that are used locally by the application.
#
import os
import sys


def forceIntValue(inValue, lower, upper):
   # Return inValue as an integer forced to be in the range from <lower> to <upper>
   result = max([lower, inValue])
   result = min([upper, result])
   result = int(result)
   return result
# end forceIntVale()


def Decimate(ts, ndown=2):
   """
   Takes a N dimensional array and decimates it by a factor of ndown, default = 2, along axis = 0
   Code adapted from analysis.binarray module: 
   http://www.astro.ucla.edu/~ianc/python/_modules/analysis.html#binarray 
   from Ian's Python Code (http://www.astro.ucla.edu/~ianc/python/index.html)
    
   Optimized for time series' with length = multiple of 2.  Will handle others, though.

   Required:
    
   ts  -  input time series

   Options:
    
   ndown  -  Factor by which to decimate time series. Default = 2.
   if ndown <= 1, returns ts       
   """
   #return a decimated array shape = x, y, z, ...) with shape =  x/ndown, y, z, ....
   if ndown <= 1:
      return ts
   n_rep = int(len(ts) / ndown)
   return np.array([ts[i::ndown][0:n_rep] for i in range(ndown)]).mean(0)
# end Decimate()

