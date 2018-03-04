import numpy as np
import matplotlib.pyplot as plt
from apputils import forceIntValue

#Low tuning frequency range
Lfcl = 0
Lfch = 4095 * windownumber
#High tuning frequency range
Hfcl = 0
Hfch = 4095 * windownumber

# Get the range of FFT indices for the low and high tunings from the command line.
szLongOpts = ['low-tuning-lower', 'low-tuning-upper', 'high-tuning-lower', 'high-tuning-upper']
szShortOpts = 'abyz'
(cmdLnOpts, cmdLnParams) = getopt.getopt(sys.argv[1:], szShortOpts, szLongOpts)[0]
for index in range(len(cmdLnOpts)):
   # Set the FFT index for the specified tuning and force it to be a non-negative integer less
   # than 4096.
   if cmdLnOpts[index] in [szShortOpts[0], szLongOpts[0]]:
      Lfcl = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
   elif cmdLnOpts[index] in [szShortOpts[1], szLongOpts[1]]:
      Lfch = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
   elif cmdLnOpts[index] in [szShortOpts[2], szLongOpts[2]]:
      Hfcl = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
   elif cmdLnOpts[index] in [szShortOpts[3], szLongOpts[3]]:
      Hfch = forceIntValue(cmdLnParams[index], 0, 4095)*windownumber
   else
      print('UNKNOWN OPTION: {opt}'.format(opt=cmdLnOpts[index]))
      exit(1)
   # end if
# end for index in range(len(cmdLnOpts))
#
   sp = np.load('waterfall.npy')
plt.plot(sp[:,0,Lfcl:Lfch].mean(0))
plt.suptitle('Bandpass Low Tuning', fontsize = 30)
plt.ylabel('Intensity',fontdict={'fontsize':16})
plt.xlabel('Frequency',fontdict={'fontsize':14})
plt.show()
plt.plot(sp[:,1,Hfcl:Hfch].mean(0))
plt.suptitle('Bandpass High Tuning', fontsize = 30)
plt.ylabel('Intensity',fontdict={'fontsize':16})
plt.xlabel('Frequency',fontdict={'fontsize':14})
plt.show()
