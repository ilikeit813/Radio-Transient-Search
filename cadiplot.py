import numpy as np
import matplotlib.pyplot as plt
import glob

def snr(a):
    return (a-a.mean() )/a.std()
tint = 0.00020897959
fn = glob.glob('ca*.npy')
print 'ploting file =',fn[0], 'offset time =', float(fn[0][35:44])/4*tint, ' sec'
offset = float(fn[0][35:44])/4*tint
sp = np.load(fn[0])


plt.suptitle('Spectrogram Low Tuning Offset Time %s sec' % str(offset), fontsize = 15)

plt.subplot(2, 2, 3)
plt.imshow(snr(sp[:,0,:]).T,origin = 'low',aspect ='auto',cmap='Greys_r')
plt.xlabel('Time (0.84 ms)')
plt.ylabel('Frequency (1.2 kHz)')
plt.colorbar().set_label('std')
plt.locator_params(nbins=8)

plt.subplot(2, 2, 1)
plt.plot(snr(sp[:,0,:]).T.mean(0))
plt.xlabel('Time series (0.84 ms)')
plt.ylabel('Arbitrary Units')
plt.locator_params(nbins=4)

plt.subplot(2, 2, 4)
plt.plot(snr(sp[:,0,:]).T.mean(1),np.arange(sp.shape[2]))
plt.xlabel('Bandpass (Arbitrary Units)')
plt.ylabel('Frequency (1.2 kHz)')
plt.locator_params(nbins=4)


plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('cadi')
