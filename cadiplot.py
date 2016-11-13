import numpy as np
import matplotlib.pyplot as plt

def snr(a):
    return (a-a.mean() )/a.std()

sp = np.load('candi057139_000656029_1_fft_offset_029943115_frames.npy')
plt.imshow(snr(sp[:,0,:]).T,origin = 'low',aspect ='auto',cmap='Greys_r')
plt.suptitle('Spectrogram Low Tuning', fontsize = 30)
plt.xlabel('Time (0.84 ms)',fontdict={'fontsize':16})
plt.ylabel('Frequency (1.2 kHz)',fontdict={'fontsize':14})
plt.colorbar().set_label('std',size=18)
plt.savefig('cadi')
