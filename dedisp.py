#script of Python script to dedisperse a file, save each DM trials time series
#Codes from Sean, modified by Jamie
#V3, super cluster technique

from mpi4py import MPI
import sys
import numpy as np
import glob
import os
import time
import smoth
import sys

def DMs(DMstart,DMend,dDM):
    """
    Calculate the number of DMs searched between DMstart and DMend, with spacing dDM * DM.

    Required:

    DMstart   - Starting Dispersion Measure in pc cm-3
    DMend     - Ending Dispersion Measure in pc cm-3
    dDM       - DM spacing in pc cm-3
    """

    NDMs = np.log10(float(DMend)/float(DMstart))/np.log10(1.0+dDM)

    return int(np.round(NDMs))


def delay2(freq, dm):
    """
    Calculate the relative delay due to dispersion over a given frequency
    range in Hz for a particular dispersion measure in pc cm^-3.  Return
    the dispersive delay in seconds.  Same as delay, but w.r.t to highest frequency.
    ***Used to simulate a dispersed pulse.***

    Required:

    freq - 1-D array of frequencies in MHz
    dm   - Dispersion Measure in pc cm-3
    """
    # Dispersion constant in MHz^2 s / pc cm^-3
    _D = 4.148808e3
    # Delay in s
    tDelay = dm*_D*((1/freq)**2 - (1/freq.max())**2)

    return tDelay


def Threshold(ts, thresh, clip=3, niter=1):
    """
    Wrapper to scipy threshold a given time series using Scipy's threshold function (in
    scipy.stats.stats).  First it calculates the mean and rms of the given time series.  It then
    makes the time series in terms of SNR.  If a given SNR value is less than the threshold, it is
    set to "-1".  Returns a SNR array with values less than thresh = -1, all other values = SNR.
    Also returns the mean and rms of the timeseries.

    Required:

    ts   -  input time series.

    Options:

    thresh  -  Time series signal-to-noise ratio threshold.  default = 5.
    clip    -  Clipping SNR threshold for values to leave out of a mean/rms calculation.  default = 3.
    niter   -  Number of iterations in mean/rms calculation.  default = 1.

    Usage:
    >>sn, mean, rms = Threshold(ts, *options*)

    """
    #  Calculate, robustly, the mean and rms of the time series.  Any values greater than 3sigma are left
    #  out of the calculation.  This keeps the mean and rms free from sturation due to large deviations.

    mean = np.mean(ts)
    std  = np.std(ts)
    #print mean,std

    if niter > 0:
        for i in range(niter):
            ones = np.where((ts-mean)/std < clip)[0]  # only keep the values less than 3sigma
            mean = np.mean(ts[ones])
            std  = np.std(ts[ones])
    SNR = (ts-mean)/std
    SNR[SNR<thresh]=-1
    #SNR = st.threshold((ts-mean)/std, threshmin=thresh, newval=-1)

    return SNR, mean, std

def Decimate_ts(ts, ndown=2):
    """
    Takes a 1-D timeseries and decimates it by a factore of ndown, default = 2.
    Code adapted from analysis.binarray module:
      http://www.astro.ucla.edu/~ianc/python/_modules/analysis.html#binarray
    from Ian's Python Code (http://www.astro.ucla.edu/~ianc/python/index.html)

    Optimized for time series' with length = multiple of 2.  Will handle others, though.

    Required:

    ts  -  input time series

    Options:

    ndown  -  Factor by which to decimate time series. Default = 2.
              if ndown = 1, returns ts
    """

    if ndown==1:
       return ts

    ncols = len(ts)
    n_rep = ncols / ndown
    ts_ds = np.array([ts[i::ndown][0:n_rep] for i in range(ndown)]).mean(0)

    return ts_ds


class OutputSource():

      pulse = None  # Pulse Number
      SNR   = None  # SNR of pulse
      DM    = None  # DM (pc/cm3) of pulse
      time  = None  # Time at which pulse ocurred
      dtau  = None  # Temporal resolution of time series
      dnu   = None  # Spectral resolution
      nu    = None  # Central Observing Frequency
      mean  = None  # Mean in the time series
      rms   = None  # RMS in the time series

      formatter = "{0.pulse:07d}    {0.SNR:10.6f}     {0.DM:10.4f}     {0.time:10.6f} "+\
                 "     {0.dtau:10.6f}     {0.dnu:.4f}     {0.nu:.4f}    {0.mean:.5f}"+\
                 "    {0.rms:0.5f}\n "

      def __str__(self):
          return self.formatter.format(self)

def bandpass(sp):
    ws = 151 # smooth windows size
    for pol in range(2):
        for i in range(3):
            a=smoth.smooth(np.median(sp[:,pol,:],0),ws)
            test=np.median(sp[:,pol,:],0)*1.0
            for i in range(1):
                test[np.where(test/a>1)]=a[np.where(test/a>1)]
                a=smoth.smooth(test,ws)
            sp[:,pol,:]/=a
    sp[:,0,:]=np.mean(sp[:,0:2,:],1) #X,Y poles are merged here
    return(sp[:,0,:])

def baseline(sp):
    """
    ws = 367 # smooth windows size
    sp = sp.transpose()
    for pol in range(2):
        for i in range(3):
            a=smoth.smooth(np.median(sp[:,pol,:],0),ws)
            test=np.median(sp[:,pol,:],0)*1.0
            for i in range(1):
                test[np.where(test/a>1)]=a[np.where(test/a>1)]
                a=smoth.smooth(test,ws)
            sp[:,pol,:]/=a
    return sp.transpose()
    """
    return sp

def RFI(sp):
    #for pol in range(2):
      #if pol == 0:
      #  sp[:,pol,np.where( np.abs(sp[:,pol,:].mean(0)-medianLf) > maxLf)] = medianLf
      #  sp[np.where( np.abs(sp[:,pol,:].mean(1)-medianLt) > maxLt),pol,:] = medianLt
      #else:
        #sp[:,pol,np.where( np.abs(sp[:,pol,:].mean(0)-medianHf) > maxHf)] = medianHf
        #sp[np.where( np.abs(sp[:,pol,:].mean(1)-medianHt) > maxHt),pol,:] = medianHt
    sp[:,np.where( np.abs(sp[:,:].mean(0)-medianHf) > maxHf)] = medianHf
    sp[np.where( np.abs(sp[:,:].mean(1)-medianHt) > maxHt),:] = medianHt
    return sp


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    fpp=2#spectrogram per processer
    nodes=34
    pps=10
    numberofFiles=fpp*nodes*pps
    maxpw=10 #Maximum pulse width to search in seconds. default = 1 s.
    dDM = 1.0/(3700-360) #1.0/len(freq)
    thresh=5.0 #SNR cut off
    #maxLf = 1.0451 #1.0481 #1.0451 #1.057 #1.062 #1.06  #1.071 #1.061 #1.061 #RFI cut off in frequency, at low, RFI is mostly in timeseries
    #maxHf = 1.0478 #1.0516 #1.0478 #1.057 #1.062 #1.045 #1.061 #1.061 #1.137 #RFI cut off in time domains, at high f, RFI is mostly in bandapss
    #maxLt = 1.0451 #1.0481 #1.0451 #1.057 #1.065 #1.06  #1.091 #1.061 #1.048 #RFI cut off in frequency
    #maxHt = 1.0478 #1.0516 #1.0478 #1.057 #1.065 #1.045 #1.061 #1.061 #1.061 #RFI cut off in time domains
    fn=sorted(glob.glob('056*.npy'))
    tInt=np.load('tInt.npy')

    DMstart=100.0#1.0 #initial DM trial
    DMend  =5000.0#90.0 #finial  DM trial

    fcl=360#low frequency cut off
    fch=3700#high frequency cut off

    npws = int(np.round(np.log2(maxpw/tInt)))+1 # +1 Due to in range(y) it goes to y-1 only


    spect=np.load(fn[0],mmap_mode='r')[:,:,fcl:fch]
    #spectarray = np.zeros((fpp,spect.shape[0],2,spect.shape[2])) # X and Y are merged already after bandpass
    spectarray = np.zeros((fpp,spect.shape[0],spect.shape[2])) # X and Y are merged already after bandpass

    #stdall=np.zeros((2,2,2,numberofFiles)) # 0:L, 1:H, 0,bp,1:ts, 0:std,1:median
    stdall=np.zeros((2,2,numberofFiles)) # 0,bp,1:ts, 0:std,1:median

    for i in range(fpp):
        #spectarray[i] = baseline( bandpass(np.asarray(np.load(fn[rank*fpp+i],mmap_mode='r')[:,:,fcl:fch],dtype=np.float64)))
        spectarray[i] = baseline( bandpass( np.asarray(np.load(fn[rank*fpp+i])[:,:,fcl:fch],dtype=np.float64)))
        for j in range(2):
            #for pol in range(2):
                #stdall[pol,j,0,rank*fpp+i] = spectarray[i,:,pol,:].mean(j).std()
                #stdall[pol,j,1,rank*fpp+i] = np.median(spectarray[i,:,pol,:].mean(j))
                stdall[j,0,rank*fpp+i] = spectarray[i,:,:].mean(j).std()
                stdall[j,1,rank*fpp+i] = np.median(spectarray[i,:,:].mean(j))


    stdalltmp=stdall*0
    comm.Allreduce(stdall,stdalltmp,op=MPI.SUM) #merge the 4 hour timeseries from all processor
    if rank ==0:
        np.save('jamieallstd',stdalltmp)

    #maxLf    = stdalltmp[0,:,0,:].min()*3
    #maxLt    = maxLf #stdalltmp[0,:,0,:].min()*3
    #medianLf = np.median(stdalltmp[0,0,1,:])
    #medianLt = medianLf #np.median(stdalltmp[0,1,1,:])
    #maxHf    = stdalltmp[1,:,0,:].min()*3
    #maxHt    = maxHf #stdalltmp[1,:,0,:].min()*3
    #medianHf = np.median(stdalltmp[1,0,1,:])
    #medianHt = medianHf #np.median(stdalltmp[1,1,1,:])
    maxHf    = stdalltmp[:,0,:].min()*3
    maxHt    = maxHf #stdalltmp[1,:,0,:].min()*3
    medianHf = np.median(stdalltmp[0,1,:])
    medianHt = medianHf #np.median(stdalltmp[1,1,1,:])



    #spect=np.append( bandpass(np.asarray((np.load(fn[rank*fpp],mmap_mode='r')[:,:,fcl:fch]),dtype=np.float64)) , bandpass(np.asarray((np.load(fn[rank*fpp+1],mmap_mode='r')[:,:,fcl:fch]),dtype=np.float64)),0)

    for i in range(fpp):
        spectarray[i] = RFI(spectarray[i])


    """#=====================================================================================================
    for pol in range(2):
        ts=np.zeros((numberofFiles*np.load(fn[0],mmap_mode='r').shape[0]))
        #ts[rank*spect.shape[0]:rank*spect.shape[0]+spect.shape[0]] += spect[:,pol,:].mean(1)
        for i in range(fpp):
            ts[ (rank*fpp+i)*spect.shape[0] : (rank*fpp+i+1)*spect.shape[0] ] += spectarray[i,:,pol,:].mean(1)

        tstotal=ts*0 #initiate a 4 hour blank time series
        comm.Allreduce(ts,tstotal,op=MPI.SUM) #merge the 4 hour timeseries from all processor
        if rank == 0 :
            np.save('ts_p%.1d_v10' % pol, tstotal) #cut the dispersed time lag

        bp=np.zeros((fch-fcl))
        bptotal=bp*0
        for i in range(fpp):
            bp+=spectarray[i,:,pol,:].mean(0)
        comm.Allreduce(bp,bptotal,op=MPI.SUM)
        if rank == 0:
            np.save('bp_p%.1d_v10'% pol, bptotal)


    sys.exit()
    """#======================================================================================================
    pol=0
    if pol==0:
    #for pol in (0,1):
        if pol==0:
            #freq=np.load('freq1.npy')
        #else:
            freq=np.load('freq2.npy')
            pol=1
        freq=freq[fcl:fch]/10**6
        cent_freq = np.median(freq)
        nDMs = DMs(DMstart,DMend,dDM)
        BW   = freq.max()-freq.min()
        DM=DMstart

        txtsize=np.zeros((npws,2),dtype=np.int32) #fileno = txtsize[ranki,0], pulse number = txtsize[ranki,1],ranki is the decimated order of 2
        txtsize[:,0]=1 #fileno star from 1

        tbmax=0 #repeat control, if dedispersion time series are identical, skip dedispersion calculation
        for iDM in range(0,nDMs):
            tb=np.round((delay2(freq,DM)/tInt)).astype(np.int32)

            if tb.max()-tbmax==0:#identical dedispersion time series checker
                tbmax=tb.max()
                DM+=dDM*DM
                #if rank ==0:
                #    print 'DM',DM,'skipped'
                continue
            tbmax=tb.max()

            ts=np.zeros((tb.max()+numberofFiles*np.load(fn[0],mmap_mode='r').shape[0]))
            for freqbin in range(len(freq)):
                for i in range(fpp):
                    #ts[tb.max()-tb[freqbin] + (rank*fpp+i)*spect.shape[0] :tb.max()-tb[freqbin] + (rank*fpp+i+1)*spect.shape[0] ] += spectarray[i,:,pol,freqbin]
                    ts[tb.max()-tb[freqbin] + (rank*fpp+i)*spect.shape[0] :tb.max()-tb[freqbin] + (rank*fpp+i+1)*spect.shape[0] ] += spectarray[i,:,freqbin]

            tstotal=ts*0#initiate a 4 hour blank time series
            comm.Allreduce(ts,tstotal,op=MPI.SUM)#merge the 4 hour timeseries from all processor
            tstotal = tstotal[tb.max():len(tstotal)-tb.max()]#cut the dispersed time lag

            #if rank==0:
            #   print '1',tb.max()*tInt
            #   print '2',tstotal.shape, tstotal.max(), tstotal.min(), tstotal.std()
            #   print 'DM',DM

            if rank<npws:#timeseries is ready for signal search
                ranki=rank
                filename = "pp_SNR_pol_%.1i_td_%.2i_no_%.05i.txt" % (pol,ranki,txtsize[ranki,0])
                outfile = open(filename,'a')
                ndown = 2**ranki #decimate the time series
                sn,mean,rms = Threshold(Decimate_ts(tstotal,ndown),thresh,niter=0)
                ones = np.where(sn!=-1)[0]
                #print '# of SNR>5',ones.shape,'decimated time',ndown
                #"""
                for one in ones:# Now record all pulses above threshold
                    pulse = OutputSource()
                    txtsize[ranki,1] += 1
                    pulse.pulse = txtsize[ranki,1]
                    pulse.SNR = sn[one]
                    pulse.DM = DM
                    pulse.time = one*tInt*ndown
                    pulse.dtau = tInt*ndown
                    pulse.dnu = freq[1]-freq[0]
                    pulse.nu = cent_freq
                    pulse.mean = mean
                    pulse.rms = rms
                    outfile.write(pulse.formatter.format(pulse)[:-1])
                    if txtsize[ranki,1] >200000*txtsize[ranki,0]:
                        outfile.close()
                        txtsize[ranki,0]+=1
                        filename = "pp_SNR_pol_%.1i_td_%.2i_no_%.05d.txt" % (pol,ranki,txtsize[ranki,0])
                        outfile = open(filename,'a')
                 #"""
            DM+=dDM*DM # End of DM loop
