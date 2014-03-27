import os
import numpy as np
import matplotlib.pyplot as plt
import glob

DMticks1=2.0   # ticks Number of DM trial in pulse number\SNR v.s. DM trial
DMticks2=4.0   # ticks Number of DM in time v.s. DM trial

fn = sorted(glob.glob('pp_p1*.txt'),key=os.path.getsize)
totaltime = 680*10000*0.0020897959
print totaltime/60/60.0,'Hr',totaltime,'sec',totaltime/20
for fname in range(len(fn)):
    DMticks3count=1

    DMticks3=1
    torg=np.loadtxt(fn[fname])
    sections=20
    for dec in range(sections):

    #for DMticks3 in (1.0/1000,1.0/100,1.0/10,1):
      #t=np.loadtxt(fn[fname])
      #for dec in range(5-2):
        #t=torg[torg[:,3].argsort()][torg.shape[0]/sections*dec:torg.shape[0]/sections*(dec+1)]
        t=torg[torg[:,3].argsort()][np.where((torg[torg[:,3].argsort()][:,3]> totaltime/sections*dec) & (torg[torg[:,3].argsort()][:,3]< totaltime/sections*(dec+1)))]
        savename='%.5s_%.1i_test' % (fn[fname],DMticks3count)
        DMticks3count+=1
        if t.shape ==0:
            continue
        t1=np.zeros((t.shape[0],3))
        for i in range(t.shape[0]):
                for j in range(3):
                    t1[i][j]=t[i][j+1]
        t1=t1[t1[:,1].argsort()]

        DMtrials=np.unique(np.copy(t1.T[1]))
        print 'number of DM trials recorded', DMtrials.shape
        dmhistorgm=t1.T[1][:]*0
        ii=0
        for i in range(1,len(t1.T[1])):
                if t1.T[1][i]==t1.T[1][i-1]:
                    dmhistorgm[i]=ii
                else:
                    ii+=1
                    dmhistorgm[i]=ii

        t1=t1.T
        t1[1]=dmhistorgm
        t1=t1.T
        t1=t1[t1[:,0].argsort()]
        t1=t1.T


        fig=plt.figure()
        #up left#
        fig.add_subplot(231).set_yscale('log')
        plt.ylabel('Number of Pulses')
        plt.xlabel('Signal to Noise')
        snhistog=np.histogram(t1[0][:],int((t1[0].max()-5)*10))
        snhistog[0][np.where(snhistog[0]==0)[0]]=1
        plt.plot( snhistog[1][1:], snhistog[0][:], color='black' )

        #Gaussian
        x0=1
        y=np.zeros((len(snhistog[1][x0:])))
        #print len(snhistog[1][x0:])-1
        for i in range(len(snhistog[1][x0:])-1):
            #print i,len(snhistog[1][x0:])
            y[i]=1.8*13072*(680*10000.0/sections)*(
                   np.exp(snhistog[1][x0:][i  ]**2/-2)/np.sqrt(2*np.pi)-
                   np.exp(snhistog[1][x0:][i+1]**2/-2)/np.sqrt(2*np.pi) )
            if y[i]<1: y[i]=0

        plt.plot(snhistog[1][x0:],y,'--')
        #Gaussian
        #plt.xticks(np.around(np.arange(int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()) - int(snhistog[1].min()))*(5)**-1 ),0))
        #print int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()+1) - int(snhistog[1].min()))*(1.0)
        plt.xticks(np.around(np.arange(int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()+1) - int(snhistog[1].min()))*(0.2) ),0))

        #up middle#
        fig.add_subplot(232).set_yscale('log')
        plt.ylabel('Number of Pulses')
        plt.xlabel(r'DM trials (pc cm$^{-3}$)')
        dmchannelhistog=np.histogram(t1[1][:],50)
        #dmchannelhistog=np.histogram(t1[1][:],len(DMtrials)/len(dmchannelhistog[1][x0:]) )
        #print len(DMtrials)/len(dmchannelhistog[1][x0:])
        #print "dmchannelhistog[0].shape",dmchannelhistog[0].shape
        #print "dmchannelhistog[1].shape",dmchannelhistog[1].shape
        #print "dmchannelhistog",dmchannelhistog
        #print 'DMtrials',DMtrials.shape


        #"""#Gaussian
        x0=0
        y=np.zeros((len(dmchannelhistog[1][x0:])))
        #print DMtrials.shape
        #print i,len(snhistog[1][x0:])
        #print len(dmchannelhistog[1][x0:])
        for i in range(len(dmchannelhistog[1][x0:])):
            yy=0
            for j in range(len(DMtrials)/len(dmchannelhistog[1][x0:])):
                #for j in range(13072/len(dmchannelhistog[1][x0:])):
                #if j+i*len(DMtrials)/len(dmchannelhistog[1][x0:])+1>DMtrials.shape[0]:
                #    continue
                #yy+=(313.4685-DMtrials[j+i*10]*2.37378)/313.4685*1.8*(680*10000)*(np.exp(5**2/-2)/np.sqrt(2*np.pi))
                yy+=1.8*(680*10000.0/sections)*(np.exp((5.0**2)/-2)/np.sqrt(2*np.pi))
            y[i]=yy
        #print y
        plt.plot(dmchannelhistog[1][x0:],y,'--')
        #"""#Gaussian

        plt.plot(dmchannelhistog[1][1:],dmchannelhistog[0],color='black')
        plt.xticks(np.around(np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max())),np.around(DMtrials[np.around(np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max())).astype('i')],1))

        #down#
        fig.add_subplot(212)
        plt.ylabel(r'DM trials(pc cm$^{-3}$)')
        plt.xlabel('time(sec)')
        cuts=t1.shape[1]*DMticks3
        #plt.scatter(t1[2][t1.shape[1]-cuts:t1.shape[1]],t1[1][t1.shape[1]-cuts:t1.shape[1]],s=t1[0][t1.shape[1]-cuts:t1.shape[1]]**6/3400,marker='o',facecolors='none')
        plt.scatter(t1[2][t1.shape[1]-cuts:t1.shape[1]],t1[1][t1.shape[1]-cuts:t1.shape[1]],s=t1[0][t1.shape[1]-cuts:t1.shape[1]],marker='o',facecolors='none')
        plt.yticks(np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks2**-1),dmhistorgm.max()),np.around(DMtrials[np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks2**-1),dmhistorgm.max()).astype('i')],1))

        t1=t1.T
        t1=t1[t1[:,1].argsort()]
        t1=t1.T

        #up right#
        fig.add_subplot(233)
        plt.ylabel('Signal to Noise')
        plt.xlabel(r'DM trials (pc cm$^{-3}$)')
        cuts=t1.shape[1]
        plt.plot(t1[1][t1.shape[1]-cuts:t1.shape[1]],t1[0][t1.shape[1]-cuts:t1.shape[1]],color='black')
        plt.yticks(np.around(np.arange(6,t1[0].max(),(t1[0].max()-6)/DMticks2),0))
        plt.xticks(
          np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max()),
          np.around(DMtrials[np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max()).astype('i')],1))

        plt.tight_layout()

        plt.savefig(savename)
        #plt.show()
