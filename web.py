import web
from web import form
import numpy as np
import shutil
import os
import sys
import glob
import sqlite3

#def guplot(timeL,timeH,SNRmin,mjd,DML,DMH,resolL,resolH):
def guplot(timeL,timeH,SNRmin,DML,DMH,resolL,resolH):
    import tempfile
    os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()    
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    mpl.rc('xtick',labelsize=8)
    mpl.rc('ytick',labelsize=8)


    fftresol = 0.0020897959183673468*16

    con = sqlite3.connect("2016.sql")
    cur = con.cursor()
    cur.execute("SELECT snr, dm, time_tag FROM 'lwa1' WHERE (time_tag> '%.3f' and time_tag< '%.3f' and snr > '%.3f' and dm > '%.4f' and dm < '%.4f' and resolution > '%.4f' and resolution < '%.4f')" %(timeL, timeH, SNRmin, DML, DMH, resolL, resolH) )
    #cur.execute("SELECT snr, dm, time_tag FROM 'lwa1' WHERE (time_tag> '%.3f' and time_tag< '%.3f' and mjd = '%.f' and snr > '%.3f' and dm > '%.4f' and dm < '%.4f' and resolution > '%.4f' and resolution < '%.4f')" %(timeL, timeH, mjd, SNRmin, DML, DMH, resolL, resolH) )
    #cur.execute("SELECT snr, dm, time_tag FROM 'lwa1' WHERE (time_tag> '%.f' and time_tag< '%.f'  and snr > %.f)" %(timeL,timeH,SNRmin) )
    t=np.array(cur.fetchall())
    print t.shape,t.max(),t.min()
    if len(t) ==0:
        shutil.copyfile('error.png', 'image.png')
        return

    torg=np.zeros(( len(t),3))
    for i in range(len(t)):
        for j in range(3):
            torg[i,j]=t[i][j]
    #torg=np.array(cur.fetchall())

    #a=[r[0] for r in cur.fetchall()]
    #qresult = np.array(cur.fetchall())

    DMticks1=2.0   # ticks Number of DM trial in pulse number\SNR v.s. DM trial
    DMticks2=4.0   # ticks Number of DM in time v.s. DM trial

    for fname in range(1):
        DMticks3count=1
        DMticks3=1
        t1=torg[torg[:,2].argsort()]
        DMticks3count+=1
        t1=t1[t1[:,1].argsort()]

        DMtrials=np.unique(np.copy(t1.T[1]))
        #print 'number of DM trials recorded', DMtrials.shape
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


        """#Gaussian
        x0=1
        y=np.zeros((len(snhistog[1][x0:])))
        #print len(snhistog[1][x0:])-1
        for i in range(len(snhistog[1][x0:])-1):
            #print i,len(snhistog[1][x0:])
            y[i]=1.8*13072*( (timeH-timeL)/fftresol )*(
                   np.exp(snhistog[1][x0:][i  ]**2/-2)/np.sqrt(2*np.pi)-
                   np.exp(snhistog[1][x0:][i+1]**2/-2)/np.sqrt(2*np.pi) )
            if y[i]<1: y[i]=0

        plt.plot(snhistog[1][x0:],y,'--')
        """#Gaussian

        #plt.xticks(np.around(np.arange(int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()) - int(snhistog[1].min()))*(5)**-1 ),0))
        #print int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()+1) - int(snhistog[1].min()))*(1.0)
        plt.xticks(np.around(np.arange(int(snhistog[1].min()),int(snhistog[1].max())+1, (int(snhistog[1].max()+1) - int(snhistog[1].min()))*(0.2) ),0))

        #up middle#
        fig.add_subplot(232).set_yscale('log')
        plt.ylabel('Number of Pulses')
        plt.xlabel(r'DM trials (pc cm$^{-3}$)')
        if t1[1][:].shape[0]>50:
            dmchannelhistog=np.histogram(t1[1][:],50)
        else:
            dmchannelhistog=np.histogram(t1[1][:],1)
        #dmchannelhistog=np.histogram(t1[1][:],len(DMtrials)/len(dmchannelhistog[1][x0:]) )
        #print len(DMtrials)/len(dmchannelhistog[1][x0:])
        #print "dmchannelhistog[0].shape",dmchannelhistog[0].shape
        #print "dmchannelhistog[1].shape",dmchannelhistog[1].shape
        #print "dmchannelhistog",dmchannelhistog
        #print 'DMtrials',DMtrials.shape


        """#Gaussian
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
                yy+=1.8*( (timeH-timeL)/fftresol )*(np.exp((5.0**2)/-2)/np.sqrt(2*np.pi))
            y[i]=yy
        #print y
        plt.plot(dmchannelhistog[1][x0:],y,'--')
        """#Gaussian
        plt.plot(dmchannelhistog[1][1:],dmchannelhistog[0],color='black')
        plt.xticks(np.around(np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max())),np.around(DMtrials[np.around(np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max())).astype('i')],1))

        #down#
        fig.add_subplot(212)
        plt.ylabel(r'DM trials(pc cm$^{-3}$)')
        plt.xlabel('time(sec)')
        cuts=t1.shape[1]*DMticks3
        #plt.scatter(t1[2][t1.shape[1]-cuts:t1.shape[1]],t1[1][t1.shape[1]-cuts:t1.shape[1]],s=t1[0][t1.shape[1]-cuts:t1.shape[1]]**6/3400,marker='o',facecolors='none')
        plt.scatter(t1[2][t1.shape[1]-cuts:t1.shape[1]],t1[1][t1.shape[1]-cuts:t1.shape[1]],s=(t1[0][t1.shape[1]-cuts:t1.shape[1]]-4.9)*5,marker='o',facecolors='none')
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
        print np.around(np.arange(6,t1[0].max(),(t1[0].max()-6)/DMticks2),0)
        #plt.yticks(np.around(np.arange(6,t1[0].max(),(t1[0].max()-6)/DMticks2),0))
        plt.xticks(
          np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max()),
          np.around(DMtrials[np.append(np.arange(0,dmhistorgm.max(),dmhistorgm.max()*DMticks1**-1),dmhistorgm.max()).astype('i')],1))

        plt.suptitle('%.4f<time(s)<%.4f with SNR>%.2f, %.1f<resol(ms)<%.1f' % (timeL, timeH,SNRmin,resolL*1000,resolH*1000))

        plt.subplots_adjust(wspace=0.3)
        #plt.subplots_adjust(left=0.1)
        #plt.subplots_adjust(right=0.8)
        #plt.tick_params(labelsize=4)
        #plt.show()
        plt.savefig('image')
        return


render = web.template.render('templates/')
urls = ('/', 'index')
app = web.application(urls, globals())

myform = form.Form( 
    form.Textbox("First time frame (in seconds >  0hr)",form.notnull), 
    form.Textbox("Final time frame (in seconds <  4hr)",form.notnull), 
    form.Textbox("S/N minimum (> 5)",
        form.notnull,
        form.regexp('\d+', 'Must be a digit'),
        form.Validator('Must be => 5', lambda x:np.round(float(x))>4)
                ),
    form.Textbox("DM minimum trials (in pc cm^-3 and >    0)",form.notnull),
    form.Textbox("DM maximum trials (in pc cm^-3 and < 5000)",form.notnull),
    form.Textbox("Temporal min. decimation (in second and > 0.0)",form.notnull),
    form.Textbox("Temporal max. decimation (in second and < 10.0)",form.notnull)#,
    #form.Textbox("mjd",form.notnull)

                  )

class index: 
    def GET(self): 
        form = myform()
        # make sure you create a copy of the form by calling it (line above)
        # Otherwise changes will appear globally
        return render.formtest(form)

    def POST(self): 
        form = myform()
        #guplot(float(-5),float(5),float(5),float(200))
        if not form.validates(): 
            return render.formtest(form)
        else:
            # form.d.boe and form['boe'].value are equivalent ways of
            # extracting the validated arguments from the form.
            # guplot(timeL, timeH,SNRmin,mjd)
            #guplot(float(form['time (in seconds > -1hr) Low'].value),float(form['time (in seconds <  4hr) High'].value),float(form['SNR minimum (> 5)'].value), float(form['mjd'].value), float(form['DM trials (in pc cm^-3 >  100)  Low'].value),float(form['DM trials (in pc cm^-3 < 5000) High'].value),float(form['resolL (in s >    0)  Low'].value),float(form['resolH (in s < 10.0) High'].value))
            guplot(float(form['First time frame (in seconds >  0hr)'].value),float(form['Final time frame (in seconds <  4hr)'].value),float(form['S/N minimum (> 5)'].value), float(form['DM minimum trials (in pc cm^-3 and >    0)'].value),float(form['DM maximum trials (in pc cm^-3 and < 5000)'].value),float(form['Temporal min. decimation (in second and > 0.0)'].value),float(form['Temporal max. decimation (in second and < 10.0)'].value))
            return open('image.png',"rb").read()
            #return "Success! TimeL: %s, TimeH: %s, , SNR_L: %s" % (form['time (seconds) Low'].value,form['time (seconds) High'].value,form['SNR minimum'].value)

if __name__=="__main__":
    web.internalerror = web.debugerror
    app.run()
