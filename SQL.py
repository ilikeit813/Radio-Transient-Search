import numpy as np
import glob
import sqlite3

con = sqlite3.connect("liu2.sql")
cur = con.cursor()
cur.execute("CREATE TABLE 'lwa1' (snr FLOAT, dm FLOAT, time_tag FLOAT, resolution FLOAT, mjd TEXT, freq_tuning TEXT, target TEXT)")
fn = sorted(glob.glob('*.txt'))
for i in fn:
    t=np.loadtxt(i)
    mjd='%.3s' % i
    freq_tuning='h'
    if i == '201.txt' or i == '365.txt' or i == '517.txt' or i == '631.txt' or i ==  '894.txt':
        target='M31'
    else:
        target='B0031-07'
    print 'doing', i, 'target=',target
    [cur.execute("""INSERT INTO 'lwa1' (snr, dm, time_tag, resolution, mjd, freq_tuning, target) VALUES (?,?,?,?,?,?,?)""", (j[1], j[2], j[3]-j[2]*0.6184458528265389, j[4], mjd, freq_tuning, target)) for j in t]

con.commit()#Do not forget this!!
