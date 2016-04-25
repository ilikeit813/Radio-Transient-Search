import numpy as np
import glob
import sqlite3

con = sqlite3.connect("2016.sql") #2016.sql is an example, change it if necessary
cur = con.cursor()
cur.execute("CREATE TABLE 'lwa1' (snr FLOAT, dm FLOAT, time_tag FLOAT, resolution FLOAT)")
fn = sorted(glob.glob('*.txt'))
for i in fn:
    t=np.loadtxt(i)
    print 'doing', i
    [cur.execute("""INSERT INTO 'lwa1' (snr, dm, time_tag, resolution) VALUES (?,?,?,?)""", (j[1], j[2], j[3], j[4])) for j in t]

con.commit()#Do not forget this!!
