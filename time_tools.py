import astropy.time as ast
from datetime import datetime
from pytz import timezone
import time



def getMJD():
    return ast.Time(time.time(), format = 'unix').mjd

def getISO():
    return ast.Time(time.time(), format='unix').iso

def getLocalTime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime() )

def MJD2local(mjd, outstring="%Y-%m-%d %H:%M:%S"):
    u = ast.Time(mjd, format='mjd').unix
    lt = time.localtime(u)
    return time.strftime(outstring, lt)

def local2MJD(loc_str):
    ttab = loc_str.replace(':',' ').replace('-',' ').split()
    wawa = timezone('Europe/Warsaw')
    ttab = [int(x) for x in ttab]
    dt = datetime(ttab[0],ttab[1],ttab[2],
                  ttab[3],ttab[4],ttab[5])
    d = wawa.localize(dt, is_dst=False)
    unix = datetime.timestamp(d)
    mjd = ast.Time(unix, format='unix').mjd
    return mjd

def MJD2UTC(mjd, strfmt='%Y-%m-%d %H:%M:%S'):
    t = ast.Time(mjd, format='mjd')
    return t.strftime(strfmt)

def MJD2GPS(mjd, strfmt='%Y-%m-%d %H:%M:%S'):
    t = ast.Time(mjd, format='mjd') + ast.TimeDelta(18.0, format='sec')
    return t.strftime(strfmt)

def GPS2MJD(gpsstr,strfmt='%Y-%m-%d %H:%M:%S'):
    t = ast.Time.strptime(gpsstr,strfmt) - ast.TimeDelta(18.0,format='sec')
    return t.mjd

def GPSrange(fgps,tgps,grid,strfmt ):
    t = ast.Time.strptime(fgps,strfmt)
    tt = ast.Time.strptime(tgps,strfmt)
    gt = ast.TimeDelta(float(grid),format='sec')
    tab=[]
    while t<tt:
        tab.append(t.strftime(strfmt))
        t = t+gt
    return tab

