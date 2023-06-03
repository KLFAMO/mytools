import sys
sys.path.insert(1, '../mytools/')
import tools as tls
import sqldata as sqd
import time_tools as tit
import FAMO_tools as ftl

def report(day_mjd):
    fSrtheor = 429228066418007.0 

    fmjd = int(day_mjd)-18/(24*60*60)
    tmjd = fmjd+1
    path = '/home/stront/lutoslawskiSVN/Sr/programy_korytarz/anda/'
    fname = path+(tit.MJD2UTC(day_mjd,strfmt='%Y-%m-%d')+
            '_UMK_Sr1-HMAOS.dat')
    file = open(fname,'w')
    file.write('#Optical frequency of UMK\'s 88Sr lattice clock '+
               '(Sr1) against the maser (HMAOS) \n')
    file.write('#File generated at '+
            tit.MJD2UTC(tit.getMJD(),'%H:%M:%S')+
            ' on the '+
            tit.MJD2UTC(tit.getMJD(),'%Y-%m-%d')+
            ' by Piotr Morzynski  \n')
    file.write('#Script used to generate this file: report.py\n')
    file.write('#The fractional frequency is defined as: '+
                '((r/r0)-1)/1e-18\n')
    file.write('## r = absolute frequency of Sr1\n')
    file.write('## r0=429228066418007.0 (BIPM2017)\n')
    file.write('#Gravitational red shift correction '+
               '(yet to be applied): y_GRS_Sr1 = %.3f e-15\n'%((-2.34/fSrtheor)*1e15))
    file.write('#Time of first obs '+
            tit.MJD2GPS(fmjd,'%Y/%m/%d '+'00:00:00') +
                ' GPS (ahead of UTC by 18s)\n')
    file.write('#\n')
    file.write('# yyyy/mm/dd hh:mm:ss (GPS) | '+
               '((r / r0) - 1) /1e-18 | '+
               'confidence | '+
               'systematic frequency correction (yet to be applied) [Hz] - will be applied later\n')
    file.write('\n')

    fb = sqd.getdata('comb_beat_698',fmjd,tmjd)
    fr = sqd.getdata('comb_f_rep',fmjd,tmjd)
    f0 = sqd.getdata('comb_f_zero',fmjd,tmjd)
    fSr1 = sqd.getdata('sr1_lock_f',fmjd,tmjd)
    fSr1 *= 4
    fSr1 -= 80e6
    fb.split(min_gap = 10)
    fr.split(min_gap = 10)
    f0.split(min_gap = 10)
    fSr1.split(min_gap = 10)
    fb.plot(show=0, color='b') 

    el = sqd.get_err_logs(fmjd,tmjd,['sr1','comb','j'])
    for e in el:
        for x in [fb,fr,f0,fSr1]:
            if e[1]>e[0]:
                x.rmrange(e[0],e[1])

    fb.plot(show=1,color='g')

    timestamps = tit.GPSrange(tit.MJD2GPS(fmjd,'%Y/%m/%d')+' 00:00:00',
                          tit.MJD2GPS(fmjd+1,'%Y/%m/%d')+ ' 00:00:00',
                          30.0,
                          '%Y/%m/%d %H:%M:%S') 

    for ts in timestamps:
        file.write(ts)
        s = fcomb_ts(ts,fb,fr,f0, fSrtheor+fSr1.mean() )
        if s!=None:
            sr = s-fSr1.mean()
            file.write('\t%#.9g\t1'% ( ((sr/fSrtheor) - 1)*1e18 ))
            file.write('\t%.3f\t%.5f'% (sr-fSrtheor, tit.GPS2MJD(ts,'%Y/%m/%d %H:%M:%S')
) )
        else:
            file.write('\t%#.9g\t0' % (0.0 ) )
            file.write('\t%.3f\t%.5f'% (0,tit.GPS2MJD(ts,'%Y/%m/%d %H:%M:%S')
) )
        file.write('\t%.f\n' % (0)  )
    
    file.close()

def calcN(fb,fr,f0,fth, fbsign='-'):
    fbs = int(fbsign+'1')
    return int( round( (fth-(2*f0 + fbs*fb))/fr , 0 ))

def fabs(fb,fr,f0,fth,N=None, fbsign='-'):
    fbs = int(fbsign+'1')
    if N==None: N = calcN(fb,fr,f0,fth, fbsign) 
    return N*fr+2*f0+fbs*fb
    
def fcomb_ts(ts, fbs, frs, f0s, fth):
    fcomb=[0,0]
    t_start = tit.GPS2MJD(ts,'%Y/%m/%d %H:%M:%S')
    t_end = t_start + 30/(24*60*60)
    (fb, fr, f0)  = [x.getrange(t_start,t_end) for x in [fbs, frs, f0s] ]
    if (fb==None or fr==None or f0==None): return None
    (fbm,frm,f0m) =  [ x.mean() for x in [fb, fr, f0] ]
    fcomb[0] = fabs(fbm,frm,f0m,fth,fbsign='-')
    fcomb[1] = fabs(fbm,frm,f0m,fth,fbsign='+')
    dif = abs(fcomb-fth)
    print(dif)
    if dif[0]<dif[1]: out = fcomb[0]
    else: out = fcomb[1]
    print(out-fth)
    return out


