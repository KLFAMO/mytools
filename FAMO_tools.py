import sys

sys.path.insert(1, "../mytools/")
import sqldata as sqd
import numpy as np
import tools as tls

fthSr88 = 429228066418007.0


def calcN(fb, fr, f0, fth, fbsign=-1,f0sign=2):
    #fbs = int(fbsign + "1")
    fbs = fbsign
    #f0s = int(f0sign + "1")
    f0s = f0sign
    #print('calcN_ftheor: ', fth, '  sign: ',fbsign)
    return int(round((fth - (f0 * f0s + fbs * fb)) / fr, 0))


def fcomb(fb, fr, f0, fth, N=None, fbsign=-1, f0sign=2):
    #print('fcomb: theor: ',fth)
    if fbsign == None:
        x = [0, 0]
        x[0] = fcomb(fb, fr, f0, fth, N, fbsign=-1,f0sign=f0sign)
        x[1] = fcomb(fb, fr, f0, fth, N, fbsign=1,f0sign=f0sign)
        a = abs(x - fth)
        if a[0] < a[1]:
            return x[0]
        else:
            return x[1]
    #fbs = int(fbsign + "1")
    fbs = fbsign
    #f0s = int(f0sign + "1")
    f0s = f0sign
    if N == None:
        N = calcN(fb, fr, f0, fth, fbsign, f0sign=f0sign)
    return N * fr + f0 * f0s + fbs * fb


def fabsSr1(fb, fr, f0, fAOMSr1, N=None, fbsign=None, comb=1):
    if comb == 1:
        f0sign = 2
        fsr = fAOMSr1 * 4 - 80e6
    if comb == 2:
        f0sign = 2
        fsr = fAOMSr1 * 4 - 80e6 +80e6
    fc = fcomb(fb, fr, f0, fthSr88 + fsr, N, fbsign, f0sign= f0sign)
    return fc - fsr


def fabsSr2(fb, fr, f0, fAOMSr2, N=None, fbsign=None):
    fsr = fAOMSr2 * 2
    fc = fcomb(fb, fr, f0, fthSr88 + fsr, N, fbsign)
    return fc - fsr


def fabsSr1_range(fmjd, tmjd, grid=1):

    fb = sqd.getdata("comb_beat_698", fmjd, tmjd)
    fb.split()
    tfb = fb.getTimePeriods()
    fr = sqd.getdata("comb_f_rep", fmjd, tmjd)
    fr.split()
    tfr = fr.getTimePeriods()
    f0 = sqd.getdata("comb_f_zero", fmjd, tmjd)
    f0.split()
    tf0 = f0.getTimePeriods()
    fSr1 = sqd.getdata("sr1_lock_f", fmjd, tmjd)
    fSr1.split()
    tfSr1 = fSr1.getTimePeriods()

    t = tfb.commonPart(tfr)
    t = t.commonPart(tf0)
    t = t.commonPart(tfSr1)
    print(t)

    mean = fr.mean()
    print(mean)
    fr = fr*0
    print(fr)
    fr = fr + mean
    print(fr)


    tpar = [fb, fr, f0, fSr1]
    isetab = [x.isempty() for x in tpar]
    if sum(isetab) > 0:
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[], val=[]))

    for i in range(0,len(tpar)):
        tpar[i]=sqd.rmerr(tpar[i],['sr1','comb','j'])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for tp in t.periods:
        for mjd in np.arange(tp.start, tp.stop, g):
            p = [x.mjd2val(mjd) for x in [fb, fr, f0, fSr1]]
            if not None in p:
                a = fabsSr1(*p)
                if abs(a - fthSr88) < 100e9:
                    mjd_t.append(mjd)
                    abs_t.append(a - fthSr88)
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=10)
    return ts

def fabsSr1_ml_range(fmjd, tmjd, ml=1, grid=1, fr_dds=None, f0_dds=None,rmerr=True, min_gap=16, comb=1):
    if comb==1:
        combstr = ''
    if comb==2:
        combstr = '2'

    #print('start fabsSr1_ml_range')
    fb = sqd.getdata("comb%s_beat_698"%(combstr), fmjd, tmjd)
    fb.split(min_gap=min_gap)
    tfb = fb.getTimePeriods()
    #print('fb 1: ',fb)

    if fr_dds==None:
        fr = sqd.getdata("comb%s_f_rep"%(combstr), fmjd, tmjd)
        fr.split(min_gap=min_gap)
        #fr.plot()
        tfr = fr.getTimePeriods()
    else:
        fr = fb * 0 + fr_dds
    #print(fr)

    if f0_dds==None:
        f0 = sqd.getdata("comb%s_f_zero"%(combstr), fmjd, tmjd)
        f0.split(min_gap=min_gap)
        #f0.plot()
        tf0 = f0.getTimePeriods()
    else:
        f0 = fb *0 + f0_dds
    #print('fb 2: ',fb)
    #print('fr: ',fr)
    #print('f0: ',f0)

    if comb==2:
        print('fbeat = ',fb.mean())
        print('frep  = ',fr.mean())
        print('fzero = ',f0.mean())

    fSr1 = sqd.getdata("sr1_ml%d_f"%(ml), fmjd, tmjd)
    fSr1.split(min_gap=min_gap)
    tfSr1 = fSr1.getTimePeriods()
    
    try:
        t = tfb.commonPart(tfSr1)
        if f0_dds==None:
            t = t.commonPart(tf0)
        if fr_dds==None:
            t = t.commonPart(tfr)
    except:
        print('No data')
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[],val=[]))

    if t==None:
        print('No data')
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[],val=[]))
        

    print('time periods: ',t)

    tpar = [fb, fr, f0, fSr1]
    isetab = [x.isempty() for x in tpar]
    if sum(isetab) > 0:
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[], val=[]))

    if rmerr==True:
        for i in range(0,len(tpar)):
            tpar[i]=sqd.rmerr(tpar[i],['sr1','comb'])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for tp in t.periods:
        for mjd in np.arange(tp.start, tp.stop, g):
            p = [x.mjd2val(mjd) for x in [fb, fr, f0, fSr1]]
            if not None in p:
                a = fabsSr1(*p, comb=comb)
                #print(a)
                if abs(a - fthSr88) < 100e6:
                    mjd_t.append(mjd)
                    abs_t.append(a - fthSr88)
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=min_gap)
    return ts


def fabsSr1_ml1_range(fmjd, tmjd, grid=1, fr_dds=None, f0_dds=None):

    fb = sqd.getdata("comb_beat_698", fmjd, tmjd)
    fb.split()
    tfb = fb.getTimePeriods()
    print('fb 1: ',fb)

    if fr_dds==None:
        fr = sqd.getdata("comb_f_rep", fmjd, tmjd)
        fr.split()
        fr.plot()
        tfr = fr.getTimePeriods()
    else:
        fr = fb * 0 + fr_dds
    #print(fr)

    if f0_dds==None:
        f0 = sqd.getdata("comb_f_zero", fmjd, tmjd)
        f0.split()
        f0.plot()
        tf0 = f0.getTimePeriods()
    else:
        f0 = fb *0 + f0_dds
    print('fb 2: ',fb)
    print('fr: ',fr)
    print('f0: ',f0)

    fSr1 = sqd.getdata("sr1_ml1_f", fmjd, tmjd)
    fSr1.split()
    tfSr1 = fSr1.getTimePeriods()
    
    try:
        t = tfb.commonPart(tfSr1)
        if f0_dds==None:
            t = t.commonPart(tf0)
        if fr_dds==None:
            t = t.commonPart(tfr)
    except:
        print('No data')
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[],val=[]))

    print('time periods: ',t)

    tpar = [fb, fr, f0, fSr1]
    isetab = [x.isempty() for x in tpar]
    if sum(isetab) > 0:
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[], val=[]))

    for i in range(0,len(tpar)):
        tpar[i]=sqd.rmerr(tpar[i],['sr1','comb','j'])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for tp in t.periods:
        for mjd in np.arange(tp.start, tp.stop, g):
            p = [x.mjd2val(mjd) for x in [fb, fr, f0, fSr1]]
            if not None in p:
                a = fabsSr1(*p)
                #print(a)
                if abs(a - fthSr88) < 100e9:
                    mjd_t.append(mjd)
                    abs_t.append(a - fthSr88)
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=10)
    return ts

def fabsSr1_ml1_range_old(fmjd, tmjd, grid=1):

    fb = sqd.getdata("comb_beat_698", fmjd, tmjd)
    fb.split()
    tfb = fb.getTimePeriods()

    fr = sqd.getdata("comb_f_rep", fmjd, tmjd)
    fr.split()
    tfr = fr.getTimePeriods()

    f0 = sqd.getdata("comb_f_zero", fmjd, tmjd)
    f0.split()
    tf0 = f0.getTimePeriods()

    fSr1 = sqd.getdata("sr1_ml1_f_lock", fmjd, tmjd)*1e6
    fSr1.split()
    tfSr1 = fSr1.getTimePeriods()
    
    try:
        t = tfb.commonPart(tfr)
        t = t.commonPart(tf0)
        t = t.commonPart(tfSr1)
    except:
        print('No data')
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[],val=[]))

    print('time periods: ',t)

    tpar = [fb, fr, f0, fSr1]
    isetab = [x.isempty() for x in tpar]
    if sum(isetab) > 0:
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[], val=[]))

    for i in range(0,len(tpar)):
        tpar[i]=sqd.rmerr(tpar[i],['sr1','comb','j'])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for tp in t.periods:
        for mjd in np.arange(tp.start, tp.stop, g):
            p = [x.mjd2val(mjd) for x in [fb, fr, f0, fSr1]]
            if not None in p:
                a = fabsSr1(*p)
                print(a)
                if abs(a - fthSr88) < 100e3:
                    mjd_t.append(mjd)
                    abs_t.append(a - fthSr88)
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=10)
    return ts

def fabsSr2_range(fmjd, tmjd, grid=1):

    fb = sqd.getdata("comb_beat_698", fmjd, tmjd)
    fb.split()
    tfb = fb.getTimePeriods()
    fr = sqd.getdata("comb_f_rep", fmjd, tmjd)
    fr.split()
    tfr = fr.getTimePeriods()
    f0 = sqd.getdata("comb_f_zero", fmjd, tmjd)
    f0.split()
    tf0 = f0.getTimePeriods()
    fSr2 = sqd.getdata("sr2_lock_f", fmjd, tmjd)
    fSr2.split()
    tfSr2 = fSr2.getTimePeriods()

    t = tfb.commonPart(tfr)
    t = t.commonPart(tf0)
    t = t.commonPart(tfSr2)
    print(t)

    tpar = [fb, fr, f0, fSr2]
    isetab = [x.isempty() for x in tpar]
    if sum(isetab) > 0:
        return tls.MTSerie(TSerie=tls.TSerie(mjd=[], val=[]))

    # for i in range(0,len(tpar)):
    #    tpar[i]=sqd.rmerr(tpar[i],['sr1','comb','j'])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for tp in t.periods:
        for mjd in np.arange(tp.start, tp.stop, g):
            print(mjd)
            p = [x.mjd2val(mjd) for x in [fb, fr, f0, fSr2]]
            if not None in p:
                a = fabsSr2(*p)
                if abs(a - fthSr88) < 10e3:
                    mjd_t.append(mjd)
                    abs_t.append(a - fthSr88)
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=10)
    return ts


def fabsSr1_range_2(fmjd, tmjd, grid=1, rel_to_th=1, min_gap=10):
    fb = sqd.getdata("comb_beat_698", fmjd, tmjd)
    fr = sqd.getdata("comb_f_rep", fmjd, tmjd)
    f0 = sqd.getdata("comb_f_zero", fmjd, tmjd)
    fSr1 = sqd.getdata("sr1_lock_f", fmjd, tmjd)

    tpar = [fb, fr, f0, fSr1]
    for i in range(0, len(tpar)):
        tpar[i] = sqd.rmerr(tpar[i], ["sr1", "comb", "j"])

    frm = fr.mean()
    f0m = f0.mean()
    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for mjd in np.arange(fmjd, tmjd, g):
        p = [x.mjd2val(mjd) for x in [fb, fSr1]]
        if not None in p:
            mjd_t.append(mjd)
            abs_t.append(fabsSr1(p[0], frm, f0m, p[1]))
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    if min_gap != None:
        ts.split(min_gap=10)
    if rel_to_th == 1:
        ts -= fthSr88
    return ts


def Sr2vsSr1(fAOMSr1, fAOMSr2):
    fsr1 = fAOMSr1 * 4 - 80e6
    fsr2 = fAOMSr2 * 2
    return fsr1 - fsr2


def Sr2vsSr1_range(fmjd, tmjd, grid=2):
    fSr1 = sqd.getdata("sr1_lock_f", fmjd, tmjd)
    fSr2 = sqd.getdata("sr2_lock_f", fmjd, tmjd)

    fSr1 = sqd.rmerr(fSr1, ["sr1", "comb"])
    fSr2 = sqd.rmerr(fSr2, ["sr2", "comb"])

    g = grid / (24 * 60 * 60)
    N = len(np.arange(fmjd, tmjd, g))
    mjd_t = np.zeros(N)
    val_t = np.zeros(N)

    for index, mjd in enumerate(np.arange(fmjd, tmjd, g)):
        p = [x.mjd2val(mjd) for x in [fSr1, fSr2]]
        if not None in p:
            mjd_t[index] = mjd
            val_t[index] = Sr2vsSr1(*p)
        else:
            mjd_t[index] = None
            val_t[index] = None

    ts = tls.MTSerie(
        TSerie=tls.TSerie(
            mjd=[x for x in mjd_t if x != None and ~np.isnan(x)],
            val=[x for x in val_t if x != None and ~np.isnan(x)],
        )
    )
    ts.split(min_gap=10)
    return ts


def fSr1vsCav(fAOMSr1, fded):
    return fAOMSr1 * 4 - 80e6 - fded


def Sr1vsCav_range(fmjd, tmjd, grid=1):
    fSr1 = sqd.getdata("sr1_lock_f", fmjd, tmjd)
    fded = sqd.getdata("cav_dds698_dedrift_f", fmjd, tmjd)

    tpar = [fSr1]
    for i in range(0, len(tpar)):
        tpar[i] = sqd.rmerr(tpar[i], ["sr1", "j"])

    g = grid / (24 * 60 * 60)
    mjd_t = []
    abs_t = []
    for mjd in np.arange(fmjd, tmjd, g):
        p = [x.mjd2val(mjd) for x in [fded, fSr1]]
        if not None in p:
            mjd_t.append(mjd)
            abs_t.append(fSr1vsCav(*p))
    ts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_t, val=abs_t))
    ts.split(min_gap=10)
    return ts


def monitor(self):
    return 0
