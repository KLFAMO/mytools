import mysql.connector as db
from mysql.connector import errorcode
import numpy as np
import matplotlib.pyplot as plt
import time
import astropy.time as ast
import pyqtgraph as pg
from PyQt5.QtWidgets import QWidget
import allantools as al
from astropy.convolution import Gaussian1DKernel , convolve
from pytz import timezone
import pytz
from datetime import datetime
from sqldata import connect, gettables, dbquery

#-------------------------------------------------------------
class Variables:
    def __init__(self,widget=0):
        self.list = {}
        self.len=0
        self.widget=widget

    def append(self,obj):
        self.list.update({ obj.label : obj})
        self.len=len(self.list)

#-----------------------------------------------------------
class TSerie:
    def __init__(self, label='', mjd=[], val=[]):
        self.label=label
        self.mjd_tab=mjd
        self.val_tab=val
        self.len = len(self.mjd_tab)
        if self.len>0:
            self.calc_tab()

    def calc_tab(self):
        self.len = len(self.mjd_tab)
        if self.len>0:
            self.mjd_start = self.mjd_tab[0]
            self.mjd_stop = self.mjd_tab[-1]
        else:
            self.mjd_start = 0
            self.mjd_stop = 0
        mjd_s = 24*60*60
        self.s_tab = (self.mjd_tab - self.mjd_start)*mjd_s
        self.t_type = 't_type'
        self.len_mjd = self.mjd_stop-self.mjd_start
        self.len_s = self.len_mjd*mjd_s

        self.mean = np.mean(self.val_tab)

    def cp(self):
        out = TSerie(label=self.label+'_cp',
                     mjd = self.mjd_tab,
                     val = self.val_tab)
        return out

    def rm_dc(self):
        self.val_tab = self.val_tab - self.mean

    def rm_drift(self):
        fit = np.polyfit(self.mjd_tab, self.val_tab,1)
        self.val_tab = ( self.val_tab -
            (self.mjd_tab*fit[0]+fit[1]) )

    def split(self, min_gap_s = 8):
        out_tab = []
        tab_i = []
        tab_i.append(0)
        for i in range(0,len(self.s_tab)-1):
            if self.s_tab[i+1]-self.s_tab[i] > min_gap_s:
                tab_i.append(i+1)
        tab_i.append(len(self.s_tab)-1)
        print(tab_i)
        for j in range(0,len(tab_i)-1):
            out_tab.append( TSerie(mjd=self.mjd_tab[tab_i[j]:tab_i[j+1]] ,
                                val=self.val_tab[tab_i[j]:tab_i[j+1]]    ))
        return out_tab

    def append(self, mjd, val):
        self.mjd_tab = np.append(self.mjd_tab,mjd)
        self.val_tab = np.append(self.val_tab,val)

    def last(self):
        out = tls.TSerie()
        out.append(self.mjd_tab[-1], self.val_tab[-1])
        return out

    def mjd2index(self, mjd, init_index=None):
        if init_index==None:
            N = int((self.len/self.len_mjd)*(mjd-self.mjd_start))
        else:
            N = init_index
        if mjd < self.mjd_start or mjd > self.mjd_stop:
            return 0
        while(1):
            if self.mjd_tab[N] > mjd:
                N = N-1
            else:
                if self.mjd_tab[N+1] > mjd:
                    return N
                else:
                    N=N+1

    def mjd2val(self, mjd, init_index=None):
        return self.val_tab[self.mjd2index(mjd,init_index=init_index)]

    def rm_outlayers_singledelta(self, max_delta):
        self.calc_tab()
        for i in range(1,self.len):
            if abs(self.val_tab[i]-
                    self.val_tab[i-1])> max_delta:
                print('outlayer detected')
                self.val_tab[i]=self.val_tab[i-1]

    def time_shift(self, sec):
        self.mjd_tab = self.mjd_tab+sec/(24*60*60)
        self.calc_tab()

    def show(self):
        N = 5
        if len(self.mjd_tab) > 2*N :
            for x in range(0,5):
                print(self.mjd_tab[x], self.val_tab[x])
            print('...')
            for x in range(self.length-6, self.length-1):
                print(self.mjd_tab[x], self.val_tab[x])
        elif len(self.mjd_tab) > 0:
            for x in range(0,2*N-1):
                print(self.mjd_tab[x], self.val_tab[x])
        else:
            print('Empty')

    def plot(self):
        plt.figure()
        plt.title(self.label)
        plt.plot(self.mjd_tab, self.val_tab)
        plt.show()

    def plot_pqg(self, minusmjd):
        w = pg.plot(self.mjd_tab-minusmjd,self.val_tab)
        w.setBackground('w')
        w.setTitle(self.label)
        w.show()
        return w

    def plot_pqg_widget(self, minusmjd=0):
        self.widget = pg.PlotWidget()
        self.widget.plot(self.mjd_tab-minusmjd,self.val_tab)
        self.widget.setBackground('w')
        self.widget.setTitle(self.label)
        return self.widget

    def plot_nice(self, fig=0, ax_site='left'):
        if fig==0:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.mjd_tab, self.val_tab)
        plt.show()
        return fig, ax

    def plot_allan_pqg(self, atom='88Sr'):
        if atom=='88Sr':
            fabs=400e12
        w = pg.PlotWidget()
        t = np.power(10, np.arange(1,int(np.log10(self.len_s))+0.1,0.1))
        #t=np.logspace(np.log10(20),np.log10(self.len_s),50)
        y=self.val_tab
        r = self.len_s/self.len
        (t2, ad, ade, adn) = al.adev(y, rate=r, data_type="freq", taus=t)
        w.plot(t2, ad/fabs, symbol='o')
        w.showGrid(True,True)
        w.setLogMode(True,True)
        w.setBackground('w')
        w.setTitle(self.label+'_ADEV')
        return w

    def plot_allan(self, atom='88Sr'):
        if atom=='88Sr':
            fabs=429228066418012
        t = np.power(10, np.arange(1,int(np.log10(self.len_s))+0.1,0.1))
        y=self.val_tab/fabs
        r = self.len_s/self.len
        a = al.Dataset(data=y, rate=r, data_type="freq", taus=t) 
        a.compute('adev')
        b = al.Plot()
        b.plot(a, errorbars=True, grid=True)
        b.show()

    def scatter(self):
        plt.scatter(self.mjd_tab, self.val_tab)
        plt.show()

    def rm_first(self,n):
        self.val_tab = np.delete(self.val_tab, np.s_[0:n-1], None)
        self.mjd_tab = np.delete(self.mjd_tab, np.s_[0:n-1], None)
        self.calc_tab()

    def rm_last(self,n):
        self.val_tab = np.delete(self.val_tab, np.s_[-n:], None)
        self.mjd_tab = np.delete(self.mjd_tab, np.s_[-n:], None)
        self.calc_tab()

    def gauss_filter(self, stddev=50):
        g = Gaussian1DKernel(stddev=stddev)
        self.val_tab = convolve(self.val_tab,g)
        self.rm_first(stddev*5)
        self.rm_last(stddev*5)

    def high_gauss_filter(self, stddev=50):
        g = Gaussian1DKernel(stddev=stddev)
        tmp = convolve(self.val_tab,g)
        self.val_tab = self.val_tab - tmp
        self.rm_first(stddev*5)
        self.rm_last(stddev*5)


#-------------------------------------------------------------
class MTSerie:
    def __init__(self, name, color='green'):
        self.name = name
        self.dtab = []
        self.color = color

    def add_mjdf_data(self, mjd, f):
        tmp = TSerie(mjd=mjd, val=f)
        self.dtab.append(tmp)

    def add_mjdf_from_file(self, file_name):
        raw = np.load(file_name, allow_pickle=True)
        self.add_mjdf_data(raw[:,0], raw[:,1])

    def add_mjdf_from_datfile(self, file_name):
        raw = np.loadtxt(file_name)
        self.add_mjdf_data(raw[:,0], raw[:,1])

    def plot(self, color='', show=1):
        for x in self.dtab:
            if color=='':
                color=self.color
            plt.plot(x.mjd_tab, x.val_tab, color=color)
        if show==1:
            plt.show()

    def split(self, min_gap = 8):
        tmp_tab = []
        for a in self.dtab:
            spl = a.split(min_gap)
            for s in spl:
                tmp_tab.append(s)
        self.dtab = tmp_tab

    def rm_dc_each(self):
        for x in self.dtab:
            x.rm_dc()

    def high_gauss_filter_each(self,stddev=50):
        i=0
        for x in self.dtab:
            i=i+1
            print('filter '+str(i))
            if x.len < stddev*10:
                print('to short')
                #self.dtab.remove(x)
            else:
                x.high_gauss_filter(stddev=stddev)

    def time_shift_each(self, sec):
        for x in self.dtab:
            x.time_shift(sec)

    def mjd2tabNo(self, mjd):
        for num, x in enumerate( self.dtab):
            if (mjd >= x.mjd_start and mjd<=x.mjd_stop):
                return num
        return -1

class MGserie:
    def __init__(self, ser1, ser2, grid_s):
        g =get_overlap_mjd(ser1, ser2)
        self.dtab=[]
        self.grid_s = grid_s
        self.grid_mjd=grid_s/(24*60*60)

        for x in g:
            s1 = []
            s2 = []
            mt = []
            for i in np.arange(x[0],x[1],self.grid_mjd):
                s1.append(ser1.dtab[x[2]].mjd2val(i) )
                s2.append(ser2.dtab[x[3]].mjd2val(i) )
                mt.append(i)
            self.dtab.append(np.array([mt, s1, s2]))

    def plot(self):
        for x in self.dtab:
            plt.plot(x[0],x[1])
            plt.plot(x[0],x[2])
        plt.show()

    def mjd2tabNo(self, mjd):
        for num, x in enumerate (self.dtab):
            if (mjd >= x[0][0] and mjd<=x[0][-1]):
                    return num
        return None

    def mjd2tabandindex(self, mjd, init_index=None):
        tabi = self.mjd2tabNo(mjd)
        if tabi != None:
            tab = self.dtab[tabi]
        else:
            return (None, None)
        leni =len(tab[0])
        len_mjd = tab[0][-1]-tab[0][0]
        if init_index==None:
            N = int((leni/len_mjd)*(mjd-tab[0][0]))
        else:
            N = init_index
        if (mjd < tab[0][0] or mjd > tab[0][-1]) :
            return (None, None)
        while(1):
            if tab[0][N] > mjd:
                N = N-1
            else:
                if tab[0][N+1] > mjd:
                    return (tabi, N)
                else:
                    N=N+1

    def corr(self, fmjd, tmjd):
        t1,n1 = self.mjd2tabandindex(fmjd)
        t2,n2 = self.mjd2tabandindex(tmjd)
        if (t1!=t2 or t1==None or t2==None):
            print('Wrong mjd range')
            return None
        dm = self.dtab[t1][0][n1:n2]
        d1 = self.dtab[t1][1][n1:n2]
        d2 = self.dtab[t1][2][n1:n2]
        plt.plot(dm,d1)
        plt.plot(dm,d2)
        plt.show()
        c = np.correlate(d1,d2,'full')
        print(c)
        plt.plot(range(0,len(c)),c)
        plt.show()
        return 1


def get_overlap_mjd(s1, s2):  #MTserie as arguments
    t = []
    for num, x in enumerate(s1.dtab):
        t.append((x.mjd_start,11, num))
        t.append((x.mjd_stop,12, num))

    for num, x in enumerate(s2.dtab):
        t.append((x.mjd_start,21, num))
        t.append((x.mjd_stop,22, num))

    t = sorted(t, key = lambda s: s[0])
    o1 = 0
    o2 = 0
    o = 0
    tt = []
    for x in t:
        if x[1]==11:  o1 = 1
        elif x[1]==12:  o1 = 0
        elif x[1]==21:  o2 = 1
        elif x[1]==22:  o2 = 0
        if (o1+o2+o==2):
            tmp = (x[0], s1.mjd2tabNo(x[0]), s2.mjd2tabNo(x[0]))
            tt.append(tmp)
            o = (o+1)%2
    t_out = []
    for x in range(0,len(tt),2):
        t_out.append((tt[x][0],tt[x+1][0], tt[x][1], tt[x][2]))
    return t_out

#-------------------------------------------------------------



def getMJD_my():
    """return actual time in mjd format
    """
    return 40587 + time.time()/86400

def getMJD():
    return ast.Time(time.time(), format='unix').mjd

def getISO():
    return ast.Time(time.time(), format='unix').iso

def getLocalTime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def MJD2local( mjd):
    u =  ast.Time(mjd, format='mjd').unix
    lt = time.localtime(u)
    return time.strftime("%Y-%m-%d %H:%M:%S", lt)

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


def getdata(name,from_mjd=57091.999,to_mjd=1000000, last=0):

    if last>0:
        to_mjd = getMJD_my()
        from_mjd = to_mjd - last/86400

    results = dbquery( "select * from " + name+
                     " where mjd>"+str(from_mjd)+
                     " and mjd<"+str(to_mjd)+" ;")
    d = np.array(results).transpose()
    if (0 in np.shape(d)):
        print('No data in database')
        out = TSerie( name, mjd=np.array([]), val=np.array([]) )
    else:
        out = TSerie( name, mjd=d[0], val=d[1])
    return out


def sendMessage(mjd, mjd2, mes, tag):
    str =( "INSERT INTO logs (mjd, mjd2, mes, tag) VALUES ("+
            "'"+mjd +"','"
            +mjd2+"','"+
            mes+"','"+tag+"');" )
    print(str)
    dbsend(str)


def dbsend(querstr):
    con = connect()
    cur=con.cursor()
    if isinstance(querstr, list):
        for x in querstr:
            cur.execute(x)
    else:
        cur.execute(querstr)
    con.commit()
    cur.close()
    con.close()

def dbsend_tmvl(tmvl):
    con = connect()
    cur=con.cursor()
    if isinstance(tmvl, list):
        for x in tmvl:
            s = "INSERT INTO %s (mjd,val) VALUES ('%s','%s');" % (
                    x[0], x[1], x[2])
            try:
                cur.execute(s)
            except db.Error as err:
                if err.errno == 1146:
                    print('No table '+x[0]+' in database')
                else:
                    raise
    else:
        print('Error: dbsend_tmvl: Argument is not a list')
    con.commit()
    cur.close()
    con.close()

def db_create_table(name):
    q = ("CREATE TABLE IF NOT EXISTS "+
            name+
            " (mjd DOUBLE NOT NULL, val DOUBLE NOT NULL);" )
    dbsend(q)

    
def monitor(name='sr2',from_mjd=57091.999,to_mjd=1000000, last=0, every=1):
    if name=='sr2':
        f = getdata('sr2_lock_f',from_mjd=from_mjd,to_mjd=to_mjd,last=last)
        p = getdata('sr2_lock_prob_right',from_mjd=from_mjd,to_mjd=to_mjd,last=last)

    elif name=='sr1':
        f = getdata('sr1_lock_f',from_mjd=from_mjd,to_mjd=to_mjd,last=last)
        p = getdata('sr1_lock_prob_right',from_mjd=from_mjd,to_mjd=to_mjd,last=last)

    fig, ax1 = plt.subplots()
    ax1.scatter(p.mjd_tab[::every], p.val_tab[::every], color='#cccccc', marker='.', s=2)
    ax1.set_xlabel('mjd')
    ax1.set_ylabel('prob')
    ax1.set_ylim(0,1)
    ax1.grid()
    ax2 = ax1.twinx()
    ax2.plot(f.mjd_tab[::every], f.val_tab[::every], 'b-', linewidth=1.0)
    ax2.set_ylabel('f')
    fig.tight_layout()
    plt.show()

def Sr1vsAOS(from_mjd=58908.38, to_mjd=58909):
    ftheor = 429228066418012
    beat = getdata('comb_beat_698',from_mjd=from_mjd, to_mjd=to_mjd)
    print('mean(beat) = %.0f Hz  + 48 MHz' % (np.mean(beat.val_tab)-48e6))
    print('stdm = ', np.std(beat.val_tab)/np.sqrt(len(beat.val_tab))  )
    return beat.plot_allan_pqg()

#----main prog -------------------------
#monitor(name='sr2',from_mjd=58661,every=1,last=0)  #since last compaign
#monitor(name='sr1',from_mjd=57205.4, to_mjd=57206.7 ,every=0,last=0)

