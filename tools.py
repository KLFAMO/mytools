import numpy as np
import matplotlib.pyplot as plt
import decimal as dec
import time
import pyqtgraph as pg
from PyQt5.QtWidgets import QWidget
import allantools as al
from astropy.convolution import Gaussian1DKernel , convolve
import pytz
from time_tools import *
#import sqldata as sqd

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

class TimePeriod:
    def __init__(self, start, stop):
        self.start = min(start, stop)
        self.stop = max(start,stop)

    def __str__(self):
        return str(self.start) + " -> " + str(self.stop)
    
    def __mul__(self, b):
        if isinstance(b, TimePeriod):
            start = max(self.start, b.start)
            stop = min(self.stop, b.stop)
            if start < stop:
                return TimePeriod(start, stop)
        return None
    
    def joinIfOverlap(self, b):
        if isinstance(b, TimePeriod):
            if self*b != None:
                start = min(self.start, b.start)
                stop = max(self.stop, b.stop)
                return TimePeriod(start,stop)
        return None  

#--------------------------------------------------------------------------    

class TimePeriods:
    def __init__(self, periods=None):
        self.periods = []
        if periods!=None:
            for x in periods:
                self.append(x)

    def __str__(self):
        return "".join(['| '+x.__str__()+' ' for x in self.periods])
    
    def appendPeriod(self, periodToAdd):
        outputPeriods = []
        insertIndex = None
        for c, x in enumerate(self.periods):
            if x*periodToAdd == None:
                outputPeriods.append(x)
                if (insertIndex == None and 
                        periodToAdd.stop < x.start):
                    insertIndex = c
            else:
                if insertIndex == None:
                    insertIndex = c
                periodToAdd = periodToAdd.joinIfOverlap(x)
        if insertIndex == None:
            insertIndex = len(outputPeriods)
        outputPeriods.insert(insertIndex, periodToAdd)
        self.periods = outputPeriods
    
    def appendPeriods(self, TimePeriodsToAdd):
        for x in TimePeriodsToAdd.periods:
            self.appendPeriod(x)
    
    def append(self, b):
        if isinstance(b, TimePeriod):
            self.appendPeriod(b)
        if isinstance(b, TimePeriods):
            self.appendPeriods(b)
    
    def commonPart(self, b):
        if isinstance(b, TimePeriod):
            outPeriods=[]
            for x in self.periods:
                xb = x*b
                if xb != None:
                    outPeriods.append(xb)
            if len(outPeriods) == 0:
                return None
            return TimePeriods(periods=outPeriods)
        if isinstance(b, TimePeriods):
            outPeriods=[]
            for x in b.periods:
                cp = self.commonPart(x)
                if cp != None:
                    outPeriods.append(cp)
            if len(outPeriods) == 0:
                return None
            return TimePeriods(periods=outPeriods)
        return None
            


#-----------------------------------------------------------
class TSerie:
    def __init__(self, label='', mjd=[], val=[]):
        self.label=label
        self.mjd_tab=np.array(mjd)
        self.val_tab=np.array(val)
        self.len = len(self.mjd_tab)
        self.calc_tab()

    def __str__(self):
        if len(self.mjd_tab)==0:
            s = '\tEmpty'
        else:
            self.calc_tab()
            s = 'TSerie:\tlabel: %s\tlength: %d\tlen_mjd: %.6f\n'%(self.label,self.len,self.len_mjd)
        
            if self.len<=10:
                for i in range(0,self.len):
                    s = s+'\t%.6f\t%f\n'%(self.mjd_tab[i], self.val_tab[i])
            else:
                for i in range(0,5):
                    s = s+'\t%.6f\t%f\n'%(self.mjd_tab[i], self.val_tab[i])
                s = s+'\t...\n'
                for i in range(self.len-5,self.len):
                    s = s+'\t%.6f\t%f\n'%(self.mjd_tab[i], self.val_tab[i])
        return s

    def __add__(self, b):
        if isinstance(b, (int, float)):
            val = self.val_tab + b
            return TSerie(val = val, mjd = self.mjd_tab)
        elif isinstance(b, (TSerie)):
            print('Adding TSerie to TSerie is  not supported yet')
            return None
        else: 
            return None

    def __sub__(self, b):
        if isinstance(b, (int, float)):
            val = self.val_tab - b
            return TSerie(val = val, mjd = self.mjd_tab)
        elif isinstance(b, (TSerie)):
            print(' TSerie - TSerie is  not supported yet')
            return None
        else: 
            return None

    def __mul__(self, b):
        if isinstance(b, (int, float)):
            val = self.val_tab * b
            return TSerie(val = val, mjd = self.mjd_tab)
        elif isinstance(b, (TSerie)):
            print('TSerie * TSerie is  not supported yet')
            return None
        else: 
            return None

    def __truediv__(self, b):
        if isinstance(b, (int, float)):
            val = self.val_tab / b
            return TSerie(val = val, mjd = self.mjd_tab)
        elif isinstance(b, (TSerie)):
            print('TSerie / TSerie is  not supported yet')
            return None
        else: 
            return None

    def calc_tab(self): #if not empty
        self.len = len(self.mjd_tab)
        if self.len>0:
            self.isempty = 0
            self.mjd_start = self.mjd_tab[0]
            self.mjd_stop = self.mjd_tab[-1]
            self.mean = np.mean(self.val_tab)
        else:
            self.isempty = 1
            self.mjd_start = 0
            self.mjd_stop = 0
        mjd_s = 24*60*60
        self.s_tab = (self.mjd_tab - self.mjd_start)*mjd_s
        self.t_type = 't_type'
        self.len_mjd = self.mjd_stop-self.mjd_start
        self.len_s = self.len_mjd*mjd_s
    
    def len(self):
        return len(self.mjd_tab)

    def cp(self):
        out = TSerie(label=self.label+'_cp',
                     mjd = self.mjd_tab,
                     val = self.val_tab)
        return out

    def mean(self):
        if len(self.val_tab) > 0:
            return np.mean(self.val_tab)
        else:
            return None

    def rm_dc(self):
        self.val_tab = self.val_tab - self.mean

    def rm_drift(self):
        try:
            fit = np.polyfit(self.mjd_tab, self.val_tab,1)
            self.val_tab = ( self.val_tab -
                (self.mjd_tab*fit[0]+fit[1]) )
        except:
            pass
            # print('rm drift problem')

    def split(self, min_gap_s = 8):
        if self.len==0:
            return []
        out_tab = []
        tab_i = []
        tab_i.append(0)
        for i in range(0,len(self.s_tab)-1):
            if self.s_tab[i+1]-self.s_tab[i] > min_gap_s:
                tab_i.append(i+1)
        tab_i.append(len(self.s_tab)-1)
        #print(tab_i)
        for j in range(0,len(tab_i)-1):
            out_tab.append( TSerie(mjd=self.mjd_tab[tab_i[j]:tab_i[j+1]] ,
                                val=self.val_tab[tab_i[j]:tab_i[j+1]]    ))
        return out_tab

    def append(self, mjd, val):
        self.mjd_tab = np.append(self.mjd_tab,mjd)
        self.val_tab = np.append(self.val_tab,val)

    def last(self):
        out = TSerie()
        out.append(self.mjd_tab[-1], self.val_tab[-1])
        return out

    def mjd2index(self, mjd, init_index=None):
        """Returns index of tab corresponding to mjd
        in: mjd, init_index - initial index for searching
        out: index of the nearest mjd <= input_mjd
             None if mjd is out of table
        """
        if mjd < self.mjd_start or mjd > self.mjd_stop:
            return None
        if init_index==None:
            if self.len_mjd==0:
                return None
            N = int((self.len/self.len_mjd)*(mjd-self.mjd_start))
        else:
            N = init_index
        if N<0:
            N=0
        if N>=self.len:
            N = self.len-1
        while(1):
            if self.mjd_tab[N] > mjd:
                N = N-1
            else:
                if N+1>=self.len or self.mjd_tab[N+1] > mjd:
                    return N
                else:
                    N=N+1

    def mjd2val(self, mjd, init_index=None):
        return self.val_tab[self.mjd2index(mjd,init_index=init_index)]

    def getrange(self,fmjd,tmjd):
        if (fmjd > self.mjd_stop or tmjd < self.mjd_start):
            return None
        if fmjd<self.mjd_start:
            fmjd = self.mjd_start
        if tmjd>self.mjd_stop:
            tmjd = self.mjd_stop
        fN = self.mjd2index(fmjd)
        tN = self.mjd2index(tmjd)
        if tN==None:
            return None
        s =  TSerie( mjd = self.mjd_tab[fN:tN+1], 
                       val = self.val_tab[fN:tN+1])
        if s.mjd_tab[0]<fmjd:
            s.rm_first(1)
        if len(s.mjd_tab)==0:
            return None
        if s.mjd_tab[-1]>tmjd:
            s.rm_last(1)
        return s

    def rmrange(self, fmjd, tmjd):
        if (fmjd > self.mjd_stop or tmjd < self.mjd_start):
            return (self, 0)
        if (fmjd <= self.mjd_start and tmjd >= self.mjd_stop):
            return (None, 1)
        fN = self.mjd2index(fmjd)
        tN = self.mjd2index(tmjd)
        if (fmjd <= self.mjd_start and tmjd < self.mjd_stop):
            return (TSerie(mjd = self.mjd_tab[(tN+1):],
                          val = self.val_tab[(tN+1):])  , 2)
        if (fmjd > self.mjd_start and tmjd>= self.mjd_stop):
            if fmjd != self.mjd_tab[fN]: 
                fN = fN+1
            return (TSerie(mjd = self.mjd_tab[:fN],
                           val = self.val_tab[:fN]),   3)
        if (fmjd > self.mjd_start and tmjd < self.mjd_stop):
            if fmjd != self.mjd_tab[fN]: 
                fN = fN+1
            l = TSerie(mjd = self.mjd_tab[:fN],
                       val = self.val_tab[:fN])
            r = TSerie(mjd = self.mjd_tab[(tN+1):],
                       val = self.val_tab[(tN+1):])
            return ([l,r],4)

    def rm_outlayers_singledelta(self, max_delta):
        self.calc_tab()
        for i in range(1,self.len):
            if abs(self.val_tab[i]-
                    self.val_tab[i-1])> max_delta:
                print('outlayer detected')
                self.val_tab[i]=self.val_tab[i-1]
    
    def rmOutlayersOfTarget(self, target, maxDifference):
        i = 0
        #print('self.len: ', self.len)
        while i<self.len:
            if abs(self.val_tab[i]-target)>maxDifference:
                self.val_tab = np.delete(self.val_tab, i)
                self.mjd_tab = np.delete(self.mjd_tab, i)
                i = i-1
                self.len -=1
            i=i+1

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
        #plt.plot(self.mjd_tab, self.val_tab)
        plt.scatter(self.mjd_tab, self.val_tab*0, s=700, marker="|")
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
        self.val_tab = np.delete(self.val_tab, np.s_[0:n], None)
        self.mjd_tab = np.delete(self.mjd_tab, np.s_[0:n], None)
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

    def toMTSerie(self):
        return MTSerie(TSerie=self)
#---------------------------------------------------------------------------------
class MTSerie:
    def __init__(self, label='', TSerie=None, color='green', txtFileName=None):
        self.label = label
        self.dtab = []
        self.color = color
        if txtFileName != None:
            self.importFromTxtFile(txtFileName)
        if TSerie != None:
            self.add_TSerie(TSerie)

    def importFromTxtFile(self, fileName):
        f = open(fileName, 'r')
        mjd_t = []
        val_t = []
        for line in f:
            if line[0]!='#':
                mjd_t.append(float(line.split()[0]))
                val_t.append(float(line.split()[1]))
        f.close()
        self.add_TSerie(TSerie(mjd=mjd_t,val=val_t))
        self.split()

    def __str__(self):
        s = 'MTSerie:\n'
        for x in self.dtab:
            s = s+x.__str__()
        return s

    def __iadd__(self, b):
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                self.dtab[i] = self.dtab[i]+b
        return self
    
    def __add__(self, b):
        out = MTSerie()
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                out.add_TSerie( self.dtab[i]+b)
        return out


    def __sub__(self, b):
        out = MTSerie()
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                out.add_TSerie( self.dtab[i]-b)
        elif isinstance(b, (MTSerie)):
            ts = TSerie()
            tp1 = self.getTimePeriods()
            tp2 = b.getTimePeriods()
            tc = tp1.commonPart(tp2)
            for per in tc.periods:
                mjds = np.arange(per.start, per.stop, 1/(60*60*24))
                for mjd in mjds:
                    ts.append(mjd, self.mjd2val(mjd)-b.mjd2val(mjd))
            out = MTSerie()
            out.add_TSerie(TSerie(mjd=ts.mjd_tab, val=ts.val_tab))
            out.split()
        return out

    def __iadd__(self, b):
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                self.dtab[i] = self.dtab[i]+b
        return self

    def __imul__(self, b):
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                self.dtab[i] = self.dtab[i] * b
        return self

    def __mul__(self, b):
        out = MTSerie()
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                out.add_TSerie( self.dtab[i]*b)
        return out
    
    def __idiv__(self, b):
        if isinstance(b, (int, float)):
            for i in range(0,len(self.dtab)):
                self.dtab[i] = self.dtab[i] / b
        return self
    
    def isempty(self):
        for x in self.dtab:
            if x.len>0:
                return False
        return True
    
    def calc_tabs(self):
        for x in self.dtab:
            x.calc_tab()

    def rmemptyseries(self):
        self.calc_tabs()
        self.dtab = [x for x in self.dtab if x.isempty==0]

    def getTimePeriods(self):
        out = TimePeriods()
        for x in self.dtab:
            if len(x.mjd_tab)>0:
                out.appendPeriod(TimePeriod(x.mjd_tab[0],x.mjd_tab[-1]))
        return out

    def alphnorm(self,atom='88Sr'):
        d = {
            #ref: Safranova RevModPhys 2018  page 8
            '88Sr'  :{'fabs': dec.Decimal('429228066418000'), 'K':1.06},
            '87Sr'  :{'fabs': dec.Decimal('429228004229873'), 'K':1.06}, 
            '171Yb' :{'fabs': dec.Decimal('518295836590865'), 'K':1.31}, 
            '171Yb+':{'fabs': dec.Decimal('688358979309307'), 'K':2.03}, #dla E2  
                                                        # 1-5.95 dla E3
        }
        for x in d.keys():
            d[x]['factor'] = float(d[x]['fabs']*dec.Decimal(d[x]['K'])*dec.Decimal('1e-18') )
        if atom in d.keys(): 
            k = d[atom]['factor']
        else:
            return None
        for x in self.dtab:
                x.val_tab =  x.val_tab / k 

    def rm_drift_each(self):
        for x in self.dtab:
            x.rm_drift()
    
    def add_TSerie(self, ser):
        self.dtab.append(ser)

    def add_mjdf_data(self, mjd, f):
        tmp = TSerie(mjd=mjd, val=f)
        self.dtab.append(tmp)

    def add_mjdf_from_file(self, file_name):
        raw = np.load(file_name, allow_pickle=True)
        self.add_mjdf_data(raw[:,0], raw[:,1])

    def add_mjdf_from_datfile(self, file_name):
        raw = np.loadtxt(file_name)
        self.add_mjdf_data(raw[:,0], raw[:,1])

    def plot(self, color='', show=1, ax=None):
        for x in self.dtab:
            if color=='':
                color=self.color
            #plt.plot(x.mjd_tab, x.val_tab, color=color, linewidth=40)
            #plt.scatter(x.mjd_tab, x.val_tab, color=color, s=700, marker="|")
            if ax==None:
                plt.plot(x.mjd_tab, x.val_tab, color=color, marker="|")
            else:
                ax.plot(x.mjd_tab, x.val_tab, color=color, marker="|")
        if show==1:
            plt.show()

    def plot_pqg_widget(self, minusmjd=0, widget=None):
        if widget == None:
            self.widget = pg.PlotWidget()
        else:
            self.widget = widget
        for x in self.dtab:
            if len(x.mjd_tab)>0:
                print("---------")
                print(x)
                self.widget.plot(x.mjd_tab-minusmjd,x.val_tab)
        if widget == None:
            self.widget.setBackground('w')
            self.widget.setTitle(self.label)
        return self.widget

    def plot_allan(self, atom='88Sr'):
        if atom=='88Sr':
            fabs=429228066418012.0
        y=self.sew()/fabs
        #y = y.flatten()
        print('y: ',y)
        t = np.power(10, np.arange(0,int(np.log10(len(y)))+0.1,0.1))
        r = 1
        a = al.Dataset(data=y, rate=r, data_type="freq", taus=t) 
        a.compute('adev')
        b = al.Plot()
        b.plot(a, errorbars=True, grid=True)
        b.show()
        pass

    def sew(self, grid_s=1):
        g = grid_s/(24*60*60)
        out = []
        for x in self.dtab:
            for mjd in np.arange(x.mjd_tab[0],x.mjd_tab[-1],g):
                out.append( x.mjd2val(mjd) )
        return np.array(out)

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

    def rm_dc(self):
        mean = self.mean()
        for i in range(0,len(self.dtab)):
            self.dtab[i] = self.dtab[i]-mean

    def high_gauss_filter_each(self,stddev=50):
        i=0
        for x in self.dtab:
            i=i+1
            #print('filter '+str(i))
            if x.len < stddev*10:
                pass
                #print('to short')
                #self.dtab.remove(x)
            else:
                x.high_gauss_filter(stddev=stddev)

    def time_shift_each(self, sec):
        for x in self.dtab:
            x.time_shift(sec)

    def mjd2tabNo(self, mjd):
        for num in range(0,len(self.dtab)):
            if len(self.dtab[num].mjd_tab)>0 and mjd < self.dtab[num].mjd_tab[0]:
                return num-1
        return len(self.dtab)-1

    def mjd2tabNoandindex(self,mjd):
        tabNo = self.mjd2tabNo(mjd)
        if tabNo == -1:
            index = None
        else:
            index = self.dtab[tabNo].mjd2index(mjd)
        return (tabNo, index)

    def mjd2val(self, mjd, mode=0):
        (t,i) = self.mjd2tabNoandindex(mjd)
        if i == None:
            return None
        else:
            if mode == 1:
                v = self.dtab[t].val_tab[i]
                dv = self.dtab[t].val_tab[i+1]-v
                dm = self.dtab[t].mjd_tab[i+1]-self.dtab[t].mjd_tab[i]
                xm = mjd - self.dtab[t].mjd_tab[i]
                return v + dv * xm / dm
            else:
                return self.dtab[t].val_tab[i]

    def getrange(self,fmjd,tmjd):
        ft, fN = self.mjd2tabNoandindex(fmjd)
        tt, tN = self.mjd2tabNoandindex(tmjd)
        if fN != None:  fN = fN+1
        if (ft==tt and fN==None): return None
        out = MTSerie()
        if ft == -1: ft = 0
        if tt == -1: tt = 0

        for tab in self.dtab[ft:tt+1]:
            tmp = tab.getrange(fmjd,tmjd)
            if tmp != None:
                out.add_TSerie(tmp)
        if len(out.dtab)==0:
            return None
        return out

    def rmrange(self, fmjd, tmjd):
        ft, fN = self.mjd2tabNoandindex(fmjd)
        tt, tN = self.mjd2tabNoandindex(tmjd)
        if fN != None: fN = fN+1
        if (ft==tt and fN==None): return self
        if ft == -1: ft = 0
        if tt == -1: tt = 0
        i=tt
        while i>=ft:
            x = self.dtab[i]
            (a,b) = x.rmrange(fmjd,tmjd)
            if b==1:
                del self.dtab[i]
            if b==2 or b==3:
                self.dtab[i] = a
            if b==4:
                del self.dtab[i]
                self.dtab.insert(i,a[1])
                self.dtab.insert(i,a[0])
            i=i-1

    def rmoutlayers(self, max_iterations=4, target=None, maxdiff=None):
        i = 0
        no_rm = 0
        while(i<max_iterations and no_rm==0):
            i += 1
            print(i)
            input_len = len(self.mjd_tab())
            if target==None:
                target = self.mean()
            if maxdiff==None:
                maxdiff = 3*self.std()
            for x in self.dtab:
                x.rmOutlayersOfTarget(target=target, maxDifference=maxdiff)
            if input_len == len(self.mjd_tab()):
                no_rm=1

    def mean(self):
        if len(self.dtab)==0: return None
        tab = []
        totlen = 0
        for x in self.dtab:
            #print(x)
            tab.append([x.mean, x.len])
            totlen=totlen+x.len
        out=0
        for x in tab:
            out = out+ x[0]*x[1]/totlen
        return out
    
    def std(self):
        if len(self.dtab)==0: return None
        tab = []
        for x in self.dtab:
            tab = np.concatenate((tab,x.val_tab), axis=0)
        return np.std(tab)

    def stdm(self):
        if len(self.dtab)==0: return None
        tab = []
        totlen = 0
        for x in self.dtab:
            totlen = totlen + x.len
            tab = np.concatenate((tab,x.val_tab), axis=0)
        return np.std(tab)/np.sqrt(totlen)

    def len_mjd_eff(self):
        l = 0
        for x in self.dtab:
            l = l+x.len_mjd
            
    def mjd_tab(self):
        if len(self.dtab)>0:
            tmp = [np.array(x.mjd_tab) for x in self.dtab
                if np.array(x.mjd_tab).ndim>0]
            if len(tmp)>0:
                return np.concatenate( tmp )
        return []

    def val_tab(self):
        if len(self.dtab)>0:
            tmp = [np.array(x.val_tab) for x in self.dtab
                if np.array(x.val_tab).ndim>0]
            if len(tmp)>0:
                return np.concatenate( tmp )
        return []
    
    def saveToTxtFile(self, fileName):
        f = open(fileName, 'w')
        for x in self.dtab:
            i = 0
            length = len(x.mjd_tab)
            f.write('#\n')
            while i < length:
                f.write( '%f\t%f\n'%(x.mjd_tab[i], x.val_tab[i]))
                i=i+1
        f.close()


class MGserie:
    def __init__(self, ser1, ser2, grid_s, fmjd=None, tmjd=None, grid_mode=1):
        self.dtab=[]
        self.grid_s = grid_s
        self.grid_mjd=grid_s/(24*60*60)
        self.mjd_min = min(ser1.mjd_tab()[0], ser2.mjd_tab()[0])
        self.mjd_max = max(ser1.mjd_tab()[-1], ser2.mjd_tab()[-1])
        if fmjd != None:
            self.mjd_min = max(self.mjd_min, fmjd)
        if tmjd != None:
            self.mjd_max = min(self.mjd_max, tmjd)

        grid_tab = np.arange(self.mjd_min,self.mjd_max,self.grid_mjd)
        self.dtab = np.zeros((3,len(grid_tab)))
        for ind,mjd in enumerate(grid_tab):
            self.dtab[1][ind] = ser1.mjd2val(mjd, mode=grid_mode)
            self.dtab[2][ind] = ser2.mjd2val(mjd, mode=grid_mode)
        self.dtab[0] = grid_tab

    def corr(self, fmjd, tmjd, shift_s=None, norm=True):
        fi1 = self.mjd2index(fmjd)
        ti1 = self.mjd2index(tmjd)
        if shift_s==None:
            fi2=fi1
            ti2=ti1
        else:
            try:
                fi2 = fi1 + int(shift_s/self.grid_s)
                ti2 = ti1 + int(shift_s/self.grid_s)
            except:
                return np.nan
            if fi2<0 or ti2>len(self.dtab[0])-1 :
                return np.nan
            
        dm1 = self.dtab[0][fi1:ti1]
        d1 = self.dtab[1][fi1:ti1]
        dm2 = self.dtab[0][fi2:ti2]
        d2 = self.dtab[2][fi2:ti2]
        #plt.plot(dm1,d1)
        #plt.plot(dm2,d2)
        #plt.show()
        if len(d1)==len(d2):
            if norm==True:
                return np.correlate(d1,d2,'valid')[0]/len(d1)
            else:
                return np.correlate(d1,d2,'valid')[0]
        else:
            return np.nan
    


    def plot(self, ax=None, show=True, color='g'):
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        for x in [1,2]:
            d = np.array([[self.dtab[0][i],self.dtab[x][i]] for
                i in range(0,len(self.dtab[0])) 
                #if not np.isnan(self.dtab[x][i]) 
                ]).transpose()
            ax.plot(d[0],d[1], color=color)
        if show==True:
            plt.show()
        return ax

    def mjd2index(self, mjd, init_index=None):
        leni =len(self.dtab[0])
        len_mjd = self.dtab[0][-1]-self.dtab[0][0]
        if init_index==None:
            N = int((leni/len_mjd)*(mjd-self.dtab[0][0]))
        else:
            N = init_index
        if (mjd < self.dtab[0][0] or mjd > self.dtab[0][-1]) :
            return None
        while(1):
            if self.dtab[0][N] > mjd:
                N = N-1
            else:
                if self.dtab[0][N+1] > mjd:
                    return N
                else:
                    N=N+1
    
