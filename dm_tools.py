from tools import MGserie, MTSerie
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scp
import sys
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[1]
sys.path.append(str(progspath / 'mytools'))

tcol = {'UMK1':'green',
        'UMK2':'red',
        'NIST':'blue',
        'NPLSr':'cyan',
        'NPLYb':'black',
        'NICT':'gray'}

class DMpair(MGserie):
    def __init__(self,lab1, lab2, grid_s=1, te_s=10, 
                fmjd=None, tmjd=None, grid_mode=1):
        d1 = prepforcalcpair(lab1,te_s)
        d2 = prepforcalcpair(lab2,te_s)
        self.te_s = te_s
        super().__init__(ser1=d1,ser2=d2,grid_s=grid_s,
                        fmjd=fmjd, tmjd=tmjd,
                        grid_mode=grid_mode)

    def calc_pair_sh(self, te_s, gshift_s):
        out_tab = []
        out_mjd = np.arange(self.mjd_min,self.mjd_max,te_s/(24*60*60))
        for i in out_mjd:
            out = self.CFDM(i,i+2*te_s/(24*60*60), shift_amp_s=3*self.te_s, 
                    shift_step_s=1, gshift_s=gshift_s)
            out_tab.append(out)
        o = np.array(out_tab)
        return np.nanmax(o)

    def CFDM(self, fmjd,tmjd, shift_amp_s=30, shift_step_s=1, gshift_s=0):
        s = np.arange(-shift_amp_s+gshift_s,
                       shift_amp_s+gshift_s+shift_step_s/2,
                       shift_step_s)
        c = [self.corr(fmjd,tmjd, shift_s=x) for x in s ]
        if np.isnan(c).any():
            return np.nan
        else:
            #fit
            ns = s - gshift_s
            (popt, pcov) = scp.curve_fit(self.ttriang, ns, c, 
                              bounds=([0,-np.inf],[np.inf,np.inf] ) )
            print('POPT:  ', popt)
            #plt.plot(ns,c)
            #plt.plot(ns, self.ttriang(ns, *popt) )
            #plt.show()
            return popt[0]
        
    def ttriang(self, x, A, gnd):
        return [triangle(z, A, self.te_s, gnd) for z in x]

def triangle( x, A, te, gnd):
        if np.abs(x) >= te:
                return gnd
        else:
                return  -np.abs( A*x/te  ) + A + gnd

def  prepforcalcpair(lab,te):
    d = MTSerie(lab, color=tcol[lab])
    d.add_mjdf_from_file(
        str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'.npy') )   )
    d.split(min_gap=12)
    d.rm_dc_each()
    d.high_gauss_filter_each(stddev=te*2)
    return d