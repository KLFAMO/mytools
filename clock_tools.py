import sys
sys.path.insert(1,'../mytools/')
import numpy as np
import matplotlib.pyplot as pqg
from scipy.optimize import curve_fit

class shift:
    def __init__(self, label='', mltab=[], jval=None):
        self.label = label
        self.join(mltab, jval)
        self.fun = None

    def set_fun(self,fun):
        self.fun = fun

    def set_coefs(self, coefs):
        self.coefs = coefs

    def norm_ml(self,tab,pval):
        for i in range(0, len(tab)):
            if tab[i][0] == pval:
                v = tab[i][1]
        for i in range(0,len(tab)):
            tab[i][1] = tab[i][1]-v
        return tab

    def join(self,tabs, pval):
        out = []
        dim = len(np.shape(tabs))
        if dim==2:
            tabs =[tabs]
        if pval == None:
            pval = tabs[0][0][0]
        out.append([pval, 0.0])
        for x in tabs:
            tmp = self.norm_ml(x,pval)
            for i in tmp:
                if i[0] != pval:
                    out.append([i[0], i[1]])
        self.mltab =  np.array(out).transpose()

    def fit(self,p0=None, maxfev=1000):
        popt, pcov = curve_fit(self.fun, self.mltab[0], self.mltab[1],
                    p0 = p0, maxfev=maxfev)
        print(popt)
        self.mltab[1] = self.mltab[1]-self.fun(0,*popt)
        popt, pcov = curve_fit(self.fun, self.mltab[0], self.mltab[1],
                    maxfev=maxfev)
        self.popt = popt
        self.pcov = pcov
        self.sigma_par = np.sqrt(np.diag(pcov))
        print("Fitting results")
        print("popt: ", popt)
        print("pcov: ", pcov)
        print("sigma: ", np.sqrt(np.diag(pcov)))
        
    def get_shift(self,x):
        return self.fun(x, *self.popt)

    def plot(self): 
        pqg.scatter(self.mltab[0], self.mltab[1], marker='+')
        x_fit = np.arange(0,max(self.mltab[0]),max(self.mltab[0])/100 )
        y_fit = self.fun(x_fit, *self.popt)
        pqg.plot(x_fit, y_fit)
        pqg.show()
