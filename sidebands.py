import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")

sqr2 = np.sqrt(2)
lamb = 813e-9
w0 =   152e-6
m = 87.62 * 1.66054e-27
h = 6.62607e-34  # J*Hz
kB = 1.380649e-23 # J/K

def lorentz(x, Amp, w, x0):
    return Amp*np.power(w,2)/(np.power((x-x0),2)+np.power(w,2))

def left_side_s(x, AL, AR, vz, Tr):
    vz = vz*1e3
    if abs(x)>vz:
        return 0
    if x<0:
        A=AL*1e-4
    else:   
        A=AR*1e-4

    Tr = Tr*1e-6
    U = m*vz*vz*lamb*lamb
    vr = vz*lamb/(sqr2*np.pi*w0 )
    C1 = A*4*U/(h*vr)
    C2 = 1 - abs(x)/vz
    C3 = np.exp( h*vr/(kB*Tr) )
    C4 = np.exp( (-4*U/(kB*Tr) * C2  ) )
    return C1*C2*C3*C4

def left_side(x, AL, AR,  vz, Tr):
    return np.array([left_side_s(i, AL, AR, vz, Tr) for i in x])


files = os.listdir('sidebands_data')
print(files)
dat = np.loadtxt('sidebands_data/'+files[0])
x = dat[:,0]*4e6 - 4*101.06e6
y = dat[:,1]

#find position of central peak
imax = np.argmax(y)
le = int(len(x)/10)
xlor = x[ imax-le:imax+le]
ylor = y[ imax-le:imax+le]
p0 = [y[imax], 1, x[imax]]
popt_lor, pcov_lor = curve_fit(lorentz, xlor, ylor, p0=p0)
print('*** Lorentz fit ****')
print('popt: ',popt_lor)

#shift data to central posiotion
x = x - popt_lor[2]

#fit left sideband
xL = x[:]
yL = y[:]
yL[int(imax)-20:int(imax)+20] = 0

p0 = [3, 1, 14*4, 1]
mi = (0, 0,    5*4,  0.5)
ma = (100, 100,  20*4, 200)
poptL, pcovL = curve_fit(left_side, xL, yL, p0=p0, bounds=(mi,ma))
print('*** Side-bands fit *************')
print('vz:  ', poptL[2])
print('Tr:  ', poptL[3], ' uK')
Tz = h*poptL[2]*1e3/(kB * np.log(poptL[0]/poptL[1]))
print('Tz:  ', Tz*1e6, 'uK') 


plt.scatter(x,y, s=1)
plt.plot(xlor, lorentz(xlor, popt_lor[0], popt_lor[1], 0), 'r')
plt.plot(xL, left_side(xL, *poptL), 'y')
plt.plot(xL, yL, 'g')
plt.show()
