import sys
sys.path.insert(1, '../mytools/')
import numpy as np
import matplotlib.pyplot as plt
import clock_tools as ctl

def fun(x,a,b,c):
    return a*(x**2)+ b*x + c


old_ml1 = [[0.2, 29.5, 0.20708 , 3 ],
           [0.5, 906.86233, 0.23224, 3],
           [0.9, 3285.65, 0.22776, 3],
           [1.4, 8198.82795, 0.22858, 3] ]
old_ml2 = [[0.2, 29.8067, 0.2052 , 4 ],
           [0.4, 523.7828, 0.2251, 2],
           [0.7, 1918.7264, 0.2453, 2] ]

new_ml1 = [ [0.3, 7.7694*4 ], [0.5, 174.3835*4], [0.7, 426.7562*4]  ]
new_ml2 = [ [0.4, 77.77543*4],[0.8, 581.80108*4],[0.3, 5.58503*4] ]
new_ml3 = [ [0.4, 67.240968*4],[0.8, 571.51599*4],[0.3, -5.302933*4] ]

#fitumitu( join([old_ml1, old_ml2], 0.2) )
#fitumitu( join([new_ml1, new_ml2, new_ml3], 0.3) )

s = ctl.shift(mltab=[new_ml1, new_ml2, new_ml3], jval=0.3)
s.set_fun(fun)
s.fit()
print('QZeeman shift for 0.2  is %.2f Hz'%( s.get_shift(0.2)))
print('QZeeman shift for 0.25 is %.2f Hz'%( s.get_shift(0.25)))
print('QZeeman shift for 0.26 is %.2f Hz'%( s.get_shift(0.26)))
print('QZeeman shift for 0.3  is %.2f Hz'%( s.get_shift(0.3)))
s.plot()
