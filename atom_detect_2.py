import numpy as np
import matplotlib.pyplot as plt
import tools as tls
import sqldata as sqt
#--------------------------------------------------
class acc:
	def __init__(self, Len=10):
		self.Len = Len
		self.tab =np.array([])
		self.isfull = 0
		self.elements = 0
	
	def append(self, e):
		self.tab = np.append(self.tab, e)
		self.elements = len(self.tab)
		if self.elements>self.Len:
			self.tab = np.delete(self.tab,0,None)
			self.elements = len(self.tab)
	
	def print_tab(self):
		print(self.tab)
	
	def mean(self):
		return np.mean(self.tab)
	
	def std(self):
		return np.std(self.tab)
	
	def larger_than(self, th):
		out = 0
		for i in self.tab:
			if i > th:
				out=out+1
		return out
	
	def smaller_than(self, th):
		out = 0
		for i in self.tab:
			if i < th:
				out=out+1
		return out
	
	def out_of(self, th1, th2):
		min = min(th1,th2)
		max = max(th1,th2)
		out = 0
		for i in self.tab:
			if i < min or i> max:
				out=out+1
		return out
#------------------------------------------------------------

class detector:
	def __init__(self, from_mjd=57154.5, to_mjd=57158, sr=1): #from 57154.5 to 57158
		self.from_mjd = from_mjd
		self.to_mjd = to_mjd
		self.len_mjd = self.to_mjd - self.from_mjd
		self.len_s = self.len_mjd * 24*60*60
		self.sr = sr
		self.s_mjd = 1/(24*60*60)
		
	def get_data(self):

		self.aL_tab =sqt.getdata( 'sr'+str(self.sr)+'_lock_atoms_left',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
    
	def sim(self):
		self.atoms = tls.TSerie()
		self.avr_atoms = tls.TSerie()
		self.err = tls.TSerie()
		da = acc(30)
		for s in range(0, int(self.len_s), 2):
			mjd = self.from_mjd+s*self.s_mjd
		    # detection----
			
			self.avr_atoms.append(s,da.mean())
			x = self.aL_tab.mjd2val(mjd)
			if x==None: x=0
			da.append(x)

			if da.mean() > x/0.6 or da.mean() < 10000:
				self.err.append(s,x) 
			else:
				self.atoms.append(s, x)  
			# -------------
		
		plt.scatter(self.atoms.mjd_tab, self.atoms.val_tab, c='green', s=10, label='valid')
		plt.scatter(self.err.mjd_tab, self.err.val_tab, c='red', s=10 , label='error')
		plt.plot(self.avr_atoms.mjd_tab, self.avr_atoms.val_tab, c='orange', label='average of the last 30 cycles')
		plt.plot(self.avr_atoms.mjd_tab, self.avr_atoms.mjd_tab*0+10000, c='black', label='threshold for average')
		
		plt.grid()
		plt.xlabel('MJD')
		plt.ylabel('camera signal')
		plt.legend()
		plt.show()
					
x=detector(sr=1)
x.get_data()
x.sim()
