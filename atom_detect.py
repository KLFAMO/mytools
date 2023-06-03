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
	def __init__(self, from_mjd=57154.5, to_mjd=57158.0, sr=1): #to 57158
		self.from_mjd = from_mjd
		self.to_mjd = to_mjd
		self.len_mjd = self.to_mjd - self.from_mjd
		self.len_s = self.len_mjd * 24*60*60
		self.sr = sr
		self.s_mjd = 1/(24*60*60)
		
	def get_data(self):
		self.f_tab = sqt.getdata( 'sr'+str(self.sr)+'_lock_f',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
		self.aL_tab =sqt.getdata( 'sr'+str(self.sr)+'_lock_atoms_left',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
		self.aR_tab =sqt.getdata( 'sr'+str(self.sr)+'_lock_atoms_right',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
		self.pL_tab =sqt.getdata( 'sr'+str(self.sr)+'_lock_prob_left',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
		self.pR_tab =sqt.getdata( 'sr'+str(self.sr)+'_lock_prob_right',
			from_mjd=self.from_mjd, to_mjd=self.to_mjd)
	
	def sim(self):
		self.out_atoms = tls.TSerie()
		self.out_laser = tls.TSerie()
		da = acc(10)
		dl = acc(30)
		for s in range(0, int(self.len_s), 2):
			mjd = self.from_mjd+s*self.s_mjd
		    # detection----
			da.append(self.aL_tab.mjd2val(mjd))
			dl.append(self.pL_tab.mjd2val(mjd))
			if da.smaller_than(500) > 5:
				self.out_atoms.append(mjd,0)
			else:
				self.out_atoms.append(mjd,1)
            
			if ((dl.smaller_than(0.1) > 29 and dl.std()<0.1)
				or da.smaller_than(500) > 5):
				self.out_laser.append(mjd,0)
			else:
				self.out_laser.append(mjd,1)
			
			# -------------
		
		fig, axs = plt.subplots(4, sharex=True, gridspec_kw={'hspace':0, 
		                                        'height_ratios':[2,1,2,1]})
		axs[0].plot(self.aL_tab.mjd_tab, self.aL_tab.val_tab)

		axs[1].plot(self.out_atoms.mjd_tab, self.out_atoms.val_tab, 'r')

		axs[2].plot(self.pL_tab.mjd_tab, self.pL_tab.val_tab)
		axs[2].set_ylim(-0.1,0.5)

		axs[3].plot(self.out_laser.mjd_tab, self.out_laser.val_tab, 'r')

		for ax in axs:
			ax.label_outer()
			ax.grid()
		plt.show()
		

			
x=detector(sr=1)
x.get_data()
x.sim()
