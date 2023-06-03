import sys
sys.path.insert(1, '../mytools/')
sys.path.insert(1, '/home/stront/lutoslawskiSVN/Sr/programy_korytarz/mytools/')
import tools as tls
import pyqtgraph as pqg
import sqldata as sqd
import FAMO_tools as ftl

class QtMplot(pqg.GraphicsLayoutWidget):
    def __init__(self, fmjd=None, tmjd=None):
        super().__init__()
        self.fmjd = fmjd
        self.tmjd = tmjd
        self.monitor()
        self.setBackground('w')
        #self.layout.setSpacing(0.)
        #self.setContentsMargins(0.,0.,0.,0.)

    def monitor(self, commonscale=1):
        smjd = 58000.0
        self.plot1 = self.addPlot(2,0)
        self.plot1.showGrid(1,1,1)
        self.plot1.setLabel('left','<font size="1"><b>Sr1: g-gnd, r-exc</b></font>')
        self.plot1.getAxis('bottom').setStyle(showValues=False)
        self.plot2 = self.addPlot(0,0)
        self.plot2.showGrid(1,1,1)
        self.plot2.setLabel('left','<font size="1"><b>Sr2: g-gnd, r-exc</b></font>')
        self.plot2.getAxis('bottom').setStyle(showValues=False)
        if commonscale==1: self.plot2.setXLink(self.plot1)
        self.plot3 = self.addPlot(4,0)
        self.plot3.showGrid(1,1,1)
        self.plot3.setLabel('left','<font size="1"><b>Sr1 vs Sr2 </b></font>')
        self.plot3.getAxis('bottom').setStyle(showValues=False)
        if commonscale==1: self.plot3.setXLink(self.plot1)
        self.plot4 = self.addPlot(5,0)
        self.plot4.showGrid(1,1,1)
        self.plot4.setLabel('left','<font size="1"><b>Sr1 abs </b></font>')
        if commonscale==1: self.plot4.setXLink(self.plot1)
        self.plot5 = self.addPlot(3,0)
        self.plot5.showGrid(1,1,1)
        self.plot5.setLabel('left','<font size="1"><b>exc_prob_Sr1 </b></font>')
        if commonscale==1: self.plot5.setXLink(self.plot1)
        self.plot5.setYRange(-0.1,1)
        self.plot6 = self.addPlot(1,0)
        self.plot6.showGrid(1,1,1)
        self.plot6.setLabel('left','<font size="1"><b>exc_prob_Sr2 </b></font>')
        if commonscale==1: self.plot6.setXLink(self.plot1)
        self.plot6.setYRange(-0.1,1)

        d1 = sqd.getdata('sr1_lock_atoms_gnd',self.fmjd,self.tmjd)
        d2 = sqd.getdata('sr1_lock_atoms_exc',self.fmjd,self.tmjd)
        d1 = sqd.rmerr(d1,'sr1')
        d2 = sqd.rmerr(d2,'sr1')
        logs = sqd.get_err_logs(self.fmjd,self.tmjd,['sr1','j'])
        self.plot1.plot(d1.mjd_tab()-smjd
                , d1.val_tab(), pen=pqg.mkPen('g'))
        self.plot1.plot(d2.mjd_tab()-smjd
                , d2.val_tab(), pen=pqg.mkPen('r'))
        for x in logs:
            lr = pqg.LinearRegionItem([x[0]-smjd,x[1]-smjd])
            self.plot1.addItem(lr)

        d1 = sqd.getdata('sr2_lock_atoms_gnd',self.fmjd,self.tmjd)
        d2 = sqd.getdata('sr2_lock_atoms_exc',self.fmjd,self.tmjd)
        d1 = sqd.rmerr(d1,'sr2')
        d2 = sqd.rmerr(d2,'sr2')
        logs = sqd.get_err_logs(self.fmjd,self.tmjd,['sr2'])
        self.plot2.plot(d1.mjd_tab()-smjd, 
                d1.val_tab(), pen=pqg.mkPen('g'))
        self.plot2.plot(d2.mjd_tab()-smjd, 
                d2.val_tab(), pen=pqg.mkPen('r'))
        for x in logs:
            lr = pqg.LinearRegionItem([x[0]-smjd,x[1]-smjd])
            self.plot2.addItem(lr)

        d1 = ftl.Sr2vsSr1_range(self.fmjd, self.tmjd)
        d1 = sqd.rmerr(d1,['sr2','sr1'])
        d1.plot_pqg_widget(widget=self.plot3,minusmjd=smjd)
        logs = sqd.get_err_logs(self.fmjd,self.tmjd,['sr2','sr1'])
        for x in logs:
            lr = pqg.LinearRegionItem([x[0]-smjd,x[1]-smjd])
            self.plot3.addItem(lr)

        d1 = ftl.fabsSr1_range_2(self.fmjd, self.tmjd)
        d1.plot_pqg_widget(widget=self.plot4, minusmjd=smjd)

        d1 = sqd.getdata('sr1_lock_prob_left',self.fmjd,self.tmjd)
        d2 = sqd.getdata('sr1_lock_prob_right',self.fmjd,self.tmjd)
        self.plot5.plot(d1.mjd_tab()-smjd, 
                d1.val_tab(), pen=pqg.mkPen('g'))
        self.plot5.plot(d2.mjd_tab()-smjd, 
                d2.val_tab(), pen=pqg.mkPen('r'))
        logs = sqd.get_err_logs(self.fmjd,self.tmjd,['sr1'])
        for x in logs:
            lr = pqg.LinearRegionItem([x[0]-smjd,x[1]-smjd])
            self.plot5.addItem(lr)

        d1 = sqd.getdata('sr2_lock_prob_left',self.fmjd,self.tmjd)
        d2 = sqd.getdata('sr2_lock_prob_right',self.fmjd,self.tmjd)
        self.plot6.plot(d1.mjd_tab()-smjd, 
                d1.val_tab(), pen=pqg.mkPen('g'))
        self.plot6.plot(d2.mjd_tab()-smjd, 
                d2.val_tab(), pen=pqg.mkPen('r'))
        logs = sqd.get_err_logs(self.fmjd,self.tmjd,['sr2'])
        for x in logs:
            lr = pqg.LinearRegionItem([x[0]-smjd,x[1]-smjd])
            self.plot6.addItem(lr)

