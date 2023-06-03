import sys
sys.path.insert(1, '../mytools/')
import sqldata
import tools as tls
from PyQt5.QtWidgets import ( QGroupBox, QVBoxLayout, QHBoxLayout,
                QLineEdit, QPushButton,
                QTableWidget, QTableWidgetItem )
from PyQt5.QtCore import QTimer
from mQt_tools import QmjdEdit   


class logsBox(QGroupBox):
    def __init__(self, parent=None):
        super(logsBox, self).__init__(parent)
        
        self.rows = 10
        self.initUI()

    def initUI(self):
        self.blay = QVBoxLayout()
        self.blay.addLayout(self.barUI() )
        self.setLayout(self.blay)
        self.table = QTableWidget()
        self.blay.addWidget(self.table)
        self.table.setColumnCount(7)
        self.table.setRowCount(self.rows)
        #self.db2table()
        #self.refr_timer = QTimer(self)
        #self.refr_timer.timeout.connect(self.db2table)
        #self.refr_timer.start(5*60*1e3)

    def barUI(self):
        self.barlay = QHBoxLayout()
        self.fmjdEdit = QmjdEdit()
        self.tmjdEdit = QmjdEdit()
        self.refBut = QPushButton('R')
        self.refBut.clicked.connect(self.db2table)
        self.addBut = QPushButton('A')
        self.addBut.clicked.connect(self.add2db)
        self.modBut = QPushButton('M')
        self.modBut.clicked.connect(self.moddb)
        self.barlay.addWidget(self.fmjdEdit)
        self.barlay.addWidget(self.tmjdEdit)
        self.barlay.addWidget(self.refBut)
        self.barlay.addWidget(self.addBut)
        self.barlay.addWidget(self.modBut)

        return self.barlay
    
    def db2table(self):
        fmjd = self.fmjdEdit.val()
        tmjd = self.tmjdEdit.val()
        res = sqldata.get_logs(fmjd,tmjd)
        self.rows = len(res)
        self.table.setRowCount(self.rows+1)
        for num, x in enumerate(res):
            for i in [0,1,2,3,4,5,6]:
                self.table.setItem(num,i,QTableWidgetItem(str(x[i])) )
        self.table.resizeRowsToContents()
        self.table.resizeColumnToContents(0)
        self.table.resizeColumnToContents(1)
        self.table.resizeColumnToContents(3)
        self.table.resizeColumnToContents(4)
        self.table.resizeColumnToContents(5)
        self.table.resizeColumnToContents(6)

    def moddb(self):
        row = self.table.currentRow()
        s = [self.table.item(row,i).text() for i in [0,1,2,3,4] ]
        sqldata.modifyMessage(*s)

    def add2db(self):
        row = self.table.currentRow()
        s = [self.table.item(row,i).text() for i in [0,1,2,3] ]
        sqldata.sendMessage(*s)
