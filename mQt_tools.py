from PyQt5.QtWidgets import QWidget
from tools import TSerie
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtGui import QDoubleValidator
from time_tools import getMJD

class dataWidget(QWidget):
    def __init__(self, name):
        super().__init__()
        self.name = name
    
    def initUI(self):
        self.setWindowTitle(self.name)
        self.setGeometry(10,10,100,100)
        

class QtTSerie(TSerie):
    def __init__(self, label='x', mjd=[], val=[]):
        super().__init__(label,mjd,val)
        self.widget=dataWidget(name=label)


class QmjdEdit(QLineEdit):
    def __init__(self, parent=None):
        super(QmjdEdit, self).__init__(parent)
        self.valid = QDoubleValidator()
        self.setValidator(self.valid)
        self.setText('%.3f'%(int(getMJD()) ) )

    def val(self):
        v = float(self.text())
        if abs(v)<1000:
            v = v + getMJD()
        return v
