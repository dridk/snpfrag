
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *


class FragmentWidget(QTabWidget):
    def __init__(self, parent = None):
        super().__init__()



if __name__ == "__main__":
    import sys 
    
    app = QApplication(sys.argv)

    w = FragmentWidget()
    w.show()

    app.exec_()

