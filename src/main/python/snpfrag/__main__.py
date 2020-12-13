
import snpfrag
from snpfrag.gui.mainwindow import MainWindow
from snpfrag.gui.ficon import setFontPath 
from PySide2.QtWidgets import QApplication, QMainWindow
import sys
import os 
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    path = os.path.dirname(snpfrag.__file__)
    setFontPath(path+"/assets/fonts/materialdesignicons-webfont.ttf")

    w = MainWindow()
    w.show()

    app.exec_()