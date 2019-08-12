
import snpfrag
from snpfrag.gui.mainwindow import MainWindow
from snpfrag.gui.ficon import setFontPath 
from PySide2.QtWidgets import QApplication, QMainWindow
from fbs_runtime.application_context.PySide2 import ApplicationContext

import sys
import os 

if __name__ == '__main__':
    appctxt = ApplicationContext()       # 1. Instantiate ApplicationContext

    path = os.path.dirname(snpfrag.__file__)
    setFontPath(path+"/assets/fonts/materialdesignicons-webfont.ttf")


    window = MainWindow()
    window.resize(800, 600)
    window.show()
    exit_code = appctxt.app.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)