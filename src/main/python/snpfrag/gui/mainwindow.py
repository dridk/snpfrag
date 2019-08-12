from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

from snpfrag.gui.fsawidget import FsaWidget

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        super().__init__(parent)

        self.w = FsaWidget()
        self.w.set_filename("/home/schutz/Dev/snpfrag/examples/A01.fsa")

        self.setCentralWidget(self.w)


        