from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

from snpfrag.gui.ficon import FIcon
from snpfrag.gui.fsawidget import FsaWidget

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        super().__init__(parent)

        self.w = FsaWidget()
        self.w.set_filename("/home/schutz/Dev/snpfrag/examples/A01.fsa")

        # setup toolbar 
        self.menuBar = QMenuBar()
        ## add file menu
        file_menu = QMenu("&File")
        file_menu.addAction(FIcon(0xf76f), "Open", self.on_open)

        self.menuBar.addMenu(file_menu)
        self.setMenuBar(self.menuBar)

        self.setCentralWidget(self.w)




    def on_open(self):

        filename = QFileDialog.getOpenFileName(self,"open fsa file","","Fragment file (*.fsa)")
        if filename:
            self.w.set_filename(filename[0])


        