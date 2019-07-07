
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

from fsawidget import * 
from glob import glob

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        super().__init__()


        self.fsa_widget = FsaWidget()
        self.list_view = QListWidget()

        self.setCentralWidget(self.fsa_widget)

        dock = QDockWidget()
        dock.setWidget(self.list_view)

        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

        self.set_fsa_directory("data")

        self.list_view.currentItemChanged.connect(self.fsa_changed)


    def set_fsa_directory(self, dir):

        self.list_view.clear()
        for file in glob(f"{dir}/*.fsa"):
            self.list_view.addItem(file)

    def fsa_changed(self):
        self.fsa_widget.set_file(self.list_view.currentItem().text())



if __name__ == "__main__":
    import sys 
    
    app = QApplication(sys.argv)

    w = MainWindow()
    w.show()

    app.exec_()
