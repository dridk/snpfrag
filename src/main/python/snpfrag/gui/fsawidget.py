from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from PySide2.QtCharts import QtCharts as qtc

from snpfrag.gui.ficon import FIcon
from snpfrag.core.fsareader import FsaReader

class FsaWidget(QWidget):
    def __init__(self, parent = None):
        super().__init__()
        
        self.toolbar = QToolBar()
        self.view = qtc.QChartView()
        #self.view.setRenderHint(QPainter.Antialiasing)
        self.chart = qtc.QChart()
        self.reader = FsaReader()

        # build toolbar
        self.dye_menu = QMenu()
        self.dye_btn = QPushButton("Dye")
        self.dye_btn.setIcon(FIcon(0xf755))
        self.dye_btn.setFlat(True)
        self.dye_btn.setMenu(self.dye_menu)
        self.toolbar.addWidget(self.dye_btn)        

        self.toolbar.addAction(FIcon(0xf293),"rescale", self.chart.zoomReset)
        self.toolbar.addAction(FIcon(0xf6ec),"Zoom in", self.chart.zoomIn)
        self.toolbar.addAction(FIcon(0xf6eb),"Zoom out", self.chart.zoomOut)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.view)
        self.setLayout(main_layout)

        # configure view 
        self.view.setRubberBand(qtc.QChartView.RectangleRubberBand)

        

    def set_filename(self, filename):
        self.reader.set_filename(filename)

        ROX_SIZE = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]
        self.reader.compute_regression("DATA4", ROX_SIZE)
        self._load_dye_menu()
        self.plot()

    def plot(self):
        self.series = {}

        for dye in self.reader.dye_names:
            channel = self.reader.channel_from_dye_name(dye)
            serie = qtc.QLineSeries()        
            serie.setName(dye)
            for x, y in self.reader.normalize_data(channel):
                serie.append(x,y)
            self.series[dye] = serie 
            self.chart.addSeries(serie)

            pen = serie.pen()
            pen.setWidth(1)
            serie.setPen(pen)
        
        self.chart.createDefaultAxes()
        self.view.setChart(self.chart)

    def _load_dye_menu(self):
        """ read and fill combobox with dye names """ 
        self.dye_menu.clear()
        for dye in self.reader.dye_names:
            action = self.dye_menu.addAction(dye)
            action.setCheckable(True)
            action.setChecked(True)
            action.triggered.connect(self._on_dye_clicked)
            self.dye_menu.addAction(action)

    def _on_dye_clicked(self, checked):

        dye_name  = self.sender().text()
        self.series[dye_name].setVisible(checked)





if __name__ == "__main__": 
    import sys
    import os
    import snpfrag
    from snpfrag.gui.ficon import setFontPath
    app = QApplication(sys.argv)
    path = os.path.dirname(snpfrag.__file__)
    setFontPath(path+"/assets/fonts/materialdesignicons-webfont.ttf")

    w = FsaWidget()
    w.set_filename("examples/A01.fsa")
    w.show()

    app.exec_()