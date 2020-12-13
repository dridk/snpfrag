from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from PySide2.QtCharts import QtCharts as qtc

from snpfrag.gui.ficon import FIcon
from snpfrag.core.fsareader import FsaReader
from snpfrag.config import DEFAULT_CHANNEL, DEFAULT_SCALES

import sys

class CalibrationWidget(QDialog):
    def __init__(self, parent = None):
        super().__init__()

        self.reader = FsaReader()
        self.treewidget = QTreeWidget()
        self.channel_view = qtc.QChartView()
        self.regression_view = qtc.QChartView()

        # setup tree 
        self.treewidget.setColumnCount(2)
        self.treewidget.setHeaderLabels(["Rox size", "predicted"])
        self.channel_view.setChart(qtc.QChart())
        self.regression_view.setChart(qtc.QChart())

        self.channel_serie = qtc.QLineSeries()
        self.regression_serie = qtc.QScatterSeries()        
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        self.channel_view.setChart(qtc.QChart())
        self.regression_view.setChart(qtc.QChart())

        v_splitter = QSplitter(Qt.Vertical)
        v_splitter.addWidget(self.regression_view)
        v_splitter.addWidget(self.channel_view)

        h_splitter = QSplitter(Qt.Horizontal)
        h_splitter.addWidget(self.treewidget)
        h_splitter.addWidget(v_splitter)

        v_layout = QVBoxLayout()
        v_layout.addWidget(h_splitter) 
        v_layout.addWidget(self.buttons)
        self.setLayout(v_layout)


    def set_filename(self, filename):
        self.reader.set_filename(filename)
        self.reader.compute_regression(DEFAULT_CHANNEL,DEFAULT_SCALES)
        self.plot_channel("DATA4")
        self.plot_regression("DATA4")
        self.load_size()

    def plot_regression(self, channel):
        self.regression_serie.clear()   
        self.regression_view.chart().removeAllSeries()
        self.regression_serie.setName("Size Marker")
        
        for i,x in enumerate(self.reader.scales):
            self.regression_serie.append(x, self.reader.transform(self.reader.peaks[i]))
        
        self.regression_view.chart().addSeries(self.regression_serie)
        self.regression_view.chart().createDefaultAxes()



    def plot_channel(self,channel):
        self.channel_serie.clear()   
        self.channel_view.chart().removeAllSeries()
        self.channel_serie.setName("Size Marker")
        for x, y in self.reader.normalize_data(channel):
                self.channel_serie.append(x,y)
        

        # Display peaks 
        self.peaks_serie = qtc.QScatterSeries()
        for i, peak in enumerate(self.reader.peaks):
            self.peaks_serie.append(self.reader.transform(peak), 1000)

    
        self.channel_view.chart().addSeries(self.channel_serie)
        self.channel_view.chart().addSeries(self.peaks_serie)

        self.channel_view.chart().createDefaultAxes()

    def load_size(self):
        self.treewidget.clear()
        for index, size in enumerate(self.reader.scales):
            item = QTreeWidgetItem()
            item.setText(0, str(size))
            predicted = self.reader.transform(self.reader.peaks[index])
            item.setText(1, str(round(predicted,2)))
            self.treewidget.addTopLevelItem(item)

        

        

if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = CalibrationWidget()
    w.show()

    app.exec_()
