from Bio import SeqIO
from Bio.SeqIO import AbiIO

from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from PySide2.QtCharts import QtCharts as qtc
import sys 
from scipy.signal import find_peaks
from scipy.stats import linregress
import numpy as np 

from glob import glob
from copy import copy

class PeaksModel(object):
    def __init__(self):
        self.scales = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]

    def find_peaks(self, data):
        peaks = find_peaks(data, height=80, distance=20)[0]
        l_scales = len(self.scales)
        l_peaks = len(peaks)
        diff = l_peaks - l_scales 
        self.peaks = peaks[diff:]
        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = linregress(self.peaks,self.scales)


    def predict(self, value: float):
        return value * self.slope + self.intercept


class FragmentWidget(qtc.QChartView):
    def __init__(self, filename, parent = None):
        super().__init__(parent)
        self.filename = filename
        self.record = SeqIO.read(self.filename, 'abi')
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMinimumHeight(400)
        self.build()

    def build(self):
        rox_data = self.record.annotations["abif_raw"]["DATA4"]
        wt_data = self.record.annotations["abif_raw"]["DATA1"]
        mt_data = self.record.annotations["abif_raw"]["DATA2"]
        self.model =  PeaksModel()
        self.model.find_peaks(rox_data)

        self.chart = qtc.QChart()
        self.chart.setTitle(self.filename)
        self.rox_series = qtc.QLineSeries()
        self.rox_series.setColor(Qt.red)

        self.wt_series = qtc.QLineSeries()
        self.wt_series.setColor(QColor("#4e9a06"))
    
        self.mt_series = qtc.QLineSeries()
        self.mt_series.setColor(QColor("#27A4DD"))
    
        for x, y in enumerate(rox_data):
            self.rox_series.append(self.model.predict(x),y)


        for x, y in enumerate(wt_data):
            self.wt_series.append(self.model.predict(x),y)

        for x, y in enumerate(mt_data):
            self.mt_series.append(self.model.predict(x),y)        


        self.chart.addSeries(self.rox_series)
        self.chart.addSeries(self.wt_series)
        self.chart.addSeries(self.mt_series)
        self.chart.createDefaultAxes()
        self.setRubberBand(qtc.QChartView.RectangleRubberBand)
        self.chart.axisX().setRange(285,290)
        self.chart.axisY().setRange(0,2000)

        self.setChart(self.chart)




if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    

    layout = QVBoxLayout()

    
    wlist = []
    for filename in glob("olddata/*.fsa"):
        w = FragmentWidget(filename)
        wlist.append(w)
        layout.addWidget(w)

    widget = QWidget()
    widget.setLayout(layout)





    view = QScrollArea()
    view.setWidgetResizable(True)
    view.setWidget(widget)
    # rox_series = qtc.QLineSeries()
    # rox_series.setPen(QPen("red"))
    
    # rox_series = qtc.QLineSeries()
    # rox_series.setPen(QPen("red"))

    # scatter = qtc.QScatterSeries()
    # scatter.setMarkerSize(10)
    # scatter.setMarkerShape(qtc.QScatterSeries.MarkerShapeRectangle)

    # for x in peaks:
    #     scatter.append(x, 1000)



    # for x, y in enumerate(rox_data):
    #     rox_series.append(x,y)

    # chart.addSeries(rox_series)
    # chart.addSeries(scatter)






    view.show()


    app.exec_()