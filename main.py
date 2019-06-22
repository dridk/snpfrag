from Bio import SeqIO
from Bio.SeqIO import AbiIO

from PySide2.QtWidgets import *
from PySide2.QtCharts import QtCharts as qtc
import sys 

record = SeqIO.read('test.fsa', 'abi')

data = record.annotations["abif_raw"]["DATA1"]


if __name__ == "__main__":
    
    app = QApplication(sys.argv)

    chart = qtc.QChart()

    view  = qtc.QChartView(chart)

    serie = qtc.QLineSeries()

    for x, y in enumerate(data):
        serie.append(x,y)

    chart.addSeries(serie)

    view.show()


    app.exec_()