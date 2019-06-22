from Bio import SeqIO
from Bio.SeqIO import AbiIO

from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCharts import QtCharts as qtc
import sys 
from scipy.signal import find_peaks
import numpy as np 

record = SeqIO.read('test.fsa', 'abi')

data = record.annotations["abif_raw"]["DATA4"]


scales = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]

peaks = find_peaks(data, height=80, distance=20)[0]

predicted = np.diff(np.array(scales)) 
predicted = predicted / np.sum(predicted)

errors = []
for i in range(0,len(peaks) - len(scales)):
    observed = np.diff(np.array(peaks[i:len(scales)+i]))
    observed = observed / np.sum(observed)
    error = np.sum(np.abs(predicted - observed))
    errors.append(error)
    

start = errors.index(min(errors))

peaks = peaks[len(peaks)-len(scales):]

if __name__ == "__main__":
    
    app = QApplication(sys.argv)

    chart = qtc.QChart()
    


    view  = qtc.QChartView(chart)

    view.resize(800,600)

    view.setRubberBand(qtc.QChartView.RectangleRubberBand)

    serie = qtc.QLineSeries()
    serie.setPen(QPen("red"))
    
    serie = qtc.QLineSeries()
    serie.setPen(QPen("red"))

    scatter = qtc.QScatterSeries()
    scatter.setMarkerSize(10)
    scatter.setMarkerShape(qtc.QScatterSeries.MarkerShapeRectangle)

    for x in peaks:
        scatter.append(x, 1000)



    for x, y in enumerate(data):
        serie.append(x,y)

    chart.addSeries(serie)
    chart.addSeries(scatter)

    chart.createDefaultAxes()





    view.show()


    app.exec_()