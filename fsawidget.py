
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtCharts import QtCharts as qtc

from Bio import SeqIO
from Bio.SeqIO import AbiIO
from scipy.signal import find_peaks
from scipy.stats import linregress


class PeaksModel(object):
    """Detect peaks from abif and create a linear regression """
    
    def __init__(self):
        self.scales = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]

    def find_peaks(self, data):
        """detect peaks from data 
        
        Arguments:
            data {list} -- list of values 
        """
        peaks = find_peaks(data, height=80, distance=20)[0]
        l_scales = len(self.scales)
        l_peaks = len(peaks)
        diff = l_peaks - l_scales 
        self.peaks = peaks[diff:]
        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = linregress(self.peaks,self.scales)


    def predict(self, value: float):
        """Transform value according linear regression previously computed
        
        Args:
            value (float): a scalar value
        
        Returns:
            float: transformed value 
        """
        return value * self.slope + self.intercept


class FragmentWidget(qtc.QChartView):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        #self.setMinimumHeight(400)
        
    def set_file(self, filename):
        self.filename = filename
        self.record = SeqIO.read(self.filename, 'abi')

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
        self.chart.axisX().setRange(0,500)

        self.setChart(self.chart)

class FsaWidget(QTabWidget):
    def __init__(self, parent = None):
        super().__init__()

        self.viewer = FragmentWidget()
        self.addTab(self.viewer, "viewer")

        self.snps_range = [

            (280,295),
            (239,244),
            (131,133),
            (280,295),
            (186,193),
            (186,193),         
            (280,295),
            (186,193),
            (465,470),
        
        ]

        # construct grid 
        self.grid_layout = QGridLayout()
        self.grid_widget = QWidget()

        self.snp_viewers  = []
        n_cols = 3 
        for i in range(9):
            x = int(i / n_cols)
            y = int(i % n_cols)
            
            w = FragmentWidget()
            self.grid_layout.addWidget(w, x,y)

            self.snp_viewers.append(w)


        self.grid_layout.setSpacing(0)
        self.grid_layout.setContentsMargins(0,0,0,0)
        self.grid_widget.setLayout(self.grid_layout)
        self.grid_area = QScrollArea()
        self.grid_area.setWidgetResizable(True)
        self.grid_area.setWidget(self.grid_widget)
        self.addTab(self.grid_area, "snps")



    def set_file(self, filename):
        self.viewer.set_file(filename)

        for i, w in enumerate(self.snp_viewers):
            w.set_file(filename)
            the_range = self.snps_range[i]
            w.chart.axisX().setRange(*the_range)
            w.chart.axisX().hide()
            w.chart.axisY().hide()
            w.chart.legend().hide()
            w.chart.setMargins(QMargins(0,0,0,0))
            w.chart.setTitle(f"snp {i}")





if __name__ == "__main__":
    import sys 
    
    app = QApplication(sys.argv)

    w = FsaWidget()
    w.set_file("data/17-38213_F01.fsa")
    w.show()

    app.exec_()

