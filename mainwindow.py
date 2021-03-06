
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

from fsawidget import * 
from glob import glob
import config
from reader import FSAReader, VCFReader
from snpview import SnpView

import qrcode

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        super().__init__()


        self.list_view = QListWidget()

        self.tab_widget = QTabWidget()
        self.snp_view = SnpView()

        dock = QDockWidget()
        dock.setWidget(self.list_view)
        self.view = qtc.QChartView()

        # Create grids 
        gird_layout = QGridLayout()
        self.grid_widget = QWidget()
        self.grid_views = []
        for id, snp in enumerate(config.SNPS):
            x = int(id / 3 )
            y = int(id % 3 )
            g_view = qtc.QChartView()
            self.grid_views.append(g_view)
            gird_layout.addWidget(g_view,x,y)

        self.grid_widget.setLayout(gird_layout)

            
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.list_view.currentItemChanged.connect(self.fsa_changed)

        self.tab_widget.addTab(self.view, "view")
        self.tab_widget.addTab(self.grid_widget, "snps")

        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.tab_widget)
        self.splitter.addWidget(self.snp_view)


        self.setCentralWidget(self.splitter)

        self.barcode_widget = QWidget()
        self.barcode_layout = QVBoxLayout()
        self.barcode_img = QLabel()
        self.barcode_img.setPixmap(QPixmap("qrcode.png"))
        self.barcode_img.setMinimumSize(64,64)
        self.barcode_edit = QLineEdit()
        self.barcode_edit.setAlignment(Qt.AlignCenter)


        self.barcode_layout.addWidget(self.barcode_img)
        self.barcode_layout.addWidget(self.barcode_edit)
        self.barcode_layout.addStretch()
        self.barcode_widget.setLayout(self.barcode_layout)

        font = QFont()
        font.setPixelSize(20)
        self.barcode_edit.setFont(font)

        barcode_dock = QDockWidget()
        barcode_dock.setWidget(self.barcode_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea,barcode_dock)



        self.resize(800,400)

        toolbar = self.addToolBar("menu")
        run_action = toolbar.addAction("Run", self.run_analyse)
        self.readers = {}

        self.set_directory("R10")
        self.run_analyse()



    def set_directory(self, dir):
        self.dir = dir 
        self.list_view.clear()
        for file in glob(f"{dir}/*.fsa"):
            info = QFileInfo(file)
            self.list_view.addItem(info.baseName())

    def fsa_changed(self):
        basename = self.list_view.currentItem().text()
        
        if self.readers:
            vcf_reader = self.readers[basename][0]
            fsa_reader = self.readers[basename][1]
            #vcf_reader.compute()
            fsa_reader.compute()
            self.load_chart(fsa_reader)
            self.snp_view.model.set_vcf_reader(vcf_reader)
            self.snp_view.model.set_fsa_reader(fsa_reader)
            self.snp_view.model.load()

        self.create_barcode(vcf_reader.sample_id)


    def create_barcode(self,ident):
        self.barcode_edit.setText(ident)
        img = qrcode.make(ident)
        img.save("qrcode.png")
        self.barcode_img.setPixmap(QPixmap("qrcode.png"))

    def run_analyse(self):
        self.readers = {}

        for index in range(0, self.list_view.count()):

            item = self.list_view.item(index)
            basename = item.text()
            vcf_file = self.dir + "/" + basename + ".vcf.gz"
            fsa_file = self.dir + "/" + basename + ".fsa"
            print("analyse ", basename)
            self.readers[basename] = [VCFReader(vcf_file), FSAReader(fsa_file)]
            item.setIcon(QIcon.fromTheme("system-run"))

    def create_series(self, reader, m, M):
        rox_series = qtc.QLineSeries()
        rox_series.setColor(Qt.red)

        _m = reader.unscale(m)
        _M = reader.unscale(M)

        for x, y in enumerate(reader.rox_data[_m:_M]):
            rox_series.append(reader.scale(x),y)
        yield rox_series

        wt_series = qtc.QLineSeries()
        wt_series.setColor(QColor("#4e9a06"))
        for x, y in enumerate(reader.wt_data[_m:_M]):
            wt_series.append(reader.scale(x),y)
        yield wt_series


        mt_series = qtc.QLineSeries()
        mt_series.setColor(QColor("#27A4DD"))
        for x, y in enumerate(reader.mt_data[_m:_M]):
            mt_series.append(reader.scale(x),y)
        yield mt_series


    def load_chart(self, reader):
    
        self.chart = qtc.QChart()
        self.chart.setTitle(self.list_view.currentItem().text())

        for serie in self.create_series(reader,0,500):   
            self.chart.addSeries(serie)

        self.chart.createDefaultAxes()
        self.view.setRubberBand(qtc.QChartView.RectangleRubberBand)
        self.view.setChart(self.chart)

        for id, sview in enumerate(self.grid_views):
            chart = qtc.QChart()
            sview.setChart(chart)

            snp = config.SNPS[id]
            sview.chart().setTitle(snp["rsid"])
            f_min = snp["frag_min"]
            f_max = snp["frag_max"]
            
            for i, serie in enumerate(self.create_series(reader,f_min,f_max)):
                if i == 0:
                    continue   
                sview.chart().addSeries(serie)
 
                sview.chart().createDefaultAxes()
                sview.chart().axisX().hide()
                #sview.chart().axisX().setRange(f_min, f_max)

if __name__ == "__main__":
    import sys 
    
    app = QApplication(sys.argv)

    w = MainWindow()
    w.show()

    app.exec_()
