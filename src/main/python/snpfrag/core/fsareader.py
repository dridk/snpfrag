from Bio import SeqIO
from Bio.SeqIO import AbiIO
from scipy.signal import find_peaks
from scipy.stats import linregress
import numpy as np
import os
import logging
class FsaReader(object):
    def __init__(self, filename = None):
        self.regression = None 
        self.seqio = None
        self.set_filename(filename)

    def set_filename(self, filename: str):
        """set filename : it must be in abif firmat
        
        Arguments:
            filename {str} -- path to abif file
        """
        if not filename:
            return 
         
        if os.path.exists(filename):
            self.filename = filename 
            self.seqio = SeqIO.read(self.filename, 'abi')
        else:
            logging.critical("no filename")

    @property
    def description(self) -> str:
        """Return description 
        
        Returns:
            str -- 
        """
        return self.seqio.annotations["description"]

    @property
    def name(self) -> str:
        """Return name 
        
        Returns:
            str -- 
        """
        return self.seqio.annotations["name"]

    @property
    def dye_names(self) -> list:
        """Return list of dyes
        
        Returns:
            list -- list of dye name
        """
        names = []
        raw = self.seqio.annotations["abif_raw"]
        for key in raw.keys():
            if "DyeN" in key:
                names.append(raw[key])
        return names

    @property
    def channel_names(self) -> list:
        """Return list of channel names
        
        Returns:
            list -- list of channel names
        """
        names = []
        raw = self.seqio.annotations["abif_raw"]
        for key in raw.keys():
            if "DATA" in key:
                names.append(key)
        return names

    def channel_from_dye_name(self, name) -> str:
        """ Return Channel name from Dye Name """
        index = self.dye_names.index(name)
        return self.channel_names[index]

    def data(self, channel: str):
        """return Data according channel 
        
        Arguments:
            channel {str} -- channel name
        
        Returns:
            data -- list of number
        """
        return self.seqio.annotations["abif_raw"][channel]

    def normalize_data(self, channel):
        """ Return Data normalized from regression """ 
        if not self.regression:
            logging.warning("No regression has been set")
            return 
        
        data = self.data(channel)
        for index, value in enumerate(data):
            yield (self.transform(index),value)


    def compute_regression(self, channel : str, scales: list, height=80, distance=20):
        """Compute regression
        
        Arguments:
            channel {str} -- Size marker chanel name
            scales  {list} -- List of markers size. [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500]
        
    
        """
        # Get peaks data from channel 
        data = self.data(channel)

        # Detect peaks 
        peaks = find_peaks(data, height=height, distance=distance)[0]

        # Keep only last peaks according length of scales 
        # Because first peaks are probabily noise 
        l_scales = len(scales)
        l_peaks = len(peaks)
        diff = l_peaks - l_scales 
        peaks = peaks[diff:]

        # Do a linear regression 
        self.regression = linregress(peaks,scales)



    def transform(self, value):
        """Transform value according linear regression parameters
        Args:
            value (float or np.array)
        
        Returns:
            float: transformed value 
        """
        return value * self.regression.slope + self.regression.intercept

    def revert_transform(self, value):
        """Revert transformation accoding linear regression parameters 
        
        Arguments:
            value {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        return int((value - self.regression.intercept) / self.regression.slope)



if __name__ == "__main__":
    from PySide2.QtWidgets import *
    from PySide2.QtCore import *
    from PySide2.QtGui import *
    from PySide2.QtCharts import QtCharts as qtc
    import sys

    app = QApplication(sys.argv)

    view = qtc.QChartView()
    ROX_SIZE = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]

    reader = FsaReader("examples/A01.fsa")
    reader.compute_regression("DATA4",ROX_SIZE )

    chart = qtc.QChart()
    serie  = qtc.QLineSeries()

    for x, y in reader.normalize_data("DATA4"):
        serie.append(x,y)

    chart.addSeries(serie)
    view.setChart(chart)

    view.show()


    app.exec_()