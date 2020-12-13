import vcf 
import config
import iupac

from Bio import SeqIO
from Bio.SeqIO import AbiIO
from scipy.signal import find_peaks
from scipy.stats import linregress
import numpy as np


from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtCharts import QtCharts as qtc


class VCFReader(object):
    """Parse VCF to extract sample id 
    
    Args:
        filename (str): path to vcf file . It must be a bgzip indexed by tabix 
    """
    def __init__(self, filename = None):
        self.filename = filename 

    def compute(self):
        """Compute sample_id by reading vcf file
        """
        if self.filename is None:
            return
        device = open(self.filename,"rb")
        snps = []
        self.genotypes = []
        for snp in config.SNPS:
            device.seek(0)
            reader = vcf.Reader(device)

            variant = list(reader.fetch(snp["chr"], snp["pos"]-1, snp["pos"]))
            if variant:
                variant = variant[0]
                sample = variant.samples[0].sample
                gt = str(variant.genotype(sample).gt_bases).replace("/","")
                igt = iupac.genotype_to_iupac(gt)
                snps.append(igt)
                self.genotypes.append(gt)
            else:
                snps.append(snp["ref"])
                self.genotypes.append(snp["ref"]*2)


        self.sample_id = "".join(snps)




class FSAReader(object):
    def __init__(self, filename = str):
        self.filename = filename

    def compute(self):
        self.record = SeqIO.read(self.filename, 'abi')
        self.rox_data = self.record.annotations["abif_raw"][config.FRAG_ROX_CHANNEL]
        self.wt_data = self.record.annotations["abif_raw"][config.FRAG_MUT_CHANNEL]
        self.mt_data = self.record.annotations["abif_raw"][config.FRAG_WT_CHANNEL]
        self.calibrate_rox(self.rox_data)

        snps = []
        self.genotypes = []
        for snp in config.SNPS:
            fmin = self.unscale(snp["frag_min"])
            fmax = self.unscale(snp["frag_max"])

            peak_wt = find_peaks(self.wt_data[fmin:fmax], height=200, distance=20)[0]
            peak_mt = find_peaks(self.mt_data[fmin:fmax], height=200, distance=20)[0]

            genotype = str()
            # if peak_wt: 
            #     genotype += snp["ref"]
            
            # if peak_mt:
            #     genotype += snp["alt"]

            # if len(genotype) == 0:
            #     genotype == "NN"

            # if len(genotype) == 1:
            #     genotype = genotype * 2

            snps.append(iupac.genotype_to_iupac(genotype))
            self.genotypes.append(genotype)

        self.sample_id = "".join(snps)




    def calibrate_rox(self, data):
        """Detect Rox peaks and find parameter of the linear regression 
        """
        peaks = find_peaks(data, height=80, distance=20)[0]

        scales = config.FRAG_ROX_SIZE

        l_scales = len(scales)
        l_peaks = len(peaks)
        diff = l_peaks - l_scales 
        peaks = peaks[diff:]
        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = linregress(peaks,scales)



    def scale(self, value):
        """Transform value according linear regression previously computed
        
        Args:
            value (float or np.array)
        
        Returns:
            float: transformed value 
        """
        return value * self.slope + self.intercept

    def unscale(self, value):
        return int((value - self.intercept) / self.slope)


if __name__ == "__main__":
    
    reader = VCFReader("data/17-36166.vcf.gz")
    reader.compute()
    print(reader.sample_id)

    import sys 
    app = QApplication(sys.argv)

    reader = FSAReader("data/17-36166.fsa")
    reader.compute()

    print(reader.sample_id)

    view = qtc.QChartView()

    chart = qtc.QChart()
    serie  = qtc.QLineSeries()


    for x, y in enumerate(reader.wt_data[4211:4248]):
        serie.append(reader.scale(x),y)

    
    for i in serie.points():
        print(i)
        break
    

    chart.addSeries(serie)
    chart.createDefaultAxes()
    view.setRubberBand(qtc.QChartView.RectangleRubberBand)
    #chart.axisX().setRange(0,500)
    chart.axisY().setRange(0,1000)
    view.setChart(chart)

    view.show()
    
    app.exec_()