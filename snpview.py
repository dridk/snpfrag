from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *
import sys
import config

from reader import VCFReader, FSAReader

class SnpModel(QAbstractTableModel):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.headers = ["snpid","chr","pos","ref","alt","vcf","fsa", "Match"]
        self.snp_data = []


    def rowCount(self, parent = QModelIndex()):
        return len(self.snp_data)

    def columnCount(self, parent = QModelIndex()):
        return len(self.headers) 

    def headerData(self, section, orientation, role ):
        if orientation == Qt.Horizontal:
            if role == Qt.DisplayRole:
                return self.headers[section]
        return None

    def data(self, index, role ):
        if not index.isValid():
            return None 

        if role == Qt.DisplayRole:
            return self.snp_data[index.row()][index.column()]

        if role == Qt.TextColorRole and index.column() == self.columnCount() -1 :
            match = self.snp_data[index.row()][index.column()]

            if match:
                return QColor("green")
            else:
                return QColor("red")

            



        return None

    def load(self):
        self.beginResetModel()
        
        self.snp_data = []
        
        for id, snp in enumerate(config.SNPS):

            vcf_result = "?"
            fsa_result = "?"

            if self.fsa_reader:
                fsa_result = self.fsa_reader.genotypes[id]

            if self.vcf_reader:
                vcf_result = self.vcf_reader.genotypes[id]

            match = vcf_result == fsa_result

            self.snp_data.append((
                snp["rsid"],
                snp["chr"],
                snp["pos"],
                snp["ref"],
                snp["alt"],
                vcf_result,
                fsa_result,
                match
                 ))

            
        self.endResetModel()


    def set_vcf_reader(self, vcf_reader):
        self.vcf_reader = vcf_reader

    def set_fsa_reader(self,fsa_reader):
        self.fsa_reader = fsa_reader





class SnpView(QTableView):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.model = SnpModel()
        self.setModel(self.model)
        self.horizontalHeader().setStretchLastSection(True)




if __name__ == "__main__":
    
    app = QApplication(sys.argv)

    view = SnpView()
    view.show()

    vcf_reader = VCFReader("data/17-36166.vcf.gz")
    fsa_reader = FSAReader("data/17-36166.fsa")

    vcf_reader.compute()
    fsa_reader.compute()

    view.model.set_fsa_reader(fsa_reader)
    view.model.set_vcf_reader(vcf_reader)

    view.model.load()


    print(vcf_reader.genotypes)
    print(fsa_reader.genotypes)


    app.exec_()
