from Bio import SeqIO
from Bio.SeqIO import AbiIO


#  Echelle des tailles pour la calibration du ROX
# Les 16 premiers pics du signal ROX correspond à ces tailles
SCALE = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500]


# Chargement du fichier
record = SeqIO.read("test.fsa", "abi")


# 3 signals à l'interieur
#  ROX = Correspond au pic pour faire la calibration
rox_data = record.annotations["abif_raw"]["DATA4"]

# WT = WILD TYPE: Correspond au pic vert
wt_data = record.annotations["abif_raw"]["DATA1"]

#  MT = MUTANT: correspond au pic bleu
mt_data = record.annotations["abif_raw"]["DATA2"]


#  OBJECTIF 1 :
#  Afficher les signaux

#  OBJECTIF 2 :
# DETECTER LA POSITION ( le temps ) DES PEAKS avec Find PEAKS

#  OBJECTIF 3 :
#  FAIRE LA CORRELATION LINEAIRE ENTRE LES PEAKS DE ROX_DATA ET SCALE

# OBJECTIF 4 :
#  PREDIRE LA TAILLE DES FRAGMENTS SUR LES SIGNAUX WT ET MT .

# OBJECTIF 5 :
# ECRIRE LES RESULTATS DANS UN DATAFRAME ET UN FICHIER

# PEAKS    TIME    TAILLE  AMPLITUDE
# PEAK 1    43sec    2234    324234
