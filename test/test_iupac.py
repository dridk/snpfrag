from snpfrag.core import iupac
def test_genotype():
    assert iupac.genotype_to_iupac("AA") == "A"
    assert iupac.genotype_to_iupac("CCCFA") == "N"