

GENOTYPES = {

"AA": "A",
"CC": "C",
"GG": "G", 
"TT": "T",
"AG": "R",
"CT": "Y",
"CG": "S",
"AT": "W",
"GT":"K",
"AC":"M",
"N":"N"
}

def genotype_to_iupac(genotype:str):
    genotype = "".join(sorted(genotype))

    if genotype in GENOTYPES.keys():
        return GENOTYPES[genotype]
    else:
        return "N"