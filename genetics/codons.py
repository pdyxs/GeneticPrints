class Codon:
    def __repr__(self):
        return self.name

    def __init__(self, name, symbol, sequences):
        self.name = name
        self.symbol = symbol
        self.sequences = sequences.split(", ")

CODONS = [
    Codon("Isoleucine", "I", "ATT, ATC, ATA"),
    Codon("Leucine", "L", "CTT, CTC, CTA, CTG, TTA, TTG"),
    Codon("Valine", "V", "GTT, GTC, GTA, GTG"),
    Codon("Phenylalanine", "F", "TTT, TTC"),
    Codon("Methionine", "M", "ATG"),
    Codon("Cysteine", "C", "TGT, TGC"),
    Codon("Alanine", "A", "GCT, GCC, GCA, GCG"),
    Codon("Glycine", "G", "GGT, GGC, GGA, GGG"),
    Codon("Proline", "P", "CCT, CCC, CCA, CCG"),
    Codon("Threonine", "T", "ACT, ACC, ACA, ACG"),
    Codon("Serine", "S", "TCT, TCC, TCA, TCG, AGT, AGC"),
    Codon("Tyrosine", "Y", "TAT, TAC"),
    Codon("Tryptophan", "W", "TGG"),
    Codon("Glutamine", "Q", "CAA, CAG"),
    Codon("Asparagine", "N", "AAT, AAC"),
    Codon("Histidine", "H", "CAT, CAC"),
    Codon("Glutamic acid", "E", "GAA, GAG"),
    Codon("Aspartic acid", "D", "GAT, GAC"),
    Codon("Lysine", "K", "AAA, AAG"),
    Codon("Arginine", "R", "CGT, CGC, CGA, CGG, AGA, AGG"),
    Codon("Stop codons", "Stop", "TAA, TAG, TGA")
]

def codonAt(nucleotides, index):
    if (len(nucleotides) < index + 3):
        return None
    seq = nucleotides[index].character + nucleotides[index + 1].character + nucleotides[index + 2].character
    for codon in CODONS:
        if seq in codon.sequences:
            return codon
