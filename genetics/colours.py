from .nucleotides import NUCLEOTIDES

def inSequence(sequence, index):
    return NUCLEOTIDE_COLOURS[sequence.characterAt(index)]

NUCLEOTIDE_COLOURS = {
    'A': "rgb(200,0,0)",
    'T': "rgb(200,100,0)",
    'G': "rgb(0,0,200)",
    'C': "rgb(0,100,200)"
}

CODON_COLOURS = {
    "I": "rgb(20,200,0)",
    "L": "rgb(40,200,0)",
    "V": "rgb(60,200,0)",
    "F": "rgb(80,200,0)",
    "M": "rgb(0,200,200)",
    "C": "rgb(100,200,0)",
    "A": "rgb(120,200,0)",
    "G": "rgb(140,200,0)",
    "P": "rgb(0,200,20)",
    "T": "rgb(0,200,40)",
    "S": "rgb(0,200,60)",
    "Y": "rgb(0,200,80)",
    "W": "rgb(0,200,100)",
    "Q": "rgb(0,200,120)",
    "N": "rgb(120,200,120)",
    "H": "rgb(20,200,20)",
    "E": "rgb(40,200,40)",
    "D": "rgb(60,200,60)",
    "K": "rgb(80,200,80)",
    "R": "rgb(100,200,100)",
    "Stop": "rgb(200,200,0)"
}
