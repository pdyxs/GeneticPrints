class Nucleotide:
    def __repr__(self):
        return self.character

    def __init__(self, character, connectsTo):
        self.character = character
        self.connectsTo = connectsTo

NUCLEOTIDES = {
    'A': Nucleotide('A', 'T'),
    'T': Nucleotide('T', 'A'),
    'G': Nucleotide('G', 'C'),
    'C': Nucleotide('C', 'G')
}
