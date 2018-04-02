from .nucleotides import NUCLEOTIDES

class Sequence:
    'Representing a DNA sequence'

    def __init__(self, nucleotideString, genes):
        self.loadNucleotides(nucleotideString)
        self.genes = []
        for i, gene in enumerate(genes):
            if (i > 0) and (gene.transcriptionRegion[0] == genes[i-1].transcriptionRegion[0]):
                continue
            gene.initialiseSequenceData(self.nucleotides)
            self.genes.append(gene)

    def loadNucleotides(self, nucleotideString):
        self.nucleotides = []
        for char in nucleotideString:
            self.nucleotides.append(NUCLEOTIDES[char])

    def length(self):
        return len(self.nucleotides)

    def characterAt(self, index):
        return self.nucleotides[index].character

    def genesAt(self, index):
        ret = []
        for gene in self.genes:
            if gene.contains(index):
                ret += [gene]
        return ret

    def codonAt(self, index):
        genes = self.genesAt(index)
        return genes[0].codonAt(index)
