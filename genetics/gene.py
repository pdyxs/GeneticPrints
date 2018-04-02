from .codons import CODONS, codonAt

class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def initialiseCodons(self, nucleotides):
        index = self.start
        self.codons = []
        while index < self.end:
            self.codons += [codonAt(nucleotides, index)]
            index += 3

class Gene:
    'A Gene'

    def __repr__(self):
        return (self.name + ' [' + str(self.transcriptionRegion[0]) +
            ' (' + str(self.codingRegion[0]) + ', ' +
            str(self.codingRegion[1]) + ') ' + str(self.transcriptionRegion[1]) + ']');

    def __init__(self, data):
        (index, self.name, self.chromosome, self.strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds) = data.split("\t")
        self.index = int(index)
        self.transcriptionRegion = (int(txStart), int(txEnd))
        self.codingRegion = (int(cdsStart), int(cdsEnd))
        exonStarts = exonStarts.split(",")
        exonEnds = exonEnds.split(",")
        self.exons = []
        for i in range(int(exonCount)):
            self.exons.append(Exon(int(exonStarts[i]), int(exonEnds[i])))

    def contains(self, index):
        return index >= self.transcriptionRegion[0] and index < self.transcriptionRegion[1]

    def isCodingAt(self, index):
        return index >= self.codingRegion[0] and index < self.codingRegion[1]

    def exonAt(self, index):
        if not self.isCodingAt(index):
            return None
        for exon in self.exons:
            if index >= exon.start and index < exon.end:
                return exon
        return None

    def codonAt(self, index):
        exon = self.exonAt(index)
        if exon is not None:
            return exon.codons[int((index - exon.start)/3)]

    def initialiseSequenceData(self, nucleotides):
        for exon in self.exons:
            exon.initialiseCodons(nucleotides)

    def isExon(self, index):
        for exon in self.exons:
            if index < exon.start:
                return False
            if index < exon.end:
                return True
        return False

    def gaps(self):
        gaps = 0
        for i in range(len(self.exons)):
            if i > 0:
                gaps += self.exons[i].start - self.exons[i-1].end
        return gaps;

    def width(self):
        return (self.codingRegion[0] - self.transcriptionRegion[0] +
            int((self.codingRegion[1] - self.codingRegion[0] - self.gaps()) / 3) +
            self.transcriptionRegion[1] - self.codingRegion[1])
