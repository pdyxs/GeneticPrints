import os
import math
from PIL import Image, ImageDraw

GENEFILE = "ncbiGene.txt"
FILENAME = "KM034562v1.fa"
VERSION = 3
OUTPUT_FOLDER = "out/"

class Nucleotide:
    def __repr__(self):
        return self.character

    character = ""
    connectsTo = ""
    colour = "rgb(255,255,255)"

    def __init__(self, character, connectsTo, colour):
        self.character = character
        self.connectsTo = connectsTo
        self.colour = colour

NUCLEOTIDES = {
    'A': Nucleotide('A', 'T', "rgb(200,0,0)"),
    'T': Nucleotide('T', 'A', "rgb(200,100,0)"),
    'G': Nucleotide('G', 'C', "rgb(0,0,200)"),
    'C': Nucleotide('C', 'G', "rgb(0,100,200)")
}

class Sequence:
    'Representing a DNA sequence'

    def __init__(self, nucleotideString, genes):
        self.loadNucleotides(nucleotideString)
        self.genes = []
        for i, gene in enumerate(genes):
            if (i > 0) and (gene.transcriptionRegion[0] == genes[i-1].transcriptionRegion[0]):
                continue
            self.genes.append(gene)

    def loadNucleotides(self, nucleotideString):
        self.nucleotides = []
        for char in nucleotideString:
            self.nucleotides.append(NUCLEOTIDES[char])

    def length(self):
        return len(self.nucleotides)

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
            self.exons.append((int(exonStarts[i]), int(exonEnds[i])))

    def isExon(self, index):
        for (s, e) in self.exons:
            if index < s:
                return False
            if index < e:
                return True
        return False

    def gaps(self):
        gaps = 0
        for i in range(len(self.exons)):
            if i > 0:
                gaps += self.exons[i][0] - self.exons[i-1][1]
        return gaps;

    def width(self):
        return (self.codingRegion[0] - self.transcriptionRegion[0] +
            int((self.codingRegion[1] - self.codingRegion[0] - self.gaps()) / 3) +
            self.transcriptionRegion[1] - self.codingRegion[1])

def readFileLines(filename, startObject, lineFunc):
    f = open(filename, "r")
    for line in f:
        startObject = lineFunc(line, startObject)
    return startObject

pairString = readFileLines(FILENAME, "",
    lambda l, o: o if l.startswith(">") else o + l.strip());

def buildGene(line, o):
    if len(line) > 1:
        o.append(Gene(line))
    return o
genes = readFileLines(GENEFILE, [],
    lambda l, o: o if len(l) <= 1 else o + [Gene(l)])

sequence = Sequence(pairString, genes)

def getWidth(sequence):
    extents = [0,0]
    pos = 0
    for count, gene in enumerate(sequence.genes):
        width = gene.width()
        overlap = 0
        if count > 0:
            overlap = max(0, sequence.genes[count - 1].transcriptionRegion[1] - gene.transcriptionRegion[0])
            overlap = (overlap + 1) % 2
        if count % 2 == 0:
            pos += overlap
            pos = pos + width
            extents[1] = max(extents[1], pos)
        else:
            pos = pos - width
            extents[0] = min(extents[0], pos)
    return (- extents[0], (extents[1] - extents[0]) + 2)

def getHeight(seq):
    pos = 0
    height = 0
    for i, gene in enumerate(seq.genes):
        if pos < gene.transcriptionRegion[0]:
            height = height + gene.transcriptionRegion[0] - pos
        else:
            height = height + pos - gene.transcriptionRegion[0]
        pos = gene.transcriptionRegion[1]
    height += seq.length() - seq.genes[-1].transcriptionRegion[1]
    return height + 3

(startPos, IMAGE_WIDTH) = getWidth(sequence)
IMAGE_SIZE = (IMAGE_WIDTH, getHeight(sequence))
print(IMAGE_SIZE)
print(sequence.nucleotides[-1])

im = Image.new('RGBA', IMAGE_SIZE)
draw = ImageDraw.Draw(im)

def drawSequence(draw,seq,startPos):
    pos = (startPos, 0)
    spos = 0
    for i, gene in enumerate(seq.genes):
        while spos < gene.transcriptionRegion[0]:
            draw.point([pos], seq.nucleotides[spos].colour)
            spos += 1
            pos = (pos[0], pos[1] + 1)
        while spos < gene.codingRegion[0]:
            draw.point([pos], seq.nucleotides[spos].colour)
            spos += 1
            pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
        while spos < gene.codingRegion[1]:
            if gene.isExon(spos):
                draw.point([(pos[0],pos[1]-1)], seq.nucleotides[spos].colour)
                draw.point([(pos[0],pos[1])], seq.nucleotides[spos+1].colour)
                draw.point([(pos[0],pos[1]+1)], seq.nucleotides[spos+2].colour)
                spos += 3
                pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
            else:
                draw.point([pos], seq.nucleotides[spos].colour)
                spos += 1
                pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
        while spos < gene.transcriptionRegion[1]:
            draw.point([pos], seq.nucleotides[spos].colour)
            spos += 1
            if i < len(seq.genes) - 1 and spos >= seq.genes[i+1].transcriptionRegion[0]:
                mid = (seq.genes[i+1].transcriptionRegion[0] + gene.transcriptionRegion[1]) / 2
                if spos < mid:
                    pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1] + 1)
                else:
                    pos = (pos[0] + (1 if i % 2 == 1 else -1), pos[1] + 1)
            else:
                pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
    while spos < seq.length():
        draw.point([pos], seq.nucleotides[spos].colour)
        spos += 1
        pos = (pos[0], pos[1] + 1)
    print(pos)

drawSequence(draw, sequence, startPos)

print(sequence.genes)
print(sequence.length())

if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

im.save(OUTPUT_FOLDER + FILENAME + '.' + str(VERSION) + '.png')
