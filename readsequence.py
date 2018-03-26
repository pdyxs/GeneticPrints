import os
import math
from PIL import Image, ImageDraw

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
    'A': Nucleotide('A', 'T', "rgb(100,0,0)"),
    'T': Nucleotide('T', 'A', "rgb(200,0,0)"),
    'G': Nucleotide('G', 'C', "rgb(0,0,100)"),
    'C': Nucleotide('C', 'G', "rgb(0,0,200)")
}

class Sequence:
    'Representing a DNA sequence'

    def __init__(self, nucleotideString, genes):
        self.loadNucleotides(nucleotideString)
        self.genes = genes

    def loadNucleotides(self, nucleotideString):
        self.nucleotides = []
        for char in nucleotideString:
            self.nucleotides.append(NUCLEOTIDES[char])

    def length(self):
        return len(self.nucleotides)

    def draw(self, draw, width):
        for counter, nucleotide in enumerate(self.nucleotides):
            coords = (counter % width, math.floor(counter / width))
            draw.point([coords], nucleotide.colour)

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

GENEFILE = "ncbiGene.txt"
FILENAME = "KM034562v1.fa"
VERSION = 3
OUTPUT_FOLDER = "out/"

pairString = readFileLines(FILENAME, "",
    lambda l, o: o if l.startswith(">") else o + l.strip());

def buildGene(line, o):
    if len(line) > 1:
        o.append(Gene(line))
    return o
genes = readFileLines(GENEFILE, [],
    lambda l, o: o if len(l) <= 1 else o + [Gene(l)])

sequence = Sequence(pairString, genes)

IMAGE_WIDTH = math.ceil(math.sqrt(sequence.length()))
IMAGE_SIZE = (IMAGE_WIDTH, IMAGE_WIDTH)

im = Image.new('RGBA', IMAGE_SIZE)
draw = ImageDraw.Draw(im)

if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

sequence.draw(draw, IMAGE_WIDTH)

im.save(OUTPUT_FOLDER + FILENAME + '.' + str(VERSION) + '.png')
