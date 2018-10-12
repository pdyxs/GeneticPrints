import os
import math
from PIL import Image, ImageDraw
from genetics.sequence import Sequence
from genetics.gene import Gene
from utils import filereader
from genetics import colours

GENEFILE = "ncbiGene.txt"
SEQUENCE_FILE = "KM034562v1.fa"
VERSION = 5
OUTPUT_FOLDER = "out/"

def readSequenceString(filename, commentStr):
    def addLine(line, full):
        return full if line.startswith(commentStr) else full + line.strip()
    return filereader.readLinesTo(filename, "", addLine)

pairString = readSequenceString(SEQUENCE_FILE, ">")

def readGenes(filename):
    def buildGene(line, list):
        if len(line) > 1:
            list.append(Gene(line))
        return list
    return filereader.readLinesTo(filename, [], buildGene)

genes = readGenes(GENEFILE)

sequence = Sequence(pairString, genes)

def getWidth(sequence):
    extents = [0,0]
    pos = 0
    for count, gene in enumerate(sequence.genes):
        width = gene.width() * 3
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

im = Image.new('RGBA', IMAGE_SIZE)
draw = ImageDraw.Draw(im)

def drawSequence(draw,seq,startPos):
    pos = (startPos, 0)
    spos = 0
    for i, gene in enumerate(seq.genes):
        while spos < gene.transcriptionRegion[0]:
            draw.point([pos], colours.inSequence(seq, spos))
            spos += 1
            pos = (pos[0], pos[1] + 1)
        while spos < gene.codingRegion[0]:
            draw.point([pos], colours.inSequence(seq, spos))
            spos += 1
            pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
        while spos < gene.codingRegion[1]:
            if gene.isExon(spos):
                dx = (1 if i % 2 == 0 else -1)
                draw.point([(pos[0],pos[1]-1)], colours.inSequence(seq, spos))
                draw.point([(pos[0],pos[1])], colours.inSequence(seq, spos + 1))
                draw.point([(pos[0],pos[1]+1)], colours.inSequence(seq, spos + 2))
                draw.point([(pos[0]+dx,pos[1]-1)], colours.inSequence(seq, spos))
                draw.point([(pos[0]+dx,pos[1])], colours.inSequence(seq, spos + 1))
                draw.point([(pos[0]+dx,pos[1]+1)], colours.inSequence(seq, spos + 2))
                draw.point([(pos[0]+2*dx,pos[1]-1)], colours.inSequence(seq, spos))
                draw.point([(pos[0]+2*dx,pos[1])], colours.inSequence(seq, spos + 1))
                draw.point([(pos[0]+2*dx,pos[1]+1)], colours.inSequence(seq, spos + 2))
                colour = colours.CODON_COLOURS[seq.codonAt(spos).symbol]
                # draw.point([(pos[0],pos[1]-2)], colour)
                # draw.point([(pos[0],pos[1]-3)], colour)
                # draw.point([(pos[0],pos[1]+2)], colour)
                # draw.point([(pos[0],pos[1]+3)], colour)
                # draw.point([(pos[0]+dx,pos[1]+2)], colour)
                # draw.point([(pos[0]+dx,pos[1]+3)], colour)
                # draw.point([(pos[0]+2*dx,pos[1]+2)], colour)
                # draw.point([(pos[0]+2*dx,pos[1]+3)], colour)
                # draw.point([(pos[0],pos[1]-1)], colour)
                draw.point([(pos[0],pos[1])], colour)
                draw.point([(pos[0],pos[1]+1)], colour)
                draw.point([(pos[0]+dx,pos[1]-1)], colour)
                # draw.point([(pos[0]+dx,pos[1])], colour)
                draw.point([(pos[0]+dx,pos[1]+1)], colour)
                draw.point([(pos[0]+2*dx,pos[1]-1)], colour)
                draw.point([(pos[0]+2*dx,pos[1])], colour)
                # draw.point([(pos[0]+2*dx,pos[1]+1)], colour)
                spos += 3
                pos = (pos[0] + 3 * dx, pos[1])
            else:
                draw.point([pos], colours.inSequence(seq, spos))
                spos += 1
                pos = (pos[0] + (1 if i % 2 == 0 else -1), pos[1])
        while spos < gene.transcriptionRegion[1]:
            draw.point([pos], colours.inSequence(seq, spos))
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
        draw.point([pos], colours.inSequence(seq, spos))
        spos += 1
        pos = (pos[0], pos[1] + 1)

drawSequence(draw, sequence, startPos)

if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

im.save(OUTPUT_FOLDER + SEQUENCE_FILE + '.' + str(VERSION) + '.png')
