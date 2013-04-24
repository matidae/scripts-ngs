#!/usr/bin/python
import sys
import math


class Gene:
    def __init__(self, **kwargs):
        vars(self).update(kwargs)


class Contig:
    def __init__(self, **kwargs):
        vars(self).update(kwargs)


def getContigs(gff_file):
    contigs = []
    for line in gff_file:
        if "##sequence-region" in line:
            contigs.append((Contig(
                name=line.split()[1], end=int(line.split()[3].rstrip()))))
        if "CDS" in line:
            break
    return contigs


def printHeader(gff_file, contigs):
    print "#gff-version\t3\n##feature-ontology\tso.obo\n##attribute-ontology\tgff3_attributes.obo"
    for c in contigs:
        print "\t".join(map(str, ["##sequence-region", c.name, "1", c.end]))


def getGenes(gff_file):
    genes = []
    for line in gff_file:
        if "CDS" in line:
            genes.append(Gene(contig=line.split()[0], source=line.split()[1],
                              key=line.split()[2], start=int(line.split()[3]),
                              end=int(line.split()[4]), strand=line.split()[6],
                              name=line.split()[8].split(";")[0][3:].rstrip()))
    return genes


def printline(gen, utr, label, strand):
    if label == "5UTR" and strand == "+" or label == "3UTR" and strand == "-":
        print "\t".join(map(str, [gen.contig, gen.source, label, utr, g.start - 1, ".", g.strand, ".", "ID=" + label + g.name]))
    if label == "3UTR" and strand == "+" or label == "5UTR" and strand == "-":
        print "\t".join(map(str, [gen.contig, gen.source, label, g.end + 1, utr, ".", g.strand, ".", "ID=" + label + g.name]))


def getPrevNextGene(gen, genes, contig):
    n = 0
    for g in genes:
        n += 1
        if g.name == gen.name:
            if n == 1:
                return (Gene(start=0, end=0), genes[1])
            elif n == int(len(genes)):
                return (genes[n - 1], Gene(start=contig.end, end=contig.end))
            else:
                return (genes[n - 2], genes[n])


def checkOverlap(gen, prevnextgene, utr, cutoff, contig):
    prevgene = prevnextgene[0]
    nextgene = prevnextgene[1]
    if gen.start > prevgene.end:
        if gen.strand == "+":
            label = "3UTR"
        else:
            label = "5UTR"
        if gen.start - utr > prevgene.end + utr:
            printline(gen, gen.start - utr, label, gen.strand)
        else:
            if prevgene.end == 0:
                utraux = gen.start - 1
            else:
                utraux = int(math.ceil(((gen.start - prevgene.end) * .5)))
            if utraux > cutoff:
                printline(gen, gen.start - utraux, label, gen.strand)
    if not gen.end >= nextgene.start and gen.end <= nextgene.end:
        if gen.strand == "+":
            label = "3UTR"
        else: 
            label = "5UTR"
        if gen.end + utr < nextgene.start - utr:
            printline(gen, gen.end + utr, label, gen.strand)
        else:
            if nextgene.start == contig.end:
                utraux = contig.end - gen.end
            else:
                utraux = int(((nextgene.start - gen.end) * .5) - 1)
            if utraux > cutoff:
                printline(gen, gen.end + utraux, label, gen.strand)


def getUTRs(gen, utr, genes, contig, cutoff):
    if contig.name == gen.contig and gen.strand == "+":
        prevnextgene = getPrevNextGene(gen, genes, contig)
        checkOverlap(gen, prevnextgene, utr, cutoff, contig)
    if contig.name == gen.contig and gen.strand == "-":
        prevnextgene = getPrevNextGene(gen, genes, contig)
        checkOverlap(gen, prevnextgene, utr, cutoff, contig)


if __name__ == "__main__":
    gff_file = open(sys.argv[1], 'r').readlines()
    utr = 200
    cutoff = 20
    if len(sys.argv) > 2: 
        utr = int(sys.argv[2])
    if len(sys.argv) > 3:
        cutoff = int(sys.argv[3])
    contigs = getContigs(gff_file)
    printHeader(gff_file, contigs)
    genes = getGenes(gff_file)
    genespos = [g for g in genes if g.strand == "+"]
    genesneg = [g for g in genes if g.strand == "-"]
    for g in genes:
        print "\t".join(map(str, [g.contig, g.source, g.key, g.start, g.end, ".", g.strand, ".", "ID=" + g.name]))
        cont = [c for c in contigs if c.name == g.contig][0]
        if g.strand == "+":
            getUTRs(g, utr, genespos, cont, cutoff)
        else:
            getUTRs(g, utr, genesneg, cont, cutoff)
