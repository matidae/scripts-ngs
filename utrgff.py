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
        if "\tCDS\t" in line:
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
                return (Gene(start=0, end=0), genes[0])
            elif n == int(len(genes)):
                return (genes[n - 1], Gene(start=contig.end, end=contig.end))
            else:
                return (genes[n - 2], genes[n])


def checkOverlap(gen, prevnextgene, utr, cutoff, contig):
    prevgene = prevnextgene[0]
    nextgene = prevnextgene[1]
#    print vars(prevnextgene[0])
#    print vars(prevnextgene[1])
#    print "--------------------------------------"
    if gen.strand == "+":
        label5 = "5UTR"
        label3 = "3UTR"
    if gen.strand == "-":
        label3 = "5UTR"
        label5 = "3UTR"

    if gen.start - utr > prevgene.end + utr:
        printline(gen, gen.start - utr - 1, label5, gen.strand)
    elif gen.start - utr <= 0 :
        utraux = 1
        printline(gen, utraux, label5, gen.strand)
    else:
        utraux = int(math.ceil((prevgene.end - gen.start) / 2))
        if utraux > cutoff:
            printline(gen, gen.start - utraux - 1, label5, gen.strand)
        else:
            pass
    if gen.end + utr < nextgene.start - utr:
        printline(gen, gen.end + utr + 1, label3, gen.strand)
    elif gen.end + utr >= contig.end :
        utraux = contig.end
        printline(gen, utraux, label3, gen.strand)
    else:
        utraux = int(math.floor((nextgene.start - gen.end) / 2))
        if utraux > cutoff:
            printline(gen, gen.end + utraux, label3, gen.strand)
        else:
            pass


def getUTRs(gen, utr, genes, contig, cutoff):
    if contig.name == gen.contig and gen.strand == "+":
        genescontig = [g for g in genes if g.contig == contig.name ]
        prevnextgene = getPrevNextGene(gen, genescontig, contig)
        checkOverlap(gen, prevnextgene, utr, cutoff, contig)
    if contig.name == gen.contig and gen.strand == "-":
        genescontig = [g for g in genes if g.contig == contig.name ]
        prevnextgene = getPrevNextGene(gen, genescontig, contig)
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
    genespos = sorted([g for g in genes if g.strand == "+"], key = lambda x : x.start)
    genesneg = sorted([g for g in genes if g.strand == "-"], key = lambda x : x.start)

    for g in genes:
        print "\t".join(map(str, [g.contig, g.source, g.key, g.start, g.end, ".", g.strand, ".", "ID=" + g.name]))
        cont = [c for c in contigs if c.name == g.contig][0]
        if g.strand == "+":
            getUTRs(g, utr, genespos, cont, cutoff)
        else:
            getUTRs(g, utr, genesneg, cont, cutoff)
